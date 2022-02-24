#ifndef COOLING_HPP_
#define COOLING_HPP_

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include <fstream>
#include <sstream>

class Cooling {
public:
  Cooling(const std::string in);
  ~Cooling();

  Real townsend_cooling(const Real temp, const Real rho, const Real dt);
  Real single_point_cooling_time(const Real T, const Real rho);

  // Load table from file
  void read_points(const std::string in, 
                   std::vector<Real> &temperatures,
                   std::vector<Real> &cooling, 
                   std::vector<Real> &heating);

  // Accessors
  Real Get_tceil() const {return Tceil_;}
  Real Get_tfloor() const {return Tfloor_;}

private:
  // Loaded Arrays
  AthenaArray<Real> temperature_table;
  AthenaArray<Real> cooling_table;
  AthenaArray<Real> heating_table;
  int nbins_;

  // Userful physical quantities 
  Real Tceil_;
  Real Tfloor_;
  Real const_factor_;

  // Scratch Arrays
  AthenaArray<Real> net_cooling;
  AthenaArray<Real> temps;
  AthenaArray<Real> Ys; 

  // Helper
  int get_temp_index(const Real T);
  Real inst_cooling(const Real T, const Real rho);
  Real get_zero_point(const int i, const Real rho);
  Real compute_constant_factor();
  bool isClose(const Real a, const Real b, const Real tol = 1e-20);
};

// Constructor
Cooling::Cooling(const std::string in) {
    // Read in table
    std::ifstream infile;
    infile.open(in.c_str());
    std::string line;
    std::vector<Real> temperatures, cooling, heating;
    while(getline(infile,line) && line[0] != '#'){
        std::stringstream ss(line);
        Real T,c,h;
        ss >> T >> c >> h;
        temperatures.push_back(T);
        cooling.push_back(c);
        heating.push_back(h);

    }
    infile.close();

    nbins_ = temperatures.size();
    // Assert that they have the same lengths
    if (nbins_ != cooling.size() || nbins_ != heating.size()) { 
        std::stringstream msg;
        msg << "### ERROR : Heating/Cooling Tables are different sizes!" << std::endl;
        ATHENA_ERROR(msg); 
    }

    AthenaArray<Real> temperature_table(nbins_);
    AthenaArray<Real> cooling_table(nbins_);
    AthenaArray<Real> heating_table(nbins_);

    for (int i=0; i<=nbins_; ++i) {
        temperature_table(i) = temperatures.at(i);
        cooling_table(i) = cooling.at(i);
        heating_table(i) = heating.at(i);
    }

    Tceil_ = temperature_table(-1);
    Tfloor_ = temperature_table(0);
    const_factor_ = compute_constant_factor();

    net_cooling.NewAthenaArray(nbins_);
    temps.NewAthenaArray(nbins_);
    Ys.NewAthenaArray(nbins_);

  return;
}

// Destructor
Cooling::~Cooling() {
    return;
}

// Exact Integration Scheme for Radiative Cooling from Townsend (2009)
// Expanded to include heating
// Returns: Temperature(K) at the next timestep after cooling/heating
// Requires: - Input temperature, density, and timestep in cgs units
Real Cooling::townsend_cooling(const Real T, const Real rho, const Real dt) {
    // Ensure temperature is within bounds
    if (T <= Tfloor_) return Tfloor_;
    if (T >= Tceil_) return Tceil_;

    // Make a copy of temperatures that can be edited
    temps = temperature_table;

    // Calculate net cooling for a given density
    for (int i=0; i<=nbins_; ++i) {
        net_cooling(i) = cooling_table(i)-heating_table(i)/rho;
    }

    // Get index of the temperature bin
    int T_idx = get_temp_index(T);

    // Check if there is net cooling or heating
    Real net_cool_T = inst_cooling(T, rho);

    // If no net cooling/heating, temperature does not change
    if (isClose(net_cool_T,0.)) return T;

    // Flags for piece-wise linear bins - only at start or end
    int start_linear_idx = -1;
    int end_linear_idx = -1;

    // Our aim is to calculate the new temperature
    Real T_new = 0.0;

    if (net_cool_T > 0) { // Cooling
        // We save time by using the next bin as the point of reference for the TEF
        // Hence we need to check for the edge case where it is a heating bin
        if (net_cooling(T_idx+1) < 0) {
            // If so we figure out what the equilibruim temperature is
            // Note everything after this index should not be used anymore
            // We then hijack the next bin to mark this temperature
            // And set the cooling there to zero
            // The piecewise power law assumption is still here, so
            //   cooling in this bin will be automatically suppressed
            temps(T_idx+1)       = get_zero_point(T_idx,rho);
            net_cooling(T_idx+1) = 1e-20;
            start_linear_idx     = T_idx;
        }

        // Set reference values
        Real T_ref = temps(T_idx+1);
        Real C_ref = net_cooling(T_idx+1);
        
        // Compute Y from idx = temp_bin+1 downwards
        // Calculate Y_k recursively (Eq. A6) from idx = temp_bin+1 downwards    
        Ys(T_idx+1) = 0.0;

        int i = T_idx;
        bool cooling_region_end = false;

        Real T_i, T_ip1, C_i, C_ip1;
        Real step, slope, sm1, oms;

        while (i >= 0 && !cooling_region_end) {
            if (net_cooling(i) < 0) { // if this bin has net heating
                temps(i) = get_zero_point(i,rho);          
                net_cooling(i) = 1e-20;
                end_linear_idx = i;
                cooling_region_end = true;
            }
                
            T_i = temps(i);
            C_i = net_cooling(i);

            T_ip1 = temps(i+1);
            C_ip1 = net_cooling(i+1);
            
            if (i == start_linear_idx || i == end_linear_idx) {
                slope = (T_ip1-T_i)/(C_ip1-C_i);
                step = (C_ref/T_ref)*slope*std::log(C_i/C_ip1);
            } else { 
                slope = std::log(C_ip1/C_i)/std::log(T_ip1/T_i);
                if (isClose(slope,1.0)) {
                    step = (C_ref/C_i)*(T_i/T_ref)*std::log(T_i/T_ip1);
                } else {
                    sm1  = slope - 1.0;
                    step = (C_ref/C_i)*(T_i/T_ref)*(std::pow(T_i/T_ip1,sm1)-1)/sm1;
                }
            }

            Ys(i) = Ys(i+1) - step;
            i -= 1;
        }

        Real T_min = temps(i+1);
        Real Y_lim = Ys(i+1);

        // Compute current TEF Y
        Real Y, Y_i;
        
        Y_i = Ys(T_idx);
        T_i = temps(T_idx);
        C_i = net_cooling(T_idx);
        T_ip1 = temps(T_idx+1);
        C_ip1 = net_cooling(T_idx+1);


        if (T_idx == start_linear_idx || T_idx == end_linear_idx) {
            slope = (T_ip1-T_i)/(C_ip1-C_i);
            Y = Y_i + (C_ref/T_ref)*slope*std::log(C_i/(C_i+(T-T_i)/slope));
        } else {
            slope = std::log(C_ip1/C_i)/std::log(T_ip1/T_i);
            if (isClose(slope,1.0)) {
                Y = Y_i + (C_ref/C_i)*(T_i/T_ref)*std::log(T_i/T);
            } else {
                sm1  = slope - 1.0;
                Y = Y_i + (C_ref/C_i)*(T_i/T_ref)*(std::pow(T_i/T,sm1)-1)/sm1;
            }
        }

        // Compute the new TEF after a timestep (Eqn. 26)
        Real Y_new = Y + rho*C_ref*const_factor_*dt/T_ref;

        // Check if we have reached the minimum temperature
        if (Y_new > Y_lim) { return T_min; }

        // TEF is a strictly decreasing function of T and new_tef > tef
        // Check if the new TEF falls into a lower bin
        // If so, update slopes and coefficients
        while ((T_idx > 0) && (Y_new > Ys(T_idx))) { T_idx -= 1; }
        
        T_i   = temps(T_idx);
        Y_i   = Ys(T_idx);
        C_i   = net_cooling(T_idx);
        T_ip1 = temps(T_idx+1);
        C_ip1 = net_cooling(T_idx+1);
        
        // Compute the Inverse Temporal Evolution Function Y^{-1}(Y) (Eq. A7)
        if (T_idx == start_linear_idx || T_idx == end_linear_idx) {
            slope = (C_ip1-C_i)/(T_ip1-T_i);
            T_new = T_i + (C_i/slope)*(1/std::exp((Y_new-Y_i)*(T_ref/C_ref)*slope)-1);
        } else {
            slope = std::log(C_ip1/C_i)/std::log(T_ip1/T_i);
            if (isClose(slope,1.0)) {
                T_new = T_i*std::exp(-(C_i/C_ref)*(T_ref/T_i)*(Y_new-Y_i));
            } else {
                oms  = 1.0 - slope;
                T_new = T_i*std::pow(1-oms*(C_i/C_ref)*(T_ref/T_i)*(Y_new-Y_i),1/oms);
            }
        }
    } else { // Heating, net_cool_T < 0
        if (net_cooling(T_idx+1) > 0) { // If the cuurent bin is cooling and not heating
            temps(T_idx)       = get_zero_point(T_idx,rho);
            net_cooling(T_idx) = -1e-20;
            start_linear_idx   = T_idx;
        }
        
        // Set reference values
        Real T_ref = temps(T_idx);
        Real C_ref = net_cooling(T_idx);

        // Compute Y from idx = temp_bin upwards
        Ys(T_idx) = 0.;

        int i = T_idx+1;
        bool heating_region_end = false;

        Real T_i, T_ip1, C_i, C_ip1;
        Real step, slope, sm1, oms;

        while (i < nbins_ && !heating_region_end) {
            if (net_cooling(i) > 0) { // if this bin has net cooling
                temps(i) = get_zero_point(i-1,rho);          
                net_cooling(i) = -1e-20;
                end_linear_idx = i-1;
                heating_region_end = true;
            }
                
            T_i = temps(i-1);
            C_i = net_cooling(i-1);

            T_ip1 = temps(i);
            C_ip1 = net_cooling(i);
            
            if (i == start_linear_idx || i == end_linear_idx) {
                slope = (T_ip1-T_i)/(C_ip1-C_i);
                step = (C_ref/T_ref)*slope*std::log(C_i/C_ip1);
            } else { 
                slope = std::log(C_ip1/C_i)/std::log(T_ip1/T_i);
                if (isClose(slope,1.0)) {
                    step = (C_ref/C_i)*(T_i/T_ref)*std::log(T_i/T_ip1);
                } else {
                    sm1  = slope - 1.0;
                    step = (C_ref/C_i)*(T_i/T_ref)*(std::pow(T_i/T_ip1,sm1)-1)/sm1;
                }
            }

            Ys(i) = Ys(i-1) + step;
            i += 1;
        }

        Real T_max = temps(i-1);
        Real Y_lim = Ys(i-1);

        // Compute current TEF Y
        Real Y, Y_i;
        
        Y_i = Ys(T_idx);
        T_i = temps(T_idx);
        C_i = net_cooling(T_idx);
        T_ip1 = temps(T_idx+1);
        C_ip1 = net_cooling(T_idx+1);


        if (T_idx == start_linear_idx || T_idx == end_linear_idx) {
            slope = (T_ip1-T_i)/(C_ip1-C_i);
            Y = Y_i + (C_ref/T_ref)*slope*std::log(C_i/(C_i+(T-T_i)/slope));
        } else {
            slope = std::log(C_ip1/C_i)/std::log(T_ip1/T_i);
            if (isClose(slope,1.0)) {
                Y = Y_i + (C_ref/C_i)*(T_i/T_ref)*std::log(T_i/T);
            } else {
                sm1  = slope - 1.0;
                Y = Y_i + (C_ref/C_i)*(T_i/T_ref)*(std::pow(T_i/T,sm1)-1)/sm1;
            }
        }

        // Compute the new TEF after a timestep (Eqn. 26)
        Real Y_new = Y + rho*C_ref*const_factor_*dt/T_ref;

        // Check if we have reached the minimum temperature
        if (Y_new < Y_lim) { return T_max; }

        // TEF is a strictly decreasing function of T and new_tef < tef
        // Check if the new TEF falls into a higher bin
        // If so, update slopes and coefficients
        while ((T_idx < nbins_) && (Y_new < Ys(T_idx+1))) { T_idx += 1; }
        
        T_i   = temps(T_idx);
        Y_i   = Ys(T_idx);
        C_i   = net_cooling(T_idx);
        T_ip1 = temps(T_idx+1);
        C_ip1 = net_cooling(T_idx+1);
        
        // Compute the Inverse Temporal Evolution Function Y^{-1}(Y) (Eq. A7)
        if (T_idx == start_linear_idx || T_idx == end_linear_idx) {
            slope = (C_ip1-C_i)/(T_ip1-T_i);
            T_new = T_i + (C_i/slope)*(1/std::exp((Y_new-Y_i)*(T_ref/C_ref)*slope)-1);
        } else {
            slope = std::log(C_ip1/C_i)/std::log(T_ip1/T_i);
            if (isClose(slope,1.0)) {
                T_new = T_i*std::exp(-(C_i/C_ref)*(T_ref/T_i)*(Y_new-Y_i));
            } else {
                oms  = 1.0 - slope;
                T_new = T_i*std::pow(1-oms*(C_i/C_ref)*(T_ref/T_i)*(Y_new-Y_i),1/oms);
            }
        }
    }

    return T_new;
}

int Cooling::get_temp_index(const Real T) {
    int T_idx = 0; // Get index of our temperature bin
    while ((T_idx < nbins_-2) && (temperature_table(T_idx+1) < T)) { 
        T_idx += 1; 
    }
    return T_idx;
}

Real Cooling::inst_cooling(const Real T, const Real rho) {
    int T_idx = get_temp_index(T);

    Real cooling_k = (std::log(cooling_table(T_idx+1)/cooling_table(T_idx))
                     /std::log(temperature_table(T_idx+1)/temperature_table(T_idx)));
    Real heating_k = (std::log(heating_table(T_idx+1)/heating_table(T_idx))
                     /std::log(temperature_table(T_idx+1)/temperature_table(T_idx)));
    Real cool_T = cooling_table(T_idx)*std::pow((T/temperature_table(T_idx)),cooling_k);
    Real heat_T = heating_table(T_idx)*std::pow((T/temperature_table(T_idx)),heating_k);
    return cool_T - heat_T/rho;
}

Real Cooling::get_zero_point(const int i, const Real rho) {
    Real log_T_ratio = std::log(temperature_table(i+1)/temperature_table(i));
    Real sc = std::log(cooling_table(i+1)/cooling_table(i))/log_T_ratio;
    Real sh = std::log(heating_table(i+1)/heating_table(i))/log_T_ratio;
    return temperature_table(i)*std::pow(heating_table(i)/(rho*cooling_table(i)),1/(sc-sh));
}

// Compute (1.0e-23)*(g-1.0)*mu/(kb*mu_e*mu_h*mh);
Real Cooling::compute_constant_factor() {
  Real mh  = 1.6605e-24;       // Atomic mass unit (g)
  Real kb  = 1.380648e-16;     // Boltzmann constant (erg/K)
  Real g  = 5.0/3.0;           // Adiabatic index
  Real X = 0.7; Real Z = 0.02; // Fully ionized, solar abundances
  Real mu   = 1.0/(2.0*X + 0.75*(1.0-X-Z) + Z/2.0);
  Real mu_e = 2.0/(1.0+X);
  Real mu_h = 1.0/X;

  return (1.0e-23)*(g-1.0)*mu/(kb*mu_e*mu_h*mh);
}

bool Cooling::isClose(const Real a, const Real b, const Real tol) {
    return std::fabs(a-b) <= tol;
}

Real Cooling::single_point_cooling_time(const Real T, const Real rho) {
    return T/(const_factor_*rho*inst_cooling(T,rho));
}

#endif // COOLING_HPP_
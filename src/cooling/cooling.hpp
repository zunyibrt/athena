#ifndef COOLING_HPP_
#define COOLING_HPP_

#include "../athena.hpp"
#include <fstream>
#include <sstream>

class Cooling {
public:
    Cooling();
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
    Real Get_tceil()  const {return Tceil_;}
    Real Get_tfloor() const {return Tfloor_;}

private:
    // Loaded Arrays
    std::vector<Real> temperature_table;
    std::vector<Real> cooling_table;
    std::vector<Real> heating_table;
    int nbins_;

    // Userful physical quantities 
    Real Tceil_;
    Real Tfloor_;
    Real const_factor_;

    // Scratch Arrays
    std::vector<Real> net_cooling;
    std::vector<Real> temps;
    std::vector<Real> Ys; 

    // Helper
    int get_temp_index(const Real T);
    Real inst_cooling(const Real T, const Real rho);
    Real get_zero_point(const int i, const Real rho);
    Real compute_constant_factor();
    bool isClose(const Real a, const Real b, const Real tol = 1e-20);
};

// Constructor
Cooling::Cooling() {
    nbins_ = 1;

    temperature_table = {0.};
    cooling_table     = {0.};
    heating_table     = {0.};

    Tceil_        = 0.;
    Tfloor_       = 0.;
    const_factor_ = 0.;

    net_cooling = {0.};
    temps       = {0.};
    Ys          = {0.};

    return;
}

Cooling::Cooling(const std::string in) {
    // Read in table
    std::ifstream infile;
    infile.open(in.c_str());
    std::string line;

    Real X    = 0.7;
    Real mu_e = 2.0/(1.0+X);
    Real mh   = 1.6605e-24;

    while (std::getline(infile,line)) {
        if (line[0] != '#') {
            std::stringstream ss(line);
            Real T,c,h;
            ss >> T >> c >> h;
            temperature_table.push_back(T);
            cooling_table.push_back(std::max(c/1.0e-23,1.0e-30));
            heating_table.push_back(std::max(mu_e*mh*h/1.0e-23,1.0e-30));
        }
    }
    infile.close();

    nbins_ = temperature_table.size();
    // Assert that they have the same lengths
    if (nbins_ != cooling_table.size() || nbins_ != heating_table.size()) { 
        std::stringstream msg;
        msg << "### ERROR : Heating/Cooling Tables are different sizes!" << std::endl;
        ATHENA_ERROR(msg); 
    }

    Tceil_        = temperature_table.at(nbins_-1);
    Tfloor_       = temperature_table.at(0);
    const_factor_ = compute_constant_factor();

    net_cooling.resize(nbins_,0);
    temps.resize(nbins_,0);
    Ys.resize(nbins_,0);

    return;
}

// Destructor
Cooling::~Cooling() {
    return;
}

// Exact Integration Scheme for Radiative Cooling from Townsend (2009)
// Expanded to include heating
// Requires: Input temperature, density, and timestep in cgs units
// Returns: Temperature(K) at the next timestep after cooling/heating
Real Cooling::townsend_cooling(const Real T, const Real rho, const Real dt) {
    // Ensure temperature is within bounds
    if (T <= Tfloor_) return Tfloor_;
    if (T >= Tceil_)  return Tceil_;

    // Make a copy of temperatures that can be edited
    temps = temperature_table;

    // Calculate net cooling for a given density
    for (int i=0; i<nbins_; ++i) {
        net_cooling.at(i) = cooling_table.at(i)-heating_table.at(i)/rho;
    }

    // Get index of the temperature bin
    int T_idx = get_temp_index(T);

    // Check if there is net cooling or heating
    Real net_cool_T = inst_cooling(T, rho);

    // If no net cooling/heating, temperature does not change
    if (isClose(net_cool_T,0.)) return T;

    // Flags for piece-wise linear bins - only at start or end
    int start_linear_idx = -1;
    int end_linear_idx   = -1;

    // Our aim is to calculate the new temperature
    Real T_new = 0.0;

    if (net_cool_T > 0) { // Cooling
        // We save time by using the next bin as the reference for the TEF
        // Hence we need to check for the edge case where it is a heating bin
        if (net_cooling.at(T_idx+1) < 0) {
            // If so we figure out what the equilibruim temperature is
            // Note everything after this index should not be used anymore
            // We then hijack the next bin to mark this temperature
            // And set the cooling there to zero
            // The piecewise power law assumption is still here, so
            //   cooling in this bin will be automatically suppressed
            temps.at(T_idx+1)       = get_zero_point(T_idx,rho);
            net_cooling.at(T_idx+1) = 1e-10;
            start_linear_idx        = T_idx;
        }

        // Set reference values
        Real T_ref = temps.at(T_idx+1);
        Real C_ref = net_cooling.at(T_idx+1);
        
        // Compute Y from idx = temp_bin+1 downwards
        // Calculate Y_k recursively (Eq. A6) from idx = temp_bin+1 downwards    
        Ys.at(T_idx+1) = 0.0;

        int i = T_idx;
        bool cooling_region_end = false;

        Real T_i, T_ip1, C_i, C_ip1;
        Real step, slope, sm1, oms;

        while (i >= 0 && !cooling_region_end) {
            if (net_cooling.at(i) < 0) { // if this bin has net heating
                temps.at(i) = get_zero_point(i,rho);          
                net_cooling.at(i) = 1e-10;
                end_linear_idx = i;
                cooling_region_end = true;
            }
                
            T_i = temps.at(i);
            C_i = net_cooling.at(i);

            T_ip1 = temps.at(i+1);
            C_ip1 = net_cooling.at(i+1);
            
            if (i == start_linear_idx || i == end_linear_idx) {
                slope = (T_ip1-T_i)/(C_ip1-C_i); // for our use case, c_ip1 != c_i
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

            Ys.at(i) = Ys.at(i+1) - step;
            i -= 1;
        }

        Real T_min = temps.at(i+1);
        Real Y_lim = Ys.at(i+1);

        // Compute current TEF Y
        Real Y, Y_i;
        
        Y_i = Ys.at(T_idx);
        T_i = temps.at(T_idx);
        C_i = net_cooling.at(T_idx);
        T_ip1 = temps.at(T_idx+1);
        C_ip1 = net_cooling.at(T_idx+1);


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
        while ((T_idx > 0) && (Y_new > Ys.at(T_idx))) { T_idx -= 1; }
        
        T_i   = temps.at(T_idx);
        Y_i   = Ys.at(T_idx);
        C_i   = net_cooling.at(T_idx);
        T_ip1 = temps.at(T_idx+1);
        C_ip1 = net_cooling.at(T_idx+1);
        
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
        if (net_cooling.at(T_idx) > 0) { // If the cuurent bin is cooling and not heating
            temps.at(T_idx)       = get_zero_point(T_idx,rho);
            net_cooling.at(T_idx) = -1e-10;
            start_linear_idx      = T_idx;
        }
        
        // Set reference values
        Real T_ref = temps.at(T_idx);
        Real C_ref = net_cooling.at(T_idx);

        // Compute Y from idx = temp_bin upwards
        Ys.at(T_idx) = 0.;

        int i = T_idx+1;
        bool heating_region_end = false;

        Real T_i, T_ip1, C_i, C_ip1;
        Real step, slope, sm1, oms;

        while (i < nbins_ && !heating_region_end) {
            if (net_cooling.at(i) > 0) { // if this bin has net cooling
                temps.at(i) = get_zero_point(i-1,rho);          
                net_cooling.at(i) = -1e-10;
                end_linear_idx = i-1;
                heating_region_end = true;
            }
                
            T_i = temps.at(i-1);
            C_i = net_cooling.at(i-1);

            T_ip1 = temps.at(i);
            C_ip1 = net_cooling.at(i);
            
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

            Ys.at(i) = Ys.at(i-1) + step;
            i += 1;
        }

        Real T_max = temps.at(i-1);
        Real Y_lim = Ys.at(i-1);

        // Compute current TEF Y
        Real Y, Y_i;
        
        Y_i = Ys.at(T_idx);
        T_i = temps.at(T_idx);
        C_i = net_cooling.at(T_idx);
        T_ip1 = temps.at(T_idx+1);
        C_ip1 = net_cooling.at(T_idx+1);

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
        while ((T_idx < nbins_) && (Y_new < Ys.at(T_idx+1))) { T_idx += 1; }
        
        T_i   = temps.at(T_idx);
        Y_i   = Ys.at(T_idx);
        C_i   = net_cooling.at(T_idx);
        T_ip1 = temps.at(T_idx+1);
        C_ip1 = net_cooling.at(T_idx+1);
        
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
    while ((T_idx < nbins_-2) && (temperature_table.at(T_idx+1) < T)) { 
        T_idx += 1; 
    }
    return T_idx;
}

Real Cooling::inst_cooling(const Real T, const Real rho) {
    int T_idx = get_temp_index(T);

    Real cooling_k = (std::log(cooling_table.at(T_idx+1)/cooling_table.at(T_idx))
                     /std::log(temperature_table.at(T_idx+1)/temperature_table.at(T_idx)));
    Real heating_k = (std::log(heating_table.at(T_idx+1)/heating_table.at(T_idx))
                     /std::log(temperature_table.at(T_idx+1)/temperature_table.at(T_idx)));
    Real cool_T = cooling_table.at(T_idx)*std::pow((T/temperature_table.at(T_idx)),cooling_k);
    Real heat_T = heating_table.at(T_idx)*std::pow((T/temperature_table.at(T_idx)),heating_k);
    
    return cool_T - heat_T/rho;
}

Real Cooling::get_zero_point(const int i, const Real rho) {
    Real log_T_ratio = std::log(temperature_table.at(i+1)/temperature_table.at(i));
    Real sc = std::log(cooling_table.at(i+1)/cooling_table.at(i))/log_T_ratio;
    Real sh = std::log(heating_table.at(i+1)/heating_table.at(i))/log_T_ratio;
    return temperature_table.at(i)*std::pow(heating_table.at(i)/(rho*cooling_table.at(i)),1/(sc-sh));
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
    return std::abs(T/(const_factor_*rho*inst_cooling(T,rho)));
}

#endif // COOLING_HPP_
// Problem generator for studying clouds in SN-driven winds

// C++ headers
#include <algorithm>  // min()
#include <cmath>      // abs(), pow(), sqrt()
#include <sstream>    // stringstream
#include <stdexcept>  // runtime_error
#include <string>     // string

// Athena++ headers
#include "../mesh/mesh.hpp"
#include "../athena.hpp"                   // macros, enums, declarations
#include "../athena_arrays.hpp"            // AthenaArray
#include "../globals.hpp"                  // Globals
#include "../parameter_input.hpp"          // ParameterInput
#include "../coordinates/coordinates.hpp"  // Coordinates
#include "../eos/eos.hpp"                  // EquationOfState
#include "../field/field.hpp"              // Field
#include "../hydro/hydro.hpp"              // Hydro
#include "../fft/athena_fft.hpp"           // FFT
#include "../utils/utils.hpp"              // Utilities
#include "../utils/sn_injection.hpp"       // SNInj
#include "../cooling/cooling.hpp"          // Cooling

// Code Units (cgs)
const Real amu       = 1.6605e-24;   // atomic mass unit (g)
const Real kb        = 1.380648e-16; // boltzmann constant (erg/K)
const Real mu        = 0.6173;       // mean molecular weight (Solar)
                                         // X = 0.7; Z = 0.02

const Real unit_len  = 3.086e20;     // 100 pc
const Real unit_temp = 1.0e4;        // 10^4 K
const Real unit_n    = 1.0;          // cm^-3

const Real unit_rho  = unit_n * amu * mu;
const Real unit_mass = unit_rho * std::pow(unit_len,3);
const Real unit_pres = kb * unit_n * unit_temp;
const Real unit_engy = unit_pres * std::pow(unit_len,3);
const Real unit_vel  = std::sqrt(unit_pres / unit_rho);
const Real unit_time = unit_len / unit_vel;
const Real unit_kap  = (unit_pres * unit_len * unit_len) /
                       (unit_time * unit_temp);
const Real unit_gam  = (unit_pres/unit_time)/(unit_n*unit_n);

// Global variables
static Real rho0, pgas0;
static Real vcir, R0;

static Cooling cooler;

static SNInj injector;
static std::vector<Real> sn_times;
static int next_sn_idx;
static Real r_inj,e_sn,m_ej;



// User defined boundary conditions 
void NoInflowInnerX3(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &a,
                     FaceField &b, AthenaArray<Real> &r,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, 
                     int kl, int ku, int ngh);
void NoInflowOuterX3(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &a,
                     FaceField &b, AthenaArray<Real> &r,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, 
                     int kl, int ku, int ngh);

// User defined source functions
void SourceFunctions(MeshBlock *pmb, const Real time, const Real dt,
                     const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &prim_scalar,
                     const AthenaArray<Real> &bcc,
                     AthenaArray<Real> &cons,
                     AthenaArray<Real> &cons_scalar);
void GravitySource(MeshBlock *pmb, const Real dt, 
                   const AthenaArray<Real> &prim, 
                   AthenaArray<Real> &cons);
Real GravAccel(Real z);
void CoolingSource(MeshBlock *pmb, const Real dt, 
                   const AthenaArray<Real> &prim, 
                   AthenaArray<Real> &cons,
                   const AthenaArray<Real> &bcc);
void SNSource(MeshBlock *pmb, const Real dt, 
              const AthenaArray<Real> &prim, 
              AthenaArray<Real> &cons);

// User defined history functions
Real CalculateSNEnergyInjection(MeshBlock *pmb, int iout);
Real CalculateSNMassInjection(MeshBlock *pmb, int iout);
Real CalculateColdGasMass(MeshBlock *pmb, int iout);

// Misc
void AssertCondition(bool condition, std::string msg);
static Real CoolingTimestep(MeshBlock *pmb); // User defined time step

//===========================================================================//
//                             Initializations                               //
//===========================================================================//
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Check for intended configuration
  AssertCondition(NON_BAROTROPIC_EOS,"Use adiabatic eos");
  AssertCondition(mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user"),
                  "Use user-defined inner boundary condition along x3");
  AssertCondition(mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user"),
                  "Use user-defined outer boundary condition along x3");
  AssertCondition(mesh_size.nx2 > 1 && mesh_size.nx3 > 1,"Run in 3D");
#ifndef FFT
  AssertCondition(false,"Compile with FFT");
#endif

  // Read in initial conditions
  rho0  = pin->GetReal("problem", "rho0");
  pgas0 = pin->GetReal("problem", "pgas0");
  vcir  = pin->GetReal("problem", "vcir");
  R0    = pin->GetReal("problem", "R0");

  // Initialize cooling, SN injection, and turbulence
  cooler    = Cooling(pin->GetString("cooling","cooling_table"));
  injector  = SNInj();
  turb_flag = pin->GetInteger("problem","turb_flag");

  // Generate SN times
  next_sn_idx = 0;
  const Real cluster_mass = pin->GetOrAddReal("SN","M_cluster",1e4);
  sn_times = injector.GetSNTimes(cluster_mass,
                                 pin->GetReal("SN","tstart"), 
                                 pin->GetReal("time","tlim"),
                                 unit_time);

  // Compute energy and mass injection densities
  r_inj = pin->GetReal("SN","r_inj"); // Input in code unis
  const Real sphere_vol = (4.0/3.0)*PI*std::pow(r_inj,3);

  const Real E_def = 1e51/unit_engy; // Default 10^51 ergs
  const Real M_def = 8.4*1.988e33/unit_mass; // Default 8.4 solar masses
  e_sn  = pin->GetOrAddReal("SN","E_sn",E_def)/sphere_vol; // Input in ergs
  m_ej  = pin->GetOrAddReal("SN","M_ej",M_def)/sphere_vol; // Input in solar mass

  // Enroll user-defined physical source terms
  EnrollUserExplicitSourceFunction(SourceFunctions);

  // Enroll user-defined boundary conditions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, NoInflowInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, NoInflowOuterX3);

  // Enroll timestep so that dt <= min t_cool
  EnrollUserTimeStepFunction(CoolingTimestep);

  // Enroll user-defined history outputs
  AllocateUserHistoryOutput(3);
  EnrollUserHistoryOutput(0,CalculateSNEnergyInjection,"SNEnergyInjection");
  EnrollUserHistoryOutput(1,CalculateSNMassInjection,"SNMassInjection");
  EnrollUserHistoryOutput(2,CalculateColdGasMass,"CalculateColdGasMass");

  // Output initialization information
  if (Globals::my_rank == 0) {
    std::cout << "==============================================" << std::endl;
    std::cout << "Units                                         " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Unit Length         : " << unit_len*3.241e-19  
                                          << " pc"    << std::endl;
    std::cout << "Unit Temperature    : " << unit_temp 
                                          << " K"     << std::endl;     
    std::cout << "Unit Number Density : " << unit_n 
                                          << " cm^-3" << std::endl;        
    std::cout << "Unit Velocity       : " << unit_vel*1e-5 
                                          << " km/s"  << std::endl;
    std::cout << "Unit Time           : " << unit_time*3.171e-14 
                                          << " Myr"   << std::endl;
    std::cout << std::endl;

    std::cout << "==============================================" << std::endl;
    std::cout << "Initialization                                " << std::endl;
    std::cout << "==============================================" << std::endl;
    std::cout << "Midplane Density     : " << rho0  
                                           << " cm^-3" << std::endl;  
    std::cout << "Midplane Temperature : " << unit_temp*pgas0/rho0 
                                           << " K"     << std::endl; 
    std::cout << "Circular Velocity    : " << vcir*unit_vel*1e-5   
                                           << " km/s"  << std::endl; 
    std::cout << "Galactic Radius      : " << R0*unit_len*3.241e-19    
                                           << " pc"    << std::endl; 
    std::cout << std::endl;
    std::cout << "Turbulence Flag  : " << turb_flag << std::endl;
    std::cout << "Cooling Table    : "
              << pin->GetString("cooling","cooling_table")     << std::endl;
    std::cout << "Cluster Mass     : " << cluster_mass        
                                       << " Solar Masses"      << std::endl;
    std::cout << "First SN at      : " << sn_times.at(0)       
                                       << " Code Units"        << std::endl;
    std::cout << "Last SN at       : " << sn_times.rbegin()[1] 
                                       << " Code Units"        << std::endl;
    std::cout << "No. of SN        : " << sn_times.size()-1    << std::endl;
    std::cout << "Energy per SN    : " << e_sn*sphere_vol*unit_engy      
                                       << " ergs"              << std::endl;
    std::cout << "Mass per SN      : " << m_ej*sphere_vol*unit_mass/1.988e33      
                                       << " Solar Masses"      << std::endl;  
    std::cout << "Injection Radius : " << r_inj*unit_len*3.241e-19
                                       << " pc"                << std::endl;
    std::cout << std::endl;
  }

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(4);

  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  const Real cs2 = pgas0/rho0; // Isothermal sound speed
  const Real H = std::sqrt(cs2)/(vcir/R0); // Scale height = cs/omega

  const Real amp = 0.01;
  std::int64_t iseed = -1 - gid;


  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real z = pcoord->x3v(k);

        // Hydrostatic solution
        Real rho = rho0*std::exp(SQR(R0/H)*
                                 (std::pow(SQR(R0)/(SQR(R0)+SQR(z)),0.5)-1));

        // Grid scale random field for initial TI
        rho *= (1 + amp*(ran2(&iseed)-0.5));
        rho = std::fmax(rho,1e-8);

        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = rho*amp*(ran2(&iseed)-0.5);
        phydro->u(IM2,k,j,i) = rho*amp*(ran2(&iseed)-0.5);
        phydro->u(IM3,k,j,i) = rho*amp*(ran2(&iseed)-0.5);

        phydro->u(IEN,k,j,i) =  phydro->u(IDN,k,j,i)*cs2/(peos->GetGamma()-1);
        phydro->u(IEN,k,j,i) += 0.5*(SQR(phydro->u(IM1,k,j,i)) +
                                     SQR(phydro->u(IM2,k,j,i)) +
                                     SQR(phydro->u(IM3,k,j,i)))
                                     /phydro->u(IDN,k,j,i);
      }
       
    }
  }


  // initialize uniform interface B
  if (MAGNETIC_FIELDS_ENABLED) {
    Real beta = pin->GetReal("problem","beta");
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie+1; i++) {
      pfield->b.x1f(k,j,i) = std::sqrt(2*phydro->u(IDN,k,j,i)*cs2/beta);
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je+1; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x2f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke+1; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      pfield->b.x3f(k,j,i) = 0.0;
    }}}
    for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
    for (int i=is; i<=ie; i++) {
      phydro->u(IEN,k,j,i) += 0.5*(SQR(pfield->b.x1f(k,j,i)) +
                                   SQR(pfield->b.x2f(k,j,i)) +
                                   SQR(pfield->b.x3f(k,j,i)));
    }}}
  }


  return;
}

//===========================================================================//
//                           Boundary Conditions                             //
//===========================================================================//
// Zero gradient no inflow upper boundary condition
void NoInflowOuterX3(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim,
                     FaceField &b, AthenaArray<Real> &r,
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, 
                     int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        prim(IDN,ku+k,j,i) = prim(IDN,ku,j,i);
        prim(IPR,ku+k,j,i) = prim(IPR,ku,j,i);

        prim(IVX,ku+k,j,i) = prim(IVX,ku,j,i); 
        prim(IVY,ku+k,j,i) = prim(IVY,ku,j,i);
        prim(IVZ,ku+k,j,i) = std::fmax(0.,prim(IVZ,ku,j,i)); 
      }
    }
  }

  return;
}

// Zero gradient no inflow lower boundary condition
void NoInflowInnerX3(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &prim,
                     FaceField &b, AthenaArray<Real> &r, 
                     Real time, Real dt,
                     int il, int iu, int jl, int ju, 
                     int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        prim(IDN,kl-k,j,i) = prim(IDN,kl,j,i);
        prim(IPR,kl-k,j,i) = prim(IPR,kl,j,i);

        prim(IVX,kl-k,j,i) = prim(IVX,kl,j,i); 
        prim(IVY,kl-k,j,i) = prim(IVY,kl,j,i);
        prim(IVZ,kl-k,j,i) = std::fmin(0.,prim(IVZ,kl,j,i)); 
      }
    }
  }

  return;
}

//===========================================================================//
//                              Source Terms                                 //
//===========================================================================//
// Source terms
void SourceFunctions(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim, 
                    const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar) {
  // Vertical gravity
  GravitySource(pmb,dt,prim,cons);

  // Radiative cooling
  CoolingSource(pmb,dt,prim,cons,bcc);

  // SNe injection
  if (time > sn_times.at(next_sn_idx)) { // Step through list of SN times
    SNSource(pmb,dt,prim,cons);
    next_sn_idx++;
  }
  
  return;
}

// Gravity source term
void GravitySource(MeshBlock *pmb, const Real dt, 
                   const AthenaArray<Real> &prim, 
                   AthenaArray<Real> &cons) {
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    Real z = pmb->pcoord->x3v(k);
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real den = prim(IDN,k,j,i);
        Real vz  = prim(IVZ,k,j,i);
        
        Real dvz = dt*GravAccel(z);

        cons(IM3,k,j,i) -= den*dvz;
        cons(IEN,k,j,i) -= 0.5*den*(2*dvz*vz - SQR(dvz));
      }
    }
  }

  return;
}

// The gravitational profile
Real GravAccel(Real z) {
  // Equivalent to g = z * GM/(R^2 + z^2)^1.5
  return SQR(vcir/R0)*std::pow(SQR(R0)/(SQR(R0)+SQR(z)),1.5)*z;
}

// Cooling source term
void CoolingSource(MeshBlock *pmb, const Real dt, 
                   const AthenaArray<Real> &prim, 
                   AthenaArray<Real> &cons,
                   const AthenaArray<Real> &bcc) {
  const Real g = pmb->peos->GetGamma();

  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        // Use cons not prime because we do not need intermediate step 
        // to calculate cooling
        Real rho  = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i)
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                          + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                          + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho; 

        if (MAGNETIC_FIELDS_ENABLED) {
          eint -= 0.5 *(bcc(IB1,k,j,i) * bcc(IB1,k,j,i)
                      + bcc(IB2,k,j,i) * bcc(IB2,k,j,i)
                      + bcc(IB3,k,j,i) * bcc(IB3,k,j,i));
        }       

        // T = P/rho
        Real temp = eint * (g-1.0)/rho;
        
        Real temp_new = temp;

        // Calculate new temperature using the Townsend Algorithm
        // The inputs should be given in cgs units
        Real temp_cgs = temp * unit_temp;
        Real rho_cgs  = rho  * unit_rho;
        Real dt_cgs   = dt   * unit_time;

        temp_new = cooler.townsend_cooling(temp_cgs,rho_cgs,dt_cgs)/unit_temp;

        // Average Neighbours if near the floor
        if (temp_new < 1.1*cooler.Get_tfloor()) {
          Real sum_temp = 0.0;
          Real n_nb = 0.0;
          for (int nk=std::max(k-1,pmb->ks); nk<=std::min(k+1,pmb->ke); ++k) {
            for (int nj=std::max(j-1,pmb->js); nj<=std::min(j+1,pmb->je); ++j) {
              for (int ni=std::max(i-1,pmb->is); ni<=std::min(i+1,pmb->ie); ++i) {
                if ((nk != k) || (nj != j) || (ni != i)) {
                  auto nrho = cons(IDN,nk,nj,ni);
                  auto neint = cons(IEN,nk,nj,ni)
                              - 0.5 *(cons(IM1,nk,nj,ni)*cons(IM1,nk,nj,ni)
                                    + cons(IM2,nk,nj,ni)*cons(IM2,nk,nj,ni)
                                    + cons(IM3,nk,nj,ni)*cons(IM3,nk,nj,ni))/nrho;
                  if (MAGNETIC_FIELDS_ENABLED) {
                    neint -= 0.5 *(bcc(IB1,nk,nj,ni) * bcc(IB1,nk,nj,ni)
                                 + bcc(IB2,nk,nj,ni) * bcc(IB2,nk,nj,ni)
                                 + bcc(IB3,nk,nj,ni) * bcc(IB3,nk,nj,ni));
                  }

                  // T = P/rho
                  auto ntemp = neint * (g-1.0)/nrho;

                  n_nb += 1;
                  sum_temp += ntemp;
                }
          }}}
          temp_new = sum_temp/n_nb;

          std::cout << " Temperature Floor Hit... Averaging Neigbours... T = " << temp_new << std::endl;
        }

        // Update energy based on change in temperature
        cons(IEN,k,j,i) += (temp_new - temp) * (rho/(g-1.0));

        // Store the change in energy/time in a user defined output variable
        pmb->user_out_var(2,k,j,i) = (temp_new - temp) * (rho/(g-1.0)) / dt;
      }
    }
  }

  return;
}

// SN source term
void SNSource(MeshBlock *pmb, const Real dt, 
              const AthenaArray<Real> &prim, 
              AthenaArray<Real> &cons) {
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    Real z = pmb->pcoord->x3v(k);
    for (int j=pmb->js; j<=pmb->je; j++) {
      Real y = pmb->pcoord->x2v(j);
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real x = pmb->pcoord->x1v(i);
        
        // If within injection sphere
        if (SQR(x) + SQR(y) + SQR(z) < SQR(r_inj)) {
          // Inject energy and mass into cell
          cons(IEN,k,j,i) += e_sn;
          cons(IDN,k,j,i) += m_ej;

          // Store the changes in a user defined output variable
          pmb->user_out_var(0,k,j,i) += e_sn;
          pmb->user_out_var(1,k,j,i) += m_ej;
        }        
      }
    }
  }

  return;
}

//===========================================================================//
//                                Analysis                                   //
//===========================================================================//
void MeshBlock::UserWorkBeforeOutput(ParameterInput *pin) {
  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real T = unit_temp*phydro->w(IPR,k,j,i)/phydro->w(IDN,k,j,i);
        Real rho = unit_rho*phydro->w(IDN,k,j,i);
        Real tcool = cooler.single_point_cooling_time(T, rho)/unit_time;

        user_out_var(3,k,j,i) = tcool;
      }
    }
  }
}

// History function for the total SN energy injected
Real CalculateSNEnergyInjection(MeshBlock *pmb, int iout) {
  AthenaArray<Real> vol(pmb->ncells1);
  Real sum = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
      for (int i=pmb->is; i<=pmb->ie; i++) {
       sum += pmb->user_out_var(0,k,j,i)*vol(i);
      }
    }
  }

  return sum;
}

// History function for the total SN mass injected
Real CalculateSNMassInjection(MeshBlock *pmb, int iout) {
  AthenaArray<Real> vol(pmb->ncells1);
  Real sum = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
      for (int i=pmb->is; i<=pmb->ie; i++) {
       sum += pmb->user_out_var(1,k,j,i)*vol(i);
      }
    }
  }

  return sum;
}

// History function for the total SN mass injected
Real CalculateColdGasMass(MeshBlock *pmb, int iout) {
  AthenaArray<Real> vol(pmb->ncells1);
  Real sum = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      pmb->pcoord->CellVolume(k,j,pmb->is,pmb->ie,vol);
      for (int i=pmb->is; i<=pmb->ie; i++) {
       Real T = unit_temp*pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i);
       if (T < 2e4) {
         sum += pmb->phydro->w(IDN,k,j,i)*vol(i);
       }
      }
    }
  }

  return sum;
}

//===========================================================================//
//                                  Misc                                     //
//===========================================================================//
// Wrapper for throwing error messages when a condition is not met
void AssertCondition(bool condition, std::string error_message) {
  std::stringstream msg;
  msg << "### FATAL ERROR : " << error_message << std::endl;
  if (!condition) { ATHENA_ERROR(msg); }

  return;
}

// Compute the minimum cooling timstep
Real CoolingTimestep(MeshBlock *pmb) {
  Real min_dt = 1.0e10;
  for (int k=pmb->ks; k<=pmb->ke; k++) {
    for (int j=pmb->js; j<=pmb->je; j++) {
      for (int i=pmb->is; i<=pmb->ie; i++) {
        Real T = unit_temp*pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i);
        Real rho = unit_rho*pmb->phydro->w(IDN,k,j,i);
        Real tcool = cooler.single_point_cooling_time(T, rho)/unit_time;
 
        if (tcool > 1e-10) {
          min_dt = std::fmin(min_dt, tcool);
        } else {
          std::cout << "Bad Cell, Density: " << pmb->phydro->w(IDN,k,j,i) 
                          << "  Pressure : " << pmb->phydro->w(IPR,k,j,i) 
                     << "   Cooling Time : " << tcool 
                                  << "  z: " << pmb->pcoord->x3v(k) 
                                             << std::endl;
        }
      }
    }
  }

  return 0.5*min_dt;
}

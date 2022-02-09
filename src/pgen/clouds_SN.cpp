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

// Code Units (cgs)
static Real const mh        = 1.6605e-24;   // mass (g)
static Real const kb        = 1.380648e-16; // boltzmann constant (erg/K)
static Real const mu        = 0.6173;       // mean molecular weight (Solar)
                                            // X = 0.7; Z = 0.02

static Real const unit_len  = 3.086e20;     // 100 pc
static Real const unit_temp = 1.0e4;        // 10^4 K
static Real const unit_n    = 1.0;          // cm^-3

static Real const unit_rho  = unit_n * mh * mu;
static Real const unit_pres = kb * unit_n * unit_temp;
static Real const unit_vel  = sqrt(unit_pres / unit_rho);
static Real const unit_time = unit_len / unit_vel;
static Real const unit_kap  = (unit_pres * unit_len * unit_len) /
                              (unit_time * unit_temp);
static Real const unit_gam  = (unit_pres/unit_time)/(unit_n*unit_n);

// Global variables
static Real rho_0, pgas_0;
static Real grav_strength; // g(z) = grav_strength*z [1/time^2]

// User defined boundary conditions 
// User defined source functions
void SourceFunctions(MeshBlock *pmb, const Real time, const Real dt,
                     const AthenaArray<Real> &prim,
                     const AthenaArray<Real> &prim_scalar,
                     const AthenaArray<Real> &bcc,
                     AthenaArray<Real> &cons,
                     AthenaArray<Real> &cons_scalar);
void gravitySource(MeshBlock *pmb, const Real dt, 
                   const AthenaArray<Real> &prim, 
                   AthenaArray<Real> &cons);
Real grav_accel(Real z);

// User defined history functions

// Misc
void assertCondition(bool condition, std::string msg);

//===========================================================================//
//                             Initializations                               //
//===========================================================================//
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Check for specific configuration
  assertCondition(NON_BAROTROPIC_EOS,"Use adiabatic eos");

  // Read in constants
  rho_0         = pin->GetReal("problem", "rho_0");
  pgas_0        = pin->GetReal("problem", "pgas_0");
  grav_strength = pin->GetReal("problem", "grav_strength");

  // Enroll user-defined physical source terms
  EnrollUserExplicitSourceFunction(SourceFunctions);

  // Enroll user-defined boundary conditions
  // Enroll user defined history outputs

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  // AllocateUserOutputVariables(n);

  return;
}

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  // Check for specific configuration
  assertCondition(block_size.nx2 > 1 && block_size.nx3 > 1,"Run in 3D");

  Real cs2 = pgas_0/rho_0;
  Real gamma = peos->GetGamma();

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        // Hydrostatic solution for g(z) = const*z
        Real rho = rho_0 * std::exp(-0.5*grav_strength*SQR(z)/cs2);

        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        phydro->u(IEN,k,j,i) = rho*cs2/(gamma-1);
      }
    }
  }

  return;
}

//===========================================================================//
//                           Boundary Conditions                             //
//===========================================================================//



//===========================================================================//
//                              Source Terms                                 //
//===========================================================================//
void SourceFunctions(MeshBlock *pmb, const Real time, const Real dt,
                    const AthenaArray<Real> &prim, 
                    const AthenaArray<Real> &prim_scalar,
                    const AthenaArray<Real> &bcc, AthenaArray<Real> &cons,
                    AthenaArray<Real> &cons_scalar) {
  // Vertical Gravity
  gravitySource(pmb,dt,prim,cons);

  // Cooling

  // SNe injection
  
  return;
}

void gravitySource(MeshBlock *pmb, const Real dt, 
                   const AthenaArray<Real> &prim, AthenaArray<Real> &cons) {
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real z = pmb->pcoord->x3v(k);
        Real den = prim(IDN,k,j,i);
        Real vz = prim(IVZ,k,j,i);
        
        Real dvz = dt*grav_accel(z);

        cons(IM3,k,j,i) -= den*dvz;
        cons(IEN,k,j,i) -= 0.5*den*(2*dvz*vz - SQR(dvz));
      }
    }
  }

  return;
}

Real grav_accel(Real z) {
  return grav_strength*z;
}

//===========================================================================//
//                                Analysis                                   //
//===========================================================================//

//===========================================================================//
//                                  Misc                                     //
//===========================================================================//
void assertCondition(bool condition, std::string error_message) {
  std::stringstream msg;
  msg << "### FATAL ERROR : " << error_message << std::endl;
  if (!condition) { ATHENA_ERROR(msg); }
}
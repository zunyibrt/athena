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
// User defined history functions

// Misc
void assertCondition(bool condition, string msg);

//===========================================================================//
//                             Initializations                               //
//===========================================================================//
void Mesh::InitUserMeshData(ParameterInput *pin) {
  // Check for specific configuration
  assertCondition(NON_BAROTROPIC_EOS,"Use adiabatic eos")
  assertCondition(block_size.nx2 > 1 && block_size.nx3 > 1,"Run in 3D")

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
  // Prepare index bounds (Is this needed?)
  int il = is - NGHOST; int iu = ie + NGHOST;
  int jl = js - NGHOST; int ju = je + NGHOST;
  int kl = ks - NGHOST; int ku = ke + NGHOST;

  Real cs2 = pgas_0/rho_0; // Sound speed

  if (Globals::my_rank == 0) {
        std::cout << " vc2o2r2 " << vc2o2r2 << "\n";
        std::cout << " cs2 " << cs2 << "\n";
        std::cout << " rho_0 " << rho_0 << "\n";
        std::cout << " rho_floor " << rho_floor << "\n";
  }

  // Initialize primitive values
  for (int k = kl; k <= ku; ++k) {
    for (int j = jl; j <= ju; ++j) {
      for (int i = il; i <= iu; ++i) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        Real rho = rho_0 * std::exp(-1.0*vc2o2r2*SQR(z)/cs2);
        rho = std::max(rho,rho_floor);

        phydro->w(IDN,k,j,i) = rho;
        phydro->w(IPR,k,j,i) = cs2 * rho; 
        phydro->w(IVX,k,j,i) = 0.0;
        phydro->w(IVY,k,j,i) = 0.0;
        phydro->w(IVZ,k,j,i) = 0.0;
      }
    }
  }

  // Initialize conserved values
  AthenaArray<Real> b;
  peos->PrimitiveToConserved(phydro->w, b, phydro->u, pcoord, 
                             il, iu, jl, ju, kl, ku);

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
        
        Real grav_acc = 42; // 2*vc2o2r2*z;
        Real dvz = dt*grav_acc;

        cons(IM3,k,j,i) -= den*dvz;
        cons(IEN,k,j,i) -= 0.5*den*(2*dvz*vz - SQR(dvz));
      }
    }
  }

  return;
}

//===========================================================================//
//                                Analysis                                   //
//===========================================================================//

//===========================================================================//
//                                  Misc                                     //
//===========================================================================//
void assertCondition(bool condition, string msg) {
  std::stringstream msg;
  msg << "### FATAL ERROR : " << msg << std::endl;
  if (!condition) { ATHENA_ERROR(msg); }
}
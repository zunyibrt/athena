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
static Real const amu       = 1.6605e-24;   // atomic mass unit (g)
static Real const kb        = 1.380648e-16; // boltzmann constant (erg/K)
static Real const mu        = 0.6173;       // mean molecular weight (Solar)
                                            // X = 0.7; Z = 0.02

static Real const unit_len  = 3.086e20;     // 100 pc
static Real const unit_temp = 1.0e4;        // 10^4 K
static Real const unit_n    = 1.0;          // cm^-3

// TODO: Output these in useful units
static Real const unit_rho  = unit_n * amu * mu;
static Real const unit_pres = kb * unit_n * unit_temp;
static Real const unit_vel  = sqrt(unit_pres / unit_rho);
static Real const unit_time = unit_len / unit_vel;
static Real const unit_kap  = (unit_pres * unit_len * unit_len) /
                              (unit_time * unit_temp);
static Real const unit_gam  = (unit_pres/unit_time)/(unit_n*unit_n);

// Global variables
static Real rho0, pgas0;
static Real vcir, R0;

// User defined boundary conditions 
void NoInflowInnerX3(MeshBlock *pmb, Coordinates *pco,
                     AthenaArray<Real> &a,
                     FaceField &b, Real time, Real dt,
                     int il, int iu, int jl, int ju, 
                     int kl, int ku, int ngh);
void NoInflowOuterX3(MeshBlock *pmb, Coordinates *pco,
                    AthenaArray<Real> &a,
                    FaceField &b, Real time, Real dt,
                    int il, int iu, int jl, int ju, 
                    int kl, int ku, int ngh);

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
  assertCondition(mesh_bcs[BoundaryFace::inner_x3] == GetBoundaryFlag("user"),
                  "Use user-defined inner boundary condition along x3");
  assertCondition(mesh_bcs[BoundaryFace::outer_x3] == GetBoundaryFlag("user"),
                  "Use user-defined outer boundary condition along x3");

  // Read in constants
  rho0  = pin->GetReal("problem", "rho0");
  pgas0 = pin->GetReal("problem", "pgas0");
  vcir  = pin->GetReal("problem", "vcir");
  R0    = pin->GetReal("problem", "R0");

  // Enroll user-defined physical source terms
  EnrollUserExplicitSourceFunction(SourceFunctions);

  // Enroll user-defined boundary conditions
  EnrollUserBoundaryFunction(BoundaryFace::inner_x3, NoInflowInnerX3);
  EnrollUserBoundaryFunction(BoundaryFace::outer_x3, NoInflowOuterX3);

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

  Real cs2 = pgas0/rho0; // Isothermal sound speed
  Real H = std::sqrt(cs2)/(vcir/R0); // Scale height = cs/omega
  Real gamma = peos->GetGamma();

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {
        Real x = pcoord->x1v(i);
        Real y = pcoord->x2v(j);
        Real z = pcoord->x3v(k);

        // Hydrostatic Solution
        Real norm = rho0/std::exp(SQR(R0)/SQR(H));
        Real rho = norm*std::exp(SQR(R0/H)*
                                 std::pow(SQR(R0)/(SQR(R0)+SQR(z)),0.5));

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
// Zero gradient no inflow upper boundary condition
void NoInflowOuterX3(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &prim,
                         FaceField &b, Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        prim(IDN,ku+k,j,i) = prim(IDN,ku,j,i);
        prim(IPR,ku+k,j,i) = prim(IPR,ku,j,i);

        prim(IVX,ku+k,j,i) = prim(IVX,ku,j,i); 
        prim(IVY,ku+k,j,i) = prim(IVY,ku,j,i);
        prim(IVZ,ku+k,j,i) = std::max(0.,prim(IVZ,ku,j,i)); 
      }
    }
  }
  return;
}

// Zero gradient no inflow lower boundary condition
void NoInflowInnerX3(MeshBlock *pmb, Coordinates *pco,
                         AthenaArray<Real> &prim,
                         FaceField &b, Real time, Real dt,
                         int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=1; k<=ngh; k++) {
    for (int j=jl; j<=ju; j++) {
      for (int i=il; i<=iu; i++) {
        prim(IDN,kl-k,j,i) = prim(IDN,kl,j,i);
        prim(IPR,kl-k,j,i) = prim(IPR,kl,j,i);

        prim(IVX,kl-k,j,i) = prim(IVX,kl,j,i); 
        prim(IVY,kl-k,j,i) = prim(IVY,kl,j,i);
        prim(IVZ,kl-k,j,i) = std::min(0.,prim(IVZ,kl,j,i)); 
      }
    }
  }
  return;
}

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
  // Equivalent to g = z * GM/(R^2 + z^2)^1.5
  return SQR(vcir/R0)*std::pow(SQR(R0)/(SQR(R0)+SQR(z)),1.5)*z;
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
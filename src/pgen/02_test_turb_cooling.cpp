//========================================================================================
// Athena++ astrophysical MHD code
// Copyright(C) 2014 James M. Stone <jmstone@princeton.edu> and other code contributors
// Licensed under the 3-clause BSD License, see LICENSE file for details
//========================================================================================
//! \file turb.cpp
//! \brief Problem generator for turbulence driver

// C headers

// C++ headers
#include <cmath>
#include <ctime>
#include <sstream>
#include <stdexcept>

// Athena++ headers
#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../fft/athena_fft.hpp"
#include "../field/field.hpp"
#include "../globals.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../utils/utils.hpp"
#include "../cooling/cooling.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// User defined source functions
void cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim,
             const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc,
             AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar);

// User defined time step
static Real cooling_timestep(MeshBlock *pmb);

// Cooling
static Cooling cooler;

// Code Units (cgs)
static Real const amu       = 1.6605e-24;   // atomic mass unit (g)
static Real const kb        = 1.380648e-16; // boltzmann constant (erg/K)
static Real const mu        = 0.6173;       // mean molecular weight (Solar)
                                            // X = 0.7; Z = 0.02

static Real const unit_len  = 3.086e20;     // 100 pc
static Real const unit_temp = 1.0e4;        // 10^4 K
static Real const unit_n    = 1.0;          // cm^-3

static Real const unit_rho  = unit_n * amu * mu;
static Real const unit_pres = kb * unit_n * unit_temp;
static Real const unit_vel  = sqrt(unit_pres / unit_rho);
static Real const unit_time = unit_len / unit_vel;
static Real const unit_kap  = (unit_pres * unit_len * unit_len) /
                              (unit_time * unit_temp);
static Real const unit_gam  = (unit_pres/unit_time)/(unit_n*unit_n);

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
  cooler = Cooling(pin->GetString("cooling","cooling_table"));

  // Enroll user-defined physical source terms
  EnrollUserExplicitSourceFunction(cooling);

  // Enroll timestep so that dt <= min t_cool
  EnrollUserTimeStepFunction(cooling_timestep);

  // turb_flag is initialzed in the Mesh constructor to 0 by default;
  // turb_flag = 1 for decaying turbulence
  // turb_flag = 2 for impulsively driven turbulence
  // turb_flag = 3 for continuously driven turbulence
  turb_flag = pin->GetInteger("problem","turb_flag");
  if (turb_flag != 0) {
#ifndef FFT
    std::stringstream msg;
    msg << "### FATAL ERROR in TurbulenceDriver::TurbulenceDriver" << std::endl
        << "non zero Turbulence flag is set without FFT!" << std::endl;
    ATHENA_ERROR(msg);
    return;
#endif
  }

  return;
}

//========================================================================================
//! \fn void MeshBlock::ProblemGenerator(ParameterInput *pin)
//  \brief
//========================================================================================

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
  Real rho = pin->GetReal("problem","rho");
  Real temp = pin->GetReal("problem","temp");

  for (int k=ks; k<=ke; k++) {
    for (int j=js; j<=je; j++) {
      for (int i=is; i<=ie; i++) {

        phydro->u(IDN,k,j,i) = rho;

        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0; 

        phydro->u(IEN,k,j,i) = rho*temp/((5.0/3.0)-1.0);
      }
    }
  }
}


//========================================================================================
//! \fn void Mesh::UserWorkAfterLoop(ParameterInput *pin)
//  \brief
//========================================================================================

void Mesh::UserWorkAfterLoop(ParameterInput *pin) {
}

// User Defined Cooling Function
void cooling(MeshBlock *pmb, const Real time, const Real dt,
             const AthenaArray<Real> &prim,
             const AthenaArray<Real> &prim_scalar,
             const AthenaArray<Real> &bcc,
             AthenaArray<Real> &cons,
             AthenaArray<Real> &cons_scalar) {
  Real g = pmb->peos->GetGamma();

  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        // Need to take density and temperature at time step n from cons, not
        // prim because we do not need intermediate step to calculate the
        // cooling function
        Real rho = cons(IDN,k,j,i);
        Real eint = cons(IEN,k,j,i)
                    - 0.5 *(cons(IM1,k,j,i)*cons(IM1,k,j,i)
                          + cons(IM2,k,j,i)*cons(IM2,k,j,i)
                          + cons(IM3,k,j,i)*cons(IM3,k,j,i))/rho;        

        // T = P/rho
        Real temp = eint * (g-1.0)/rho;
        
        Real temp_new = temp;

        // Calculate new temperature using the Townsend Algorithm
        // The inputs should be given in cgs units
        Real temp_cgs = temp * unit_temp;
        Real rho_cgs  = rho  * unit_rho;
        Real dt_cgs   = dt   * unit_time;

        // if (i==2 && j==2 && k==2) std:: cout << "start : " << temp_new*unit_temp << "  " << rho_cgs << "  " << dt << std::endl; 
        temp_new = cooler.townsend_cooling(temp_cgs,rho_cgs,dt_cgs)/unit_temp;
        // if (i==2 && j==2 && k==2) std:: cout << "end : " << temp_new*unit_temp << std::endl; 

        // Update energy based on change in temperature
        cons(IEN,k,j,i) += (temp_new - temp) * (rho/(g-1.0));

        // Store the change in energy/time in a user defined output variable
        // pmb->user_out_var(0,k,j,i) = (temp_new - temp) * (rho/(g-1.0)) / dt;
 
      }
    }
  }

  return;
}

Real cooling_timestep(MeshBlock *pmb)
{
  Real min_dt=1.0e10;
  for (int k=pmb->ks; k<=pmb->ke; ++k) {
    for (int j=pmb->js; j<=pmb->je; ++j) {
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real T = pmb->phydro->w(IPR,k,j,i)/pmb->phydro->w(IDN,k,j,i);
        Real rho = pmb->phydro->w(IDN,k,j,i);
        Real tcool = cooler.single_point_cooling_time(T*unit_temp, rho*unit_rho)/unit_time;
 
        min_dt = std::min(min_dt, 0.25*tcool);
      }
    }
  }

  return min_dt;
}

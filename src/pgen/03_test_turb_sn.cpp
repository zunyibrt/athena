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
#include "../utils/sn_injection.hpp"

#ifdef OPENMP_PARALLEL
#include <omp.h>
#endif

// User defined source functions
void sn_injection(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc,
                  AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar);

// User defined history functions
Real CalculateSNEnergyInjection(MeshBlock *pmb, int iout);
Real CalculateSNMassInjection(MeshBlock *pmb, int iout);
Real CalculateSNVolumeInjection(MeshBlock *pmb, int iout);

// Code Units (cgs)
static Real const amu       = 1.6605e-24;   // atomic mass unit (g)
static Real const kb        = 1.380648e-16; // boltzmann constant (erg/K)
static Real const mu        = 0.6173;       // mean molecular weight (Solar)
                                            // X = 0.7; Z = 0.02

static Real const unit_len  = 3.086e20;     // 100 pc
static Real const unit_temp = 1.0e4;        // 10^4 K
static Real const unit_n    = 1.0;          // cm^-3

static Real const unit_rho  = unit_n * amu * mu;
static Real const unit_mass = unit_rho * std::pow(unit_len,3);
static Real const unit_pres = kb * unit_n * unit_temp;
static Real const unit_engy = unit_pres * std::pow(unit_len,3);
static Real const unit_vel  = std::sqrt(unit_pres / unit_rho);
static Real const unit_time = unit_len / unit_vel;
static Real const unit_kap  = (unit_pres * unit_len * unit_len) /
                              (unit_time * unit_temp);
static Real const unit_gam  = (unit_pres/unit_time)/(unit_n*unit_n);

// Globals
static std::vector<Real> sn_times;
static int next_sn_idx;
static Real r_inj,e_sn,m_ej;

//========================================================================================
//! \fn void Mesh::InitUserMeshData(ParameterInput *pin)
//  \brief
//========================================================================================
void Mesh::InitUserMeshData(ParameterInput *pin) {
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

  EnrollUserExplicitSourceFunction(sn_injection);
 
  // Generate SN times with M_cluster and tlim
  next_sn_idx = 0;
  Real cluster_mass = pin->GetOrAddReal("SN","M_cluster",1e6);
  sn_times = get_sn_times(cluster_mass,pin->GetReal("SN","tlim"));

  // Compute E and M injection densities from injection radius
  r_inj = pin->GetReal("SN","r_inj"); // Input in code unis
  Real sphere_vol = (4.0/3.0)*PI*std::pow(r_inj,3);
  e_sn  = pin->GetOrAddReal("SN","E_sn",std::pow(10,51)/unit_engy)/sphere_vol; // Input in ergs
  m_ej  = pin->GetOrAddReal("SN","M_ej",8.4*1.988e33/unit_mass)/sphere_vol; // Input in solar mass

  std::cout << M_ej

  // Enroll user defined history outputs
  AllocateUserHistoryOutput(3);
  EnrollUserHistoryOutput(0,CalculateSNEnergyInjection,"SNEnergyInjection");
  EnrollUserHistoryOutput(1,CalculateSNMassInjection,"SNMassInjection");
  EnrollUserHistoryOutput(2,CalculateSNVolumeInjection,"SNVolumeInjection");

  return;
}

void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateUserOutputVariables(3);
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

// User Defined SN Function
void sn_injection(MeshBlock *pmb, const Real time, const Real dt,
                  const AthenaArray<Real> &prim,
                  const AthenaArray<Real> &prim_scalar,
                  const AthenaArray<Real> &bcc,
                  AthenaArray<Real> &cons,
                  AthenaArray<Real> &cons_scalar) {
  if (time > sn_times.at(next_sn_idx)) {
    for (int k=pmb->ks; k<=pmb->ke; ++k) {
      for (int j=pmb->js; j<=pmb->je; ++j) {
        for (int i=pmb->is; i<=pmb->ie; ++i) {
          Real x = pmb->pcoord->x1v(i);
          Real y = pmb->pcoord->x2v(j);
          Real z = pmb->pcoord->x3v(k);
          
          if (SQR(x) + SQR(y) + SQR(z) < SQR(r_inj)) {
            // Inject energy and mass into cell
            cons(IEN,k,j,i) += e_sn;
            cons(IDN,k,j,i) += m_ej;

            // Store the changes in a user defined output variable
            pmb->user_out_var(0,k,j,i) += e_sn;
            pmb->user_out_var(1,k,j,i) += m_ej;
            pmb->user_out_var(2,k,j,i) += 1.0;
          }        
        }
      }
    }
    next_sn_idx++;
  }

  return;
}

// User defined history functions
Real CalculateSNEnergyInjection(MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> vol(pmb->ncells1);

  Real sum = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    for (int i=is; i<=ie; i++) {
      sum += pmb->user_out_var(0,k,j,i)*vol(i);
    }
  }}

  return sum;
}

Real CalculateSNMassInjection(MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> vol(pmb->ncells1);

  Real sum = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    for (int i=is; i<=ie; i++) {
      sum += pmb->user_out_var(1,k,j,i)*vol(i);
    }
  }}

  return sum;
}

Real CalculateSNVolumeInjection(MeshBlock *pmb, int iout) {
  int is=pmb->is, ie=pmb->ie, js=pmb->js, je=pmb->je, ks=pmb->ks, ke=pmb->ke;
  AthenaArray<Real> vol(pmb->ncells1);

  Real sum = 0.0;

  for (int k=ks; k<=ke; k++) {
  for (int j=js; j<=je; j++) {
    pmb->pcoord->CellVolume(k,j,is,ie,vol);
    for (int i=is; i<=ie; i++) {
      sum += pmb->user_out_var(2,k,j,i)*vol(i);
    }
  }}

  return sum;
}
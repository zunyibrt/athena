// Local Lax-Friedrichs Riemann solver for relativistic magnetohydrodynamics

// Primary header
#include "../../hydro_integrator.hpp"

// C++ headers
#include <algorithm>  // max(), min()
#include <cmath>      // sqrt()

// Athena headers
#include "../../../hydro.hpp"                       // Hydro
#include "../../../eos/eos.hpp"                     // HydroEqnOfState
#include "../../../../athena.hpp"                   // enums, macros, Real
#include "../../../../athena_arrays.hpp"            // AthenaArray
#include "../../../../mesh.hpp"                     // MeshBlock
#include "../../../../coordinates/coordinates.hpp"  // Coordinates

// Riemann solver
// Inputs:
//   k,j: x3- and x2-indices
//   il,iu: lower and upper x1-indices
//   ivx: type of interface (IVX for x1, IVY for x2, IVZ for x3)
//   bb: 3D array of normal magnetic fields
//   prim_l, prim_r: left and right primitive states
// Outputs:
//   flux: fluxes across interface
// Notes:
//   prim_l, prim_r overwritten
//   implements LLF algorithm similar to that of fluxcalc() in step_ch.c in Harm
//   references Mignone & Bodo 2006, MNRAS 368 1040 (MB)
//   references Mignone, Ugliano, & Bodo 2009, MNRAS 393 1141 (MUB)
void HydroIntegrator::RiemannSolver(const int k, const int j, const int il,
    const int iu, const int ivx, const AthenaArray<Real> &bb, AthenaArray<Real> &prim_l,
    AthenaArray<Real> &prim_r, AthenaArray<Real> &flux)
{
  // Transform primitives to locally flat coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_hydro->pmy_block->pcoord->PrimToLocal1(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVY:
        pmy_hydro->pmy_block->pcoord->PrimToLocal2(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
      case IVZ:
        pmy_hydro->pmy_block->pcoord->PrimToLocal3(k, j, il, iu, bb, prim_l, prim_r,
            bb_normal_);
        break;
    }
  else  // SR; need to populate 1D normal B array
  {
    #pragma simd
    for (int i = il; i <= iu; ++i)
      bb_normal_(i) = bb(k,j,i);
  }

  // Calculate wavespeeds
  pmy_hydro->pf_eos->FastMagnetosonicSpeedsSR(prim_l, bb_normal_, il, iu, ivx,
      lambdas_p_l_, lambdas_m_l_);
  pmy_hydro->pf_eos->FastMagnetosonicSpeedsSR(prim_r, bb_normal_, il, iu, ivx,
      lambdas_p_r_, lambdas_m_r_);

  // Calculate cyclic permutations of indices
  int ivy = IVX + ((ivx-IVX)+1)%3;
  int ivz = IVX + ((ivx-IVX)+2)%3;

  // Extract ratio of specific heats
  const Real gamma_adi = pmy_hydro->pf_eos->GetGamma();

  // Go through each interface
  #pragma simd
  for (int i = il; i <= iu; ++i)
  {
    // Extract left primitives
    const Real &rho_l = prim_l(IDN,i);
    const Real &pgas_l = prim_l(IEN,i);
    Real u_l[4];
    if (GENERAL_RELATIVITY)
    {
      u_l[1] = prim_l(ivx,i);
      u_l[2] = prim_l(ivy,i);
      u_l[3] = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 + SQR(u_l[1]) + SQR(u_l[2]) + SQR(u_l[3]));
    }
    else  // SR
    {
      const Real &vx_l = prim_l(ivx,i);
      const Real &vy_l = prim_l(ivy,i);
      const Real &vz_l = prim_l(ivz,i);
      u_l[0] = std::sqrt(1.0 / (1.0 - SQR(vx_l) - SQR(vy_l) - SQR(vz_l)));
      u_l[1] = u_l[0] * vx_l;
      u_l[2] = u_l[0] * vy_l;
      u_l[3] = u_l[0] * vz_l;
    }
    const Real &bb2_l = prim_l(IBY,i);
    const Real &bb3_l = prim_l(IBZ,i);

    // Extract right primitives
    const Real &rho_r = prim_r(IDN,i);
    const Real &pgas_r = prim_r(IEN,i);
    Real u_r[4];
    if (GENERAL_RELATIVITY)
    {
      u_r[1] = prim_r(ivx,i);
      u_r[2] = prim_r(ivy,i);
      u_r[3] = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 + SQR(u_r[1]) + SQR(u_r[2]) + SQR(u_r[3]));
    }
    else  // SR
    {
      const Real &vx_r = prim_r(ivx,i);
      const Real &vy_r = prim_r(ivy,i);
      const Real &vz_r = prim_r(ivz,i);
      u_r[0] = std::sqrt(1.0 / (1.0 - SQR(vx_r) - SQR(vy_r) - SQR(vz_r)));
      u_r[1] = u_r[0] * vx_r;
      u_r[2] = u_r[0] * vy_r;
      u_r[3] = u_r[0] * vz_r;
    }
    const Real &bb2_r = prim_r(IBY,i);
    const Real &bb3_r = prim_r(IBZ,i);

    // Extract normal magnetic field
    const Real &bb1 = bb_normal_(i);

    // Calculate 4-magnetic field in left state
    Real b_l[4];
    b_l[0] = bb1*u_l[1] + bb2_l*u_l[2] + bb3_l*u_l[3];
    b_l[1] = (bb1 + b_l[0] * u_l[1]) / u_l[0];
    b_l[2] = (bb2_l + b_l[0] * u_l[2]) / u_l[0];
    b_l[3] = (bb3_l + b_l[0] * u_l[3]) / u_l[0];
    Real b_sq_l = -SQR(b_l[0]) + SQR(b_l[1]) + SQR(b_l[2]) + SQR(b_l[3]);

    // Calculate 4-magnetic field in right state
    Real b_r[4];
    b_r[0] = bb1*u_r[1] + bb2_r*u_r[2] + bb3_r*u_r[3];
    b_r[1] = (bb1 + b_r[0] * u_r[1]) / u_r[0];
    b_r[2] = (bb2_r + b_r[0] * u_r[2]) / u_r[0];
    b_r[3] = (bb3_r + b_r[0] * u_r[3]) / u_r[0];
    Real b_sq_r = -SQR(b_r[0]) + SQR(b_r[1]) + SQR(b_r[2]) + SQR(b_r[3]);

    // Calculate extremal wavespeeds
    Real lambda_l = std::min(lambdas_m_l_(i), lambdas_m_r_(i));  // (MB 55)
    Real lambda_r = std::max(lambdas_p_l_(i), lambdas_p_r_(i));  // (MB 55)
    Real lambda = std::max(lambda_r, -lambda_l);

    // Calculate conserved quantities in L region (MUB 8)
    Real cons_l[NWAVE];
    Real wtot_l = rho_l + gamma_adi/(gamma_adi-1.0) * pgas_l + b_sq_l;
    Real ptot_l = pgas_l + 0.5*b_sq_l;
    cons_l[IDN] = rho_l * u_l[0];
    cons_l[IEN] = wtot_l * u_l[0] * u_l[0] - b_l[0] * b_l[0] - ptot_l;
    cons_l[ivx] = wtot_l * u_l[1] * u_l[0] - b_l[1] * b_l[0];
    cons_l[ivy] = wtot_l * u_l[2] * u_l[0] - b_l[2] * b_l[0];
    cons_l[ivz] = wtot_l * u_l[3] * u_l[0] - b_l[3] * b_l[0];
    cons_l[IBY] = b_l[2] * u_l[0] - b_l[0] * u_l[2];
    cons_l[IBZ] = b_l[3] * u_l[0] - b_l[0] * u_l[3];

    // Calculate fluxes in L region (MUB 15)
    Real flux_l[NWAVE];
    flux_l[IDN] = rho_l * u_l[1];
    flux_l[IEN] = wtot_l * u_l[0] * u_l[1] - b_l[0] * b_l[1];
    flux_l[ivx] = wtot_l * u_l[1] * u_l[1] - b_l[1] * b_l[1] + ptot_l;
    flux_l[ivy] = wtot_l * u_l[2] * u_l[1] - b_l[2] * b_l[1];
    flux_l[ivz] = wtot_l * u_l[3] * u_l[1] - b_l[3] * b_l[1];
    flux_l[IBY] = b_l[2] * u_l[1] - b_l[1] * u_l[2];
    flux_l[IBZ] = b_l[3] * u_l[1] - b_l[1] * u_l[3];

    // Calculate conserved quantities in R region (MUB 8)
    Real cons_r[NWAVE];
    Real wtot_r = rho_r + gamma_adi/(gamma_adi-1.0) * pgas_r + b_sq_r;
    Real ptot_r = pgas_r + 0.5*b_sq_r;
    cons_r[IDN] = rho_r * u_r[0];
    cons_r[IEN] = wtot_r * u_r[0] * u_r[0] - b_r[0] * b_r[0] - ptot_r;
    cons_r[ivx] = wtot_r * u_r[1] * u_r[0] - b_r[1] * b_r[0];
    cons_r[ivy] = wtot_r * u_r[2] * u_r[0] - b_r[2] * b_r[0];
    cons_r[ivz] = wtot_r * u_r[3] * u_r[0] - b_r[3] * b_r[0];
    cons_r[IBY] = b_r[2] * u_r[0] - b_r[0] * u_r[2];
    cons_r[IBZ] = b_r[3] * u_r[0] - b_r[0] * u_r[3];

    // Calculate fluxes in R region (MUB 15)
    Real flux_r[NWAVE];
    flux_r[IDN] = rho_r * u_r[1];
    flux_r[IEN] = wtot_r * u_r[0] * u_r[1] - b_r[0] * b_r[1];
    flux_r[ivx] = wtot_r * u_r[1] * u_r[1] - b_r[1] * b_r[1] + ptot_r;
    flux_r[ivy] = wtot_r * u_r[2] * u_r[1] - b_r[2] * b_r[1];
    flux_r[ivz] = wtot_r * u_r[3] * u_r[1] - b_r[3] * b_r[1];
    flux_r[IBY] = b_r[2] * u_r[1] - b_r[1] * u_r[2];
    flux_r[IBZ] = b_r[3] * u_r[1] - b_r[1] * u_r[3];

    // Set fluxes
    for (int n = 0; n < NWAVE; ++n)
      flux(n,i) = 0.5 * (flux_l[n] + flux_r[n] - lambda * (cons_r[n] - cons_l[n]));

    // Set conserved quantities in GR
    if (GENERAL_RELATIVITY)
      for (int n = 0; n < NWAVE; ++n)
        cons_(n,i) = 0.5 * (cons_r[n] + cons_l[n] + (flux_l[n] - flux_r[n]) / lambda);
  }

  // Transform fluxes to global coordinates if in GR
  if (GENERAL_RELATIVITY)
    switch (ivx)
    {
      case IVX:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal1(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVY:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal2(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
      case IVZ:
        pmy_hydro->pmy_block->pcoord->FluxToGlobal3(k, j, il, iu, cons_, bb_normal_,
            flux);
        break;
    }
  return;
}
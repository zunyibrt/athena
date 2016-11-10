//======================================================================================
// Athena++ astrophysical MHD code
// Copyright (C) 2014 James M. Stone  <jmstone@princeton.edu>
//
// This program is free software: you can redistribute and/or modify it under the terms
// of the GNU General Public License (GPL) as published by the Free Software Foundation,
// either version 3 of the License, or (at your option) any later version.
//
// This program is distributed in the hope that it will be useful, but WITHOUT ANY
// WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A 
// PARTICULAR PURPOSE.  See the GNU General Public License for more details.
//
// You should have received a copy of GNU GPL in the file LICENSE included in the code
// distribution.  If not see <http://www.gnu.org/licenses/>.
//======================================================================================

// Athena headers
#include "../athena.hpp"

#include <iostream>


//======================================================================================
//! \file FourthPolyRoot.cpp
//======================================================================================

//--------------------------------------------------------------------------------------


// Exact solution for fourth order polynomical with the format
// coef4 * x^4 + x + tconst == 0

int FouthPolyRoot(const Real coef4, const Real tconst, Real &root)
{

// First, get the real root of
// z^3-4*tconst/coef4 * z - 1/coef4^2==0
  Real asquar = coef4 * coef4;
  Real acubic = coef4 * asquar;
  Real ccubic = tconst * tconst * tconst;
  Real delta1 = 0.25 - 64.0 * ccubic * coef4/27.0;
  if(delta1 < 0.0) return -1;
  else delta1 = sqrt(delta1);
  Real zroot = 0.0;
  if(delta1 > 1.e11){
    // to avoid small number cancellation
    zroot = pow(delta1,-2.0/3.0)/3.0;
  }else{
    zroot = (pow(0.5 + delta1, 1.0/3.0) -
               pow(-0.5 + delta1, 1.0/3.0));
  }

  zroot *= pow(coef4,-2.0/3.0);
  
  Real rcoef = sqrt(zroot);
  Real delta2 = -zroot + 2.0/(coef4*rcoef);
  if(delta2 < 0.0) return -1;
  else delta2 = sqrt(delta2);
  root = 0.5 * (delta2 - rcoef);
  if(root < 0.0) return -1;

  return 0;
}
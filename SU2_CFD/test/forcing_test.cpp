/*!
 * \file forcing_test.cpp
 * \brief Test the forcing terms for the hybrid RANS/LES model
 * \author C. Pederson
 * \version 5.0.0 "Raven"
 *
 * SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
 *                      Dr. Thomas D. Economon (economon@stanford.edu).
 *
 * SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
 *                 Prof. Piero Colonna's group at Delft University of Technology.
 *                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *                 Prof. Rafael Palacios' group at Imperial College London.
 *                 Prof. Edwin van der Weide's group at the University of Twente.
 *                 Prof. Vincent Terrapon's group at the University of Liege.
 *
 * Copyright (C) 2012-2016 SU2, the open-source CFD code.
 *
 * SU2 is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public
 * License as published by the Free Software Foundation; either
 * version 2.1 of the License, or (at your option) any later version.
 *
 * SU2 is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with SU2. If not, see <http://www.gnu.org/licenses/>.
 */

#ifdef BUILD_TESTS

#define BOOST_TEST_MODULE Hybrid_Model
#include "boost/test/included/unit_test.hpp"

#include "../include/hybrid_RANS_LES_forcing.hpp"

BOOST_AUTO_TEST_CASE(Coordinate_Shifting) {

  const unsigned short nDim = 3;
  const unsigned short nVar = nDim + 2;
  const su2double pi = std::atan(1.0)*4.0;

  CHybridForcing forcing(nDim);

  su2double time = 1.0;
  su2double L_m = 1.0;
  su2double T_m = 2.0;
  su2double alpha = 1.0;
  su2double x_original[nDim], x_shifted[nDim], x_exact[nDim];
  su2double prim_vars[nVar];

  x_original[0] = -1.0; x_original[1] = 0; x_original[2] = 1.0;
  x_exact[0] = 0; x_exact[1] = 2*pi; x_exact[2] = 4*pi;

  prim_vars[1] = 1.0; prim_vars[2] = 1.0; prim_vars[3] = 1.0;

  forcing.SetTurbLengthscale(L_m);
  forcing.SetTurbTimescale(T_m);
  forcing.SetHybridParameter(&alpha);
  forcing.SetPrimitive(prim_vars);

  su2double tol = 1e-8;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    x_shifted[iDim] = forcing.TransformCoords(x_original, iDim, time);
    BOOST_CHECK_SMALL(x_shifted[iDim] - x_exact[iDim], tol);
  }

}

BOOST_AUTO_TEST_CASE(Forcing_Functions) {

  const unsigned short nDim = 3;
  const unsigned short nVar = nDim + 2;
  const su2double pi = std::atan(1.0)*4.0;

  /*--- Setup ---*/

  CHybridForcing forcing(nDim);

  su2double x[nDim];
  su2double prim_vars[nVar];
  x[0] = pi/6.0; x[1] = pi/3.0; x[2] = pi/2.0;
  prim_vars[1] = 3.0; prim_vars[2] = 2.0; prim_vars[3] = 1.0;

  su2double** tau_F = new su2double*[nDim];
  su2double** tau_F_exact = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    tau_F[iDim] = new su2double[nDim];
    tau_F_exact[iDim] = new su2double[nDim];
  }

  /*--- Hand-computed test results ---*/
  tau_F_exact[0][0] = 5.0625;
  tau_F_exact[0][1] = 2.4375;
  tau_F_exact[0][2] = 0.75;
  tau_F_exact[1][0] = 2.4375;
  tau_F_exact[1][1] = 1.0625;
  tau_F_exact[1][2] = 0.25;
  tau_F_exact[2][0] = 0.75;
  tau_F_exact[2][1] = 0.25;
  tau_F_exact[2][2] = 0.0;

  /*--- Test ---*/

  forcing.SetPrimitive(prim_vars);
  forcing.BuildForcingStress(x, tau_F);

  su2double tol = 1e-8;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      BOOST_CHECK_SMALL(tau_F[iDim][iDim] - tau_F_exact[iDim][iDim], tol);
    }
  }


  /*--- Teardown ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] tau_F_exact[iDim];
  }
  delete [] tau_F_exact;

}

BOOST_AUTO_TEST_CASE(Smoke_Test_for_Forcing) {

  /*--- Setup ---*/

  const unsigned short nDim = 3;
  const unsigned short nVar = nDim + 2;

  su2double alpha = 1.0;

  su2double x[nDim];
  su2double prim_vars[nVar];
  su2double turb_vars[4];
  su2double source_terms[2];
  x[0] = 0.0; x[1] = 0.0; x[2] = 0.0;
  prim_vars[1] = 1.0; prim_vars[2] = 1.0; prim_vars[3] = 1.0;
  turb_vars[0] = 1.0; turb_vars[1] = 1.0;
  source_terms[0] = 1.0; source_terms[1] = 0.0;

  su2double** prim_var_grad = new su2double*[nVar];
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    prim_var_grad[iVar] = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      prim_var_grad[iVar][iDim] = 0.0;
  }

  /*--- Calculate forcing ---*/

  CHybridForcing forcing(nDim);

  forcing.SetCoord(x);
  forcing.SetTurbLengthscale(1.0);
  forcing.SetTurbTimescale(1.0);
  forcing.SetHybridParameter(&alpha);
  forcing.SetPrimitive(prim_vars);
  forcing.SetPrimVarGradient(prim_var_grad);
  forcing.SetTurbVar(turb_vars);
  forcing.SetHybridSourceTerms(source_terms);

  su2double** tau_F = forcing.GetStress();
  su2double P_F = forcing.GetProduction();
  su2double P_F_over_P_lim = forcing.GetForcingRatio();

  /*--- Teardown ---*/

  for (unsigned short iVar = 0; iVar < nVar; iVar++)
    delete [] prim_var_grad[iVar];
  delete [] prim_var_grad;

}

#endif

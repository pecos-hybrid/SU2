/*!
 * \file resolution_integration_test.cpp
 * \brief This test checks whether the resolution tensor is correctly set for a grid
 * of quadrilateral cells.
 * \author C. Pederson
 * \version 4.3.0 "Cardinal"
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

#define BOOST_TEST_MODULE ViscousProjFlux
#include "MPI_global_fixture.hpp"


#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

static const unsigned short nDim = 3;
static const unsigned short nVar = 5;

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

BOOST_AUTO_TEST_CASE(IsotropicRANSFlux) {

  /*---
   * SETUP
   * ---*/

  CConfig* test_config = new CConfig();

  CAvgGrad_Flow numerics(nDim, nVar, false, test_config);

  su2double** gradprimvar = new su2double*[nVar];
  for (int iVar = 0; iVar < nVar; iVar++) {
    gradprimvar[iVar] = new su2double[nDim];
    for (int jDim = 0; jDim < nDim; jDim++)
      gradprimvar[iVar][jDim] = 0.0;
  }

  gradprimvar[1][0] = -1; // dU/dx
  gradprimvar[2][1] = 1; // dV/dy
  su2double normal[nDim] = {1.0, 0.0, 0.0};
  su2double prim_var[6] = {1.0, 1.0, 1.0, 1.0, 1.0, 1.0};
  const su2double turb_ke = 3.0;
  const su2double laminar_viscosity = 0.0;
  const su2double eddy_viscosity = 1.0;

  /*---
   * TEST
   * ---*/

//  numerics.GetTau(prim_var, gradprimvar, turb_ke,
//                  laminar_viscosity, eddy_viscosity, false, tau);
//  numerics.GetViscousHeatFlux(gradprimvar, laminar_viscosity, eddy_viscosity);
//  numerics.GetViscousProjFlux(prim_var, gradprimvar, normal);
//
//  const su2double* output = numerics.Proj_Flux_Tensor;
//  const su2double correct_output[nVar] = {0.0, -4.0, 0.0, 0.0, -4.0};
//
//  const su2double machine_eps = std::numeric_limits<su2double>::epsilon();
//  for (int iVar = 0; iVar < nVar; iVar++) {
//    // Build the test info
//    std::stringstream msg;
//    msg << "iVar = " << iVar << "\n";
//    BOOST_TEST_CONTEXT(msg.str());
//    BOOST_CHECK_SMALL(output[iVar] - correct_output[iVar], machine_eps);
//  }

  /*---
   * TEARDOWN
   * ---*/

  delete test_config;
  for (int iVar = 0; iVar < nVar; iVar++) {
    delete [] gradprimvar[iVar];
  }
  delete [] gradprimvar;

}

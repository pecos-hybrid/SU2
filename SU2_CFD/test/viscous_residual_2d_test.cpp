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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

namespace viscous_residual_2d_test {

const unsigned short nDim = 2;
const unsigned short nVar = 4;
const unsigned short nPrimVar = nVar+5;

static void PrintInformation(su2double* residual_i,
                      su2double** Jacobian_i,
                      su2double** Jacobian_j) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    cout << "residual_i[" << iVar << "]: ";
    cout << residual_i[iVar] << "\n";
  }
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      cout << "expected_jacobian_i[" << iVar << "][" << jVar << "] = ";
      cout << Jacobian_i[iVar][jVar] << "\n";
    }
  }
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      cout << "Jacobian_j[" << iVar << "][" << jVar << "]: ";
      cout << Jacobian_j[iVar][jVar] << "\n";
    }
  }
}

static void WriteCfgFile(const char* filename) {

  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "TIME_DISCRE_FLOW= EULER_IMPLICIT" << std::endl;

  cfg_file.close();

}

struct ViscousResidualFixture{
  ViscousResidualFixture()
      : distance(1), area(3) {

    char cfg_filename[100] = "viscous_residual_2d_test.cfg";
    WriteCfgFile(cfg_filename);
    const unsigned short iZone = 0;
    const unsigned short nZone = 1;
    config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
    config->SetGas_ConstantND(287.058);
    std::remove(cfg_filename);

    numerics = new CAvgGrad_Flow(2, 4, false, config);

    /*--- Inputs ---*/

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      coord_i[iDim] = 0;
      coord_j[iDim] = 0;
    }
    coord_j[0] += distance;

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      normal[iDim] = coord_j[iDim]/distance*area;
    }

    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
      primvar_i[iVar] = 0.0; primvar_j[iVar] = 0.0;
    }

    primvar_grad_i = new su2double*[nPrimVar];
    primvar_grad_j = new su2double*[nPrimVar];
    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
      primvar_grad_i[iVar] = new su2double[nDim];
      primvar_grad_j[iVar] = new su2double[nDim];
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        primvar_grad_i[iVar][iDim] = 0.0;
        primvar_grad_j[iVar][iDim] = 0.0;
      }
    }

    turbvar_grad = new su2double*[1];
    for (unsigned short iVar = 0; iVar < 1; iVar++) {
      turbvar_grad[iVar] = new su2double[nDim];
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        turbvar_grad[iVar][iDim] = 0.0;
      }
    }

    numerics->SetCoord(coord_i, coord_j);
    numerics->SetNormal(normal);
    numerics->SetSecondary(NULL, NULL);

    /*--- Outputs ---*/

    Jacobian_i = new su2double*[nVar];
    Jacobian_j = new su2double*[nVar];
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double[nVar];
      Jacobian_j[iVar] = new su2double[nVar];
    }

    residual_i = new su2double[nVar];

    expected_jacobian_i = new su2double*[nVar];
    expected_jacobian_j = new su2double*[nVar];
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      expected_jacobian_i[iVar] = new su2double[nVar];
      expected_jacobian_j[iVar] = new su2double[nVar];
      for (unsigned short jVar = 0; jVar < nVar; jVar++) {
        expected_jacobian_i[iVar][jVar] = 0;
        expected_jacobian_j[iVar][jVar] = 0;
      }
    }
  }

  ~ViscousResidualFixture() {
    delete config;
    delete numerics;
    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] primvar_grad_i[iVar];
      delete[] primvar_grad_j[iVar];
    }
    delete[] primvar_grad_i;
    delete[] primvar_grad_j;

    for (unsigned short iVar = 0; iVar < 1; iVar++) {
      delete[] turbvar_grad[iVar];
    }
    delete[] turbvar_grad;

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] Jacobian_i[iVar];
      delete[] Jacobian_j[iVar];
    }
    delete[] Jacobian_i;
    delete[] Jacobian_j;
    delete[] residual_i;

    for (unsigned short iVar = 0; iVar < nVar; iVar++) {
      delete[] expected_jacobian_i[iVar];
      delete[] expected_jacobian_j[iVar];
    }
    delete[] expected_jacobian_i;
    delete[] expected_jacobian_j;
  }

  CConfig* config;
  CNumerics* numerics;
  const su2double distance;
  const su2double area;
  su2double coord_i[nDim], coord_j[nDim];
  su2double normal[nDim];
  su2double** primvar_grad_i, **primvar_grad_j;
  su2double primvar_i[nPrimVar];
  su2double primvar_j[nPrimVar];
  su2double** turbvar_grad;
  su2double** Jacobian_i, **Jacobian_j;
  su2double* residual_i;

  su2double** expected_jacobian_i, **expected_jacobian_j;
};

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE(ViscousResidual2dTest);

BOOST_FIXTURE_TEST_CASE(ViscousResidualwithEverything, ViscousResidualFixture) {

  /*---
   * SETUP
   * ---*/

  primvar_i[nDim+1] = 1.0; // pressure
  primvar_i[nDim+2] = 1.0; // density
  primvar_i[nDim+5] = 1.0; // laminar viscosity
  primvar_i[nDim+6] = 1.0; // turbulent viscosity
  for (unsigned short iVar = 1; iVar < nDim+1; iVar++) {
    primvar_i[iVar] = iVar; // Velocities
  }
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_j[iVar] = primvar_i[iVar];
  }

  primvar_grad_i[1][0] =  1.0; // du/dx
  primvar_grad_i[2][1] =  2.0; // dv/dy
  primvar_grad_i[1][1] =  1.0; // du/dy
  primvar_grad_i[2][0] =  1.0; // dv/dx
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      primvar_grad_j[iVar][iDim] = primvar_grad_i[iVar][iDim];
    }
  }

  const su2double tke = 3; // 3 cancels out 3 in denominator

  /*---
   * TEST
   * ---*/

  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  numerics->SetTurbVarGradient(turbvar_grad, turbvar_grad);
  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

  su2double expected_residual[nVar] = {0, -6, 12, 18};
  expected_jacobian_i[1][0] = 8;
  expected_jacobian_i[1][1] = -8;
  expected_jacobian_i[2][0] = 12;
  expected_jacobian_i[2][2] = -6;
  expected_jacobian_i[3][3] = -10.5;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      expected_jacobian_j[iVar][jVar] = -expected_jacobian_i[iVar][jVar];
    }
  }
  expected_jacobian_i[3][0] = 23;
  expected_jacobian_i[3][1] = -0.5;
  expected_jacobian_i[3][2] = 15;

  expected_jacobian_j[3][0] = -41;
  expected_jacobian_j[3][1] = -5.5;
  expected_jacobian_j[3][2] = -3;

  cout << "Calculated:\n";
  PrintInformation(residual_i, Jacobian_i, Jacobian_j);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    BOOST_CHECK_CLOSE_FRACTION(expected_residual[iVar], residual_i[iVar], tolerance);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      cout << "iVar: " << iVar << "\tjVar: " << jVar << "\n";
      BOOST_CHECK_CLOSE_FRACTION(expected_jacobian_i[iVar][jVar], Jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_CLOSE_FRACTION(expected_jacobian_j[iVar][jVar], Jacobian_j[iVar][jVar], tolerance);
    }
  }
}

BOOST_AUTO_TEST_SUITE_END();

} // end namespace

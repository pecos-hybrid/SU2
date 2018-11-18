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

#define BOOST_TEST_MODULE ViscousResidual
#include "MPI_global_fixture.hpp"


#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

#ifdef BUILD_TESTS

const unsigned short nDim = 3;
const unsigned short nVar = 5;
const unsigned short nPrimVar = nVar+5;

void PrintInformation(su2double* residual_i,
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

class TestConfig : public CConfig {
 public:
  TestConfig() : CConfig() {};
  unsigned short GetKind_TimeIntScheme_Flow(void) { return EULER_IMPLICIT; };
};

struct ViscousResidualFixture{
  ViscousResidualFixture()
      : distance(1), area(3) {

    config = new TestConfig();
    numerics = new CAvgGrad_Flow(3, 5, config);

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
  su2double** Jacobian_i, **Jacobian_j;
  su2double* residual_i;

  su2double** expected_jacobian_i, **expected_jacobian_j;
};

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

BOOST_FIXTURE_TEST_CASE(ViscousResidualwithRotationOnly, ViscousResidualFixture) {

  /*---
   * SETUP
   * ---*/

  primvar_i[nDim+1] = 1.0; // pressure
  primvar_i[nDim+2] = 1.0; // density
  primvar_i[nDim+5] = 1.0; // laminar viscosity
  primvar_i[nDim+6] = 1.0; // turbulent viscosity
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_j[iVar] = primvar_i[iVar];
  }

  primvar_grad_i[1][1] =  1.0; // du/dy
  primvar_grad_i[2][0] = -1.0; // dv/dx
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      primvar_grad_j[iVar][iDim] = primvar_grad_i[iVar][iDim];
    }
  }

  const su2double tke = 0;

  /*---
   * TEST
   * ---*/

  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

  su2double expected_residual[nVar] = {0, 0, 0, 0, 0};
  expected_jacobian_i[1][1] = -8;
  expected_jacobian_i[2][2] = -6;
  expected_jacobian_i[3][3] = -6;
  expected_jacobian_i[4][0] = 26.25;
  expected_jacobian_i[4][4] = -10.5;

  // cout << "Calculated:\n";
  // PrintInformation(residual_i, Jacobian_i, Jacobian_j);

  const su2double tolerance = 100*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    BOOST_CHECK_SMALL(expected_residual[iVar] - residual_i[iVar], tolerance);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      // cout << "iVar: " << iVar << "\tjVar: " << jVar << "\n";
      BOOST_CHECK_SMALL(expected_jacobian_i[iVar][jVar] - Jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_SMALL(expected_jacobian_i[iVar][jVar] + Jacobian_j[iVar][jVar], tolerance);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(ViscousResidualwithNoViscosity, ViscousResidualFixture) {

  /*---
   * SETUP
   * ---*/

  primvar_i[nDim+1] = 1.0; // pressure
  primvar_i[nDim+2] = 1.0; // density
  primvar_i[nDim+5] = 0.0; // laminar viscosity
  primvar_i[nDim+6] = 0.0; // turbulent viscosity
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_j[iVar] = primvar_i[iVar];
  }

  primvar_grad_i[1][1] =  1.0; // du/dy
  primvar_grad_i[2][0] =  1.0; // dv/dx
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      primvar_grad_j[iVar][iDim] = primvar_grad_i[iVar][iDim];
    }
  }

  const su2double tke = 0;

  /*---
   * TEST
   * ---*/

  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

  su2double expected_residual[nVar] = {0, 0, 0, 0, 0};

  // cout << "Calculated:\n";
  // PrintInformation(residual_i, Jacobian_i, Jacobian_j);

  const su2double tolerance = 100*std::numeric_limits<su2double>::epsilon();
  // Ignore Jacobian of energy flux (zero viscosity can create NaNs)
  for (unsigned short iVar = 0; iVar < nVar-1; iVar++) {
    BOOST_CHECK_SMALL(expected_residual[iVar] - residual_i[iVar], tolerance);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      // cout << "iVar: " << iVar << "\tjVar: " << jVar << "\n";
      BOOST_CHECK_SMALL(expected_jacobian_i[iVar][jVar] - Jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_SMALL(expected_jacobian_i[iVar][jVar] + Jacobian_j[iVar][jVar], tolerance);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(ViscousResidualwithTKEOnly, ViscousResidualFixture) {

  /*---
   * SETUP
   * ---*/

  primvar_i[nDim+1] = 1.0; // pressure
  primvar_i[nDim+2] = 1.0; // density
  primvar_i[nDim+5] = 1.0; // laminar viscosity
  primvar_i[nDim+6] = 1.0; // turbulent viscosity
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    primvar_j[iVar] = primvar_i[iVar];
  }

  const su2double tke = 3; // 3 cancels out 3 in denominator

  /*---
   * TEST
   * ---*/

  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

  su2double expected_residual[nVar] = {0, -6, 0, 0, 0};
  expected_jacobian_i[1][1] = -8;
  expected_jacobian_i[2][2] = -6;
  expected_jacobian_i[3][3] = -6;
  expected_jacobian_i[4][0] = 26.25;
  expected_jacobian_i[4][4] = -10.5;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      expected_jacobian_j[iVar][jVar] = -expected_jacobian_i[iVar][jVar];
    }
  }
  expected_jacobian_i[4][1] = -3;
  expected_jacobian_j[4][1] = -3;

//  cout << "Calculated:\n";
//  PrintInformation(residual_i, Jacobian_i, Jacobian_j);

  const su2double tolerance = 100*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    BOOST_CHECK_SMALL(expected_residual[iVar] - residual_i[iVar], tolerance);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      // cout << "iVar: " << iVar << "\tjVar: " << jVar << "\n";
      BOOST_CHECK_SMALL(expected_jacobian_i[iVar][jVar] - Jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_SMALL(expected_jacobian_j[iVar][jVar] - Jacobian_j[iVar][jVar], tolerance);
    }
  }
}

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
  primvar_grad_i[3][2] =  3.0; // dw/dz
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
  numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);

  su2double expected_residual[nVar] = {0, -18, 12, 0, 6};
  expected_jacobian_i[1][0] = 8;
  expected_jacobian_i[1][1] = -8;
  expected_jacobian_i[2][0] = 12;
  expected_jacobian_i[2][2] = -6;
  expected_jacobian_i[3][0] = 18;
  expected_jacobian_i[3][3] = -6;
  expected_jacobian_i[4][4] = -10.5;
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      expected_jacobian_j[iVar][jVar] = -expected_jacobian_i[iVar][jVar];
    }
  }
  expected_jacobian_i[4][0] = 41.75;
  expected_jacobian_i[4][1] = -6.5;
  expected_jacobian_i[4][2] = 15;
  expected_jacobian_i[4][3] = 13.5;

  expected_jacobian_j[4][0] = -35.75;
  expected_jacobian_j[4][1] = -11.5;
  expected_jacobian_j[4][2] = -3;
  expected_jacobian_j[4][3] = -13.5;

  cout << "Calculated:\n";
  PrintInformation(residual_i, Jacobian_i, Jacobian_j);

  const su2double tolerance = 100*std::numeric_limits<su2double>::epsilon();
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    BOOST_CHECK_SMALL(expected_residual[iVar] - residual_i[iVar], tolerance);
    for (unsigned short jVar = 0; jVar < nVar; jVar++) {
      cout << "iVar: " << iVar << "\tjVar: " << jVar << "\n";
      BOOST_CHECK_SMALL(expected_jacobian_i[iVar][jVar] - Jacobian_i[iVar][jVar], tolerance);
      BOOST_CHECK_SMALL(expected_jacobian_j[iVar][jVar] - Jacobian_j[iVar][jVar], tolerance);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(ViscousTiming, ViscousResidualFixture) {

  /*---
   * SETUP
   * ---*/

  CNumerics* numerics = new CAvgGrad_Flow(3, 5, config);

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
  primvar_grad_i[3][2] =  3.0; // dw/dz
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

  numerics->SetCoord(coord_i, coord_j);
  numerics->SetNormal(normal);
  numerics->SetSecondary(NULL, NULL);
  numerics->SetPrimitive(primvar_i, primvar_j);
  numerics->SetPrimVarGradient(primvar_grad_i, primvar_grad_j);
  numerics->SetTurbKineticEnergy(tke, tke);
  clock_t begin = clock();
  for (unsigned long i = 0; i < 1E6; i++)
    numerics->ComputeResidual(residual_i, Jacobian_i, Jacobian_j, config);
  clock_t end = clock();
  double elapsed_secs = double(end - begin) / CLOCKS_PER_SEC;
  std::cout << "Elapsed time: " << elapsed_secs << std::endl;

  delete numerics;
}

#endif

/*!
 * \file cumulative_average_test.cpp
 * \brief Tests the cumulative moving average
 * \author C. Pederson
 * \version 6.1.0 "Falcon"
 *
 * The current SU2 release has been coordinated by the
 * SU2 International Developers Society <www.su2devsociety.org>
 * with selected contributions from the open-source community.
 *
 * The main research teams contributing to the current release are:
 *  - Prof. Juan J. Alonso's group at Stanford University.
 *  - Prof. Piero Colonna's group at Delft University of Technology.
 *  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
 *  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
 *  - Prof. Rafael Palacios' group at Imperial College London.
 *  - Prof. Vincent Terrapon's group at the University of Liege.
 *  - Prof. Edwin van der Weide's group at the University of Twente.
 *  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
 *
 * Copyright 2012-2018, Francisco D. Palacios, Thomas D. Economon,
 *                      Tim Albring, and the SU2 contributors.
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

#include <limits> // used to find machine epsilon

#include "../include/solver_structure.hpp"
#include "../include/solver_structure_v2f.hpp"
#include "../include/variable_structure_v2f.hpp"

namespace cumulative_average_test {

const unsigned short nDim = 3;
const unsigned short nVar = nDim+2;
const unsigned short nPrimVar = nDim+9;
const unsigned short nSecVar = 4;

/**
 * Write a cfg file to be used in initializing the CConfig object.
 */
static void WriteCfgFile(const char* filename, const string& period) {

  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "TIME_DISCRE_FLOW= EULER_IMPLICIT" << std::endl;

  cfg_file.close();

}

class CTestGeometry : public CGeometry {
 public:
  CTestGeometry(CConfig* config) : CGeometry() {
    nPoint = 1;
    node = new CPoint*[nPoint];
    for (unsigned short iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CPoint(0, 0, 0, iPoint, config);
      node[iPoint]->SetVolume(1);
    }
  };
};

class CTestSolver : public CNSSolver {
 public:
  CTestSolver(su2double* initial_solution,
              unsigned short val_nDim, unsigned short val_nVar,
              CConfig* config)
       : CNSSolver() {

    nDim = val_nDim;
    nVar = val_nVar;
    nPoint = 1;
    nPointDomain = 1;

    Solution = new su2double[val_nVar];

    node = new CVariable*[nPoint];
    average_node = NULL;

    su2double* velocity = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      velocity[iDim] = initial_solution[iDim+1]/initial_solution[0];
    }

    const su2double density = initial_solution[0];
    const su2double energy = initial_solution[nDim+1];

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CNSVariable(density, velocity, energy, nDim, nVar, config);
    }

    delete [] velocity;
  };

  ~CTestSolver() {
  }

};


struct CMA_Fixture {
  CMA_Fixture() {

    char cfg_filename[100] = "averaging_timescale_test.cfg";
    WriteCfgFile(cfg_filename, "FLOW_TIMESCALE");
    const unsigned short iZone = 0;
    const unsigned short nZone = 1;
    config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
    std::remove(cfg_filename);

    su2double solution[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
    solver = new CTestSolver(solution, nDim, nVar, config);
    solver->UpdateCMAverage();

  }
  ~CMA_Fixture() {
    delete config;
    delete solver;
  }

  CConfig* config;
  CTestSolver* solver;
};


BOOST_AUTO_TEST_SUITE(CumulativeMovingAverage);

BOOST_FIXTURE_TEST_CASE(CumulativeAverageOfDensity, CMA_Fixture) {

  /*--- Act ---*/

  su2double solution[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
  solution[0] = 3.0;
  solver->node[0]->SetSolution(solution);
  solver->UpdateCMAverage();
  const su2double density_1 = solver->node[0]->GetCMA_variables()[0];

  solution[0] = 5.0;
  solver->node[0]->SetSolution(solution);
  solver->UpdateCMAverage();
  const su2double density_2 = solver->node[0]->GetCMA_variables()[0];

  /*--- Assert ---*/

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(density_1, 2.0, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(density_2, 3.0, tolerance);

}

BOOST_FIXTURE_TEST_CASE(CumulativeAverageOfMomentum, CMA_Fixture) {

  /*--- Act ---*/

  su2double solution[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
  solution[1] = 2.0;
  solver->node[0]->SetSolution(solution);
  solver->UpdateCMAverage();
  const su2double momentum_1 = solver->node[0]->GetCMA_variables()[1];

  solution[1] = 4.0;
  solver->node[0]->SetSolution(solution);
  solver->UpdateCMAverage();
  const su2double momentum_2 = solver->node[0]->GetCMA_variables()[1];

  /*--- Assert ---*/
  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(momentum_1, 1.0, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(momentum_2, 2.0, tolerance);

}

BOOST_FIXTURE_TEST_CASE(CumulativeAverageOfUU, CMA_Fixture) {

  /*--- Act ---*/

  su2double solution[5] = {1.0, 0.0, 0.0, 0.0, 0.0};
  solution[1] = 3.0;
  solver->node[0]->SetSolution(solution);
  solver->node[0]->SetVelocity();
  solver->UpdateCMAverage();
  const su2double uu_1 = solver->node[0]->GetCMA_variables()[5];

  solution[1] = 3.0;
  solver->node[0]->SetSolution(solution);
  solver->node[0]->SetVelocity();
  solver->UpdateCMAverage();
  const su2double uu_2 = solver->node[0]->GetCMA_variables()[5];

  /*--- Assert ---*/

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(uu_1, 4.5, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(uu_2, 6.0, tolerance);

}


BOOST_AUTO_TEST_SUITE_END();

} // end namespace

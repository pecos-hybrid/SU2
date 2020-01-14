/*!
 * \file averaging_timescale_test.cpp
 * \brief Tests the averaging timescale used
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

#define BOOST_TEST_MODULE ViscousIdealVsGeneral
#include "MPI_global_fixture.hpp"


#include <limits> // used to find machine epsilon

#include "../include/solver_structure.hpp"
#include "../include/solver_structure_v2f.hpp"
#include "../include/variable_structure_v2f.hpp"

const unsigned short nDim = 3;
const unsigned short nVar = nDim+2;
const unsigned short nPrimVar = nDim+9;
const unsigned short nSecVar = 4;

/**
 * Write a cfg file to be used in initializing the CConfig object.
 */
void WriteCfgFile(const char* filename, const string& period) {

  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "TIME_DISCRE_FLOW= EULER_IMPLICIT" << std::endl;
  cfg_file << "NUM_AVERAGING_PERIODS= 4.0" << std::endl;
  cfg_file << "AVERAGING_PERIOD= " << period << std::endl;

  cfg_file.close();

}

class CTestGeometry : public CGeometry {
 public:
  CTestGeometry(CConfig* config) : CGeometry() {
    nPoint = 3;
    node = new CPoint*[nPoint];
    for (unsigned short iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CPoint(0, 0, 0, iPoint, config);
      node[iPoint]->SetVolume(1);
    }
  };
};

class CTestSolver : public CNSSolver {
 public:
  CTestSolver(su2double density, su2double average_density,
              su2double laminar_viscosity,
              unsigned short val_nDim, unsigned short val_nVar,
              CConfig* config)
       : CNSSolver() {

    nDim = val_nDim;
    nVar = val_nVar;
    nPoint = 1;
    nPointDomain = 1;

    Solution = new su2double[val_nVar];

    node = new CVariable*[nPoint];
    average_node = new CVariable*[nPoint];

    su2double* velocity = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      velocity[iDim] = 0;
    }

    node[0] = new CNSVariable(density, velocity, 0, nDim, nVar, config);
    average_node[0] =
        new CNSVariable(average_density, velocity, 0, nDim, nVar, config);

    node[0]->SetLaminarViscosity(laminar_viscosity);
    average_node[0]->SetLaminarViscosity(laminar_viscosity);

    delete [] velocity;
  };
};

class CTestTurbSolver : public CTurbKESolver {
 public:
  CTestTurbSolver(su2double kine, su2double epsi,
                  unsigned short val_nDim, CConfig* config)
       : CTurbKESolver() {

    nDim = val_nDim;
    nVar = 4;
    nPoint = 1;
    nPointDomain = 1;

    Solution = new su2double[nVar];

    node = new CVariable*[nPoint];

    const su2double zeta = 1;
    const su2double f = 1;
    const su2double muT = 1;
    const su2double Tm = 1;
    const su2double Lm = 1;

    node[0] = new CTurbKEVariable(kine, epsi, zeta, f, muT, Tm, Lm,
                                  nDim, nVar, config);
  };
};

// Setup MPI
BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

/**
 */
BOOST_AUTO_TEST_CASE(FlowTimescale) {

  /*--- Setup ---*/

  char cfg_filename[100] = "averaging_timescale_test.cfg";
  WriteCfgFile(cfg_filename, "FLOW_TIMESCALE");
  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
  std::remove(cfg_filename);

  CGeometry* geometry = NULL;

  CSolver** solver_container = new CSolver*[1];
  const su2double inst_density = 1.0;
  const su2double avg_density = 0.0;
  const su2double mu = 0.0;
  solver_container[0] =
        new CTestSolver(inst_density, avg_density, mu, nDim, nVar, config);

  // Averaging relies on having an initialized, physical time
  config->SetCurrent_UnstTime(1.0);
  config->SetDelta_UnstTimeND(0.04);
  config->SetLength_Ref(1.0);
  config->SetVelocity_FreeStreamND(1/0.99, 0);
  config->SetVelocity_FreeStreamND(0.0,    1);
  config->SetVelocity_FreeStreamND(0.0,    2);
  config->SetModVel_FreeStream(1/0.99);
  config->SetTime_Ref(1.0);

  /*--- Test ---*/

  /*--- Weight should be 0.04 / (4*(0.99 + 0.01)) = 0.01 --*/

  solver_container[0]->SetAverages(geometry, solver_container, config);

  const su2double* average =
      solver_container[0]->average_node[0]->GetSolution();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(average[0], 0.01, tolerance);

  /*--- Teardown ---*/

  delete solver_container[0];
  delete [] solver_container;

}

BOOST_AUTO_TEST_CASE(TypicalTurbTimescale) {

  /*--- Setup ---*/

  char cfg_filename[100] = "averaging_timescale_test.cfg";
  WriteCfgFile(cfg_filename, "TURB_TIMESCALE");
  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
  std::remove(cfg_filename);

  /*--- Initialize some cfg values to needed for constructors ---*/

  const su2double vel_freestream = 1/0.99;
  config->SetModVel_FreeStreamND(vel_freestream);

  CGeometry* geometry = NULL;

  CSolver** solver_container = new CSolver*[TURB_SOL + 1];
  const su2double inst_density = 2.0;
  const su2double avg_density = 1.0;
  const su2double kine = 0.99;
  const su2double eps = 1.0;
  const su2double mu = 0.0;
  solver_container[FLOW_SOL] =
        new CTestSolver(inst_density, avg_density, mu, nDim, nVar, config);
  CTestTurbSolver* turb_solver = new CTestTurbSolver(kine, eps, nDim, config);
  turb_solver->CalculateTurbScales(solver_container, config);
  solver_container[TURB_SOL] = turb_solver;

  // Averaging relies on having an initialized, physical time
  config->SetCurrent_UnstTime(1.0);
  config->SetDelta_UnstTimeND(0.04);
  config->SetLength_Ref(1.0);
  config->SetVelocity_FreeStreamND(vel_freestream, 0);
  config->SetVelocity_FreeStreamND(0.0,    1);
  config->SetVelocity_FreeStreamND(0.0,    2);
  config->SetModVel_FreeStreamND(vel_freestream);
  config->SetTime_Ref(1.0);

  /*--- Test ---*/

  /*---
   * Turb timescale used:
   * k/eps = 0.99
   * kolmogorov timescale = 0.0
   *
   * Weight should be 1 + 0.04 / (4*(0.99 + 0.01)) = 1.01 --*/

  solver_container[0]->SetAverages(geometry, solver_container, config);

  const su2double* average =
    solver_container[0]->average_node[0]->GetSolution();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(average[0], 1.01, tolerance);

  /*--- Teardown ---*/

  delete solver_container[0];
  delete [] solver_container;

}

BOOST_AUTO_TEST_CASE(KolTurbTimescale) {

  /*--- Setup ---*/

  char cfg_filename[100] = "averaging_timescale_test.cfg";
  WriteCfgFile(cfg_filename, "TURB_TIMESCALE");
  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
  std::remove(cfg_filename);

  /*--- Initialize some cfg values to needed for constructors ---*/

  const su2double vel_freestream = 1/0.99;
  config->SetModVel_FreeStreamND(vel_freestream);

  CGeometry* geometry = NULL;

  CSolver** solver_container = new CSolver*[TURB_SOL + 1];
  const su2double inst_density = 2.0;
  const su2double avg_density = 1.0;
  const su2double kine = 0.0;
  const su2double eps = 1.0;

  /*--- Craft mu such that T = 0.99. Note that since no model-split is
   * specified in the cfg file, this will use instantaneous density.  ---*/

  const su2double desired_timescale = 0.99;
  const su2double C_T = 6.0;
  const su2double mu = inst_density * eps * pow(desired_timescale / C_T, 2);

  solver_container[FLOW_SOL] =
        new CTestSolver(inst_density, avg_density, mu, nDim, nVar, config);
  CTestTurbSolver* turb_solver = new CTestTurbSolver(kine, eps, nDim, config);
  turb_solver->CalculateTurbScales(solver_container, config);
  solver_container[TURB_SOL] = turb_solver;

  // Averaging relies on having an initialized, physical time
  config->SetCurrent_UnstTime(1.0);
  config->SetDelta_UnstTimeND(0.04);
  config->SetLength_Ref(1.0);
  config->SetVelocity_FreeStreamND(vel_freestream, 0);
  config->SetVelocity_FreeStreamND(0.0,    1);
  config->SetVelocity_FreeStreamND(0.0,    2);
  config->SetModVel_FreeStreamND(vel_freestream);
  config->SetTime_Ref(1.0);

  /*--- Test ---*/

  /*---
   * Turb timescale used:
   * k/eps = 0
   * kolmogorov timescale = 0.99
   *
   * Weight should be 1 + 0.04 / (4*(0.99 + 0.01)) = 1.01 --*/

  solver_container[0]->SetAverages(geometry, solver_container, config);

  const su2double* average =
    solver_container[0]->average_node[0]->GetSolution();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(average[0], 1.01, tolerance);

  /*--- Teardown ---*/

  delete solver_container[0];
  delete [] solver_container;

}

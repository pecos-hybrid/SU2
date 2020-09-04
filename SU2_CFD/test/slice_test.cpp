/*!
 * \file slice_test.cpp
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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include <limits> // used to find machine epsilon

#include "../include/solver_structure.hpp"
#include "../include/solver_structure_v2f.hpp"
#include "../include/variable_structure_v2f.hpp"
#include "../include/slice_file_reader.hpp"

namespace slice_test {

const unsigned short nDim = 3;

/**
 * Write a cfg file to be used in initializing the CConfig object.
 */
static void WriteCfgFile(const char* filename) {

  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file.close();

}

static void WriteSliceFile(const char* filename) {
  std::ofstream slice_file;
  slice_file.open(filename, ios::out);
  slice_file << "\"x\"	\"y\"	\"Density\"	\"Momentum_x\"	\"Momentum_y\"	\"Momentum_z\"	\"Energy\"	\"TKE\"	\"Dissipation\"	\"v2\"	\"f\"" << std::endl;
  slice_file << "-3.00000000000000000e+00	1.00000000000000000e+00	1.00000000000000000e+00	0.00000000000000000e+00	3.00000000000000000e+00	0.00000000000000000e+00	2.83000000000000000e+00	0.00000000000000000e+00	3.35000000000000000e-05	1.00000000000000000e+00	0.00000000000000000e+00" << std::endl;
  slice_file << "-2.00000000000000000e+00	2.00000000000000000e+00	1.00000000000000000e+00	1.00000000000000000e+00	4.00000000000000000e+00	0.00000000000000000e+00	2.83000000000000000e+00	0.00000000000000000e+00	3.35000000000000000e-05	2.00000000000000000e+00	0.00000000000000000e+00" << std::endl;
  slice_file << "-1.00000000000000000e+00	3.00000000000000000e+00	1.00000000000000000e+00	2.00000000000000000e+00	5.00000000000000000e+00	0.00000000000000000E+00	2.83000000000000000e+00	0.00000000000000000e+00	3.35000000000000000e-05	3.00000000000000000e+00	0.00000000000000000e+02" << std::endl;
  slice_file.close();
}

class CTestGeometry : public CGeometry {
 public:
  CTestGeometry(CConfig* config) : CGeometry() {
    nPoint = 9;
    nPointDomain = nPoint;
    Global_nPointDomain = nPoint;
    nDim = 3;
    node = new CPoint*[nPoint];
    for (unsigned short i = 0; i < 3; i++) {
      node[i*3+0] = new CPoint(-3, 1, i, 0, config);
      node[i*3+1] = new CPoint(-2, 2, i, 1, config);
      node[i*3+2] = new CPoint(-1, 3, i, 2, config);
    }
    for (unsigned short iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint]->SetVolume(1);
    }
  };

  long GetGlobal_to_Local_Point(unsigned long iPoint_Global) override {
    return iPoint_Global;
  };
  unsigned long GetGlobal_nPointDomain(void) override { return Global_nPointDomain; }
};

class CTestSolver : public CNSSolver {
 public:
  CTestSolver(CConfig* config)
       : CNSSolver() {

    nDim = 3;
    nVar = 5;
    nPoint = 9;
    nPointDomain = 9;

    Solution = new su2double[nVar];

    node = new CVariable*[nPoint];
    average_node = NULL;

    su2double* velocity = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      velocity[iDim] = 0;
    }

    for (unsigned short iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CNSVariable(0, velocity, 0, nDim, nVar, config);
    }

    delete [] velocity;

  };

};

class CTestTurbSolver : public CTurbKESolver {
 public:
  CTestTurbSolver(CConfig* config)
       : CTurbKESolver() {

    nDim = 3;
    nVar = 4;
    nPoint = 9;
    nPointDomain = 9;

    Solution = new su2double[nVar];

    node = new CVariable*[nPoint];

    const su2double kine = 0;
    const su2double epsi = 1;
    const su2double zeta = 0;
    const su2double f = 0;
    const su2double muT = 0;
    const su2double Tm = 1;
    const su2double Lm = 1;

    for (unsigned short iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CTurbKEVariable(kine, epsi, zeta, f, muT, Tm, Lm,
                                    nDim, nVar, config);
    }
  };
};

struct SliceFixture {
  SliceFixture() {

    char cfg_filename[100] = "slice_test.cfg";
    WriteCfgFile(cfg_filename);
    const unsigned short iZone = 0;
    const unsigned short nZone = 1;
    config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
    std::remove(cfg_filename);

    geometry = new CTestGeometry(config);

    solver_container = new CSolver*[TURB_SOL+1];
    solver_container[FLOW_SOL] = new CTestSolver(config);
    flow_sol = dynamic_cast<CTestSolver*>(solver_container[FLOW_SOL]);
    solver_container[TURB_SOL] = new CTestTurbSolver(config);
    turb_sol = dynamic_cast<CTestTurbSolver*>(solver_container[TURB_SOL]);

    slice_filename = "slice_test_data.dat";
    WriteSliceFile(slice_filename.c_str());
  }

  ~SliceFixture() {
    delete solver_container[FLOW_SOL];
    delete solver_container[TURB_SOL];
    delete [] solver_container;

    delete geometry;

    delete config;
  }

  std::string slice_filename;
  CConfig* config;
  CGeometry* geometry;
  CSolver** solver_container;
  CTestSolver* flow_sol;
  CTestTurbSolver* turb_sol;
};

BOOST_AUTO_TEST_SUITE(ReadSliceFile);

/**
 */
BOOST_FIXTURE_TEST_CASE (SimpleRead, SliceFixture) {
  CFileReader_Cartesian file_reader;

  file_reader.Read_SliceFile_ASCII(config, slice_filename);

  const passivedouble* restart_data = file_reader.GetRestart_Data();
  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(restart_data[0], -3.0, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(restart_data[0+11], -2.0, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(restart_data[0+22], -1.0, tolerance);

}

/**
 */
BOOST_FIXTURE_TEST_CASE(ReadFlowVarsInCartesian, SliceFixture) {
  CFileReader_Cartesian file_reader;

  file_reader.Read_SliceFile_ASCII(config, slice_filename);
  file_reader.LoadSolutionFromSlice(slice_filename, config, geometry, FLOW_SOL, flow_sol->node);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  {
    const su2double* solution = flow_sol->node[0]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[1], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 3.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[1]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 4.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[2]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 2.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 5.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }

  /*--- Points 3, 4, 5 have the same x, y as points 0, 1, 2, but
   * they have a different z.  So the solution should be identical ---*/

  {
    const su2double* solution = flow_sol->node[3]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[1], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 3.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[4]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 4.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[5]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 2.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 5.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }

}

BOOST_FIXTURE_TEST_CASE(ReadFlowVarsInCylindrical, SliceFixture) {
  CFileReader_Cylindrical file_reader;

  /*--- Instead of y = 0, 1, 2, 0, 1, 2, 0, 1, 2
   *               z = 0, 0, 0, 1, 1, 1, 2, 2, 2
   *    Set        r = 0, 1, 2, 0, 1, 2, 0, 1, 2
   *               theta = 0, 0, 0, pi/2, pi/2, pi/2, pi, pi, pi
   *---*/
  for (unsigned short i = 0; i < 3; i++) {
    for (unsigned short j = 0; j < 3; j++) {
      const su2double r = 1 + j;
      const su2double theta = M_PI*i/2.0;
      const su2double y = r*cos(theta);
      const su2double z = r*sin(theta);
      su2double* coord = geometry->node[i*3+j]->GetCoord();
      coord[0] = -3+j; coord[1] = y; coord[2] = z;
    }
  }

  file_reader.Read_SliceFile_ASCII(config, slice_filename);
  file_reader.LoadSolutionFromSlice(slice_filename, config, geometry, FLOW_SOL, flow_sol->node);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  {
    const su2double* solution = flow_sol->node[0]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[1], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[1]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 1.0, tolerance);
  }
  {
    const su2double* solution = flow_sol->node[2]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 2.0, tolerance);
  }

  /*--- points 3, 4, 5 have the same x, r as points 0, 1, 2, but
   * they have a different theta.  so the solution should be identical ---*/

  {
    const su2double* solution = flow_sol->node[3]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[1], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[4]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 1.0, tolerance);
  }
  {
    const su2double* solution = flow_sol->node[5]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[0], 1.0, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 2.0, tolerance);
  }
}

BOOST_FIXTURE_TEST_CASE(TransformFlowVarsInCylindrical, SliceFixture) {
  CFileReader_Cylindrical file_reader;

  /*--- Instead of y = 0, 1, 2, 0, 1, 2, 0, 1, 2
   *               z = 0, 0, 0, 1, 1, 1, 2, 2, 2
   *    Set        r = 0, 1, 2, 0, 1, 2, 0, 1, 2
   *               theta = 0, 0, 0, pi/2, pi/2, pi/2, pi, pi, pi
   *---*/
  for (unsigned short i = 0; i < 3; i++) {
    for (unsigned short j = 0; j < 3; j++) {
      const su2double r = 1 + j;
      const su2double theta = M_PI*i/2.0;
      const su2double y = r*cos(theta);
      const su2double z = r*sin(theta);
      su2double* coord = geometry->node[i*3+j]->GetCoord();
      coord[0] = -3+j; coord[1] = y; coord[2] = z;
    }
  }

  file_reader.Read_SliceFile_ASCII(config, slice_filename);
  file_reader.LoadSolutionFromSlice(slice_filename, config, geometry, FLOW_SOL, flow_sol->node);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  {
    const su2double* solution = flow_sol->node[0]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 3.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[1]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 4.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[2]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 5.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }

  /*--- Points 3, 4, 5 have the same x, r as points 0, 1, 2, but
   * they have a different theta.  So the momentum should be rotated ---*/

  {
    const su2double* solution = flow_sol->node[3]->GetSolution();
    BOOST_CHECK_SMALL(solution[2], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[3], 3.0, tolerance);
  }
  {
    const su2double* solution = flow_sol->node[4]->GetSolution();
    BOOST_CHECK_SMALL(solution[2], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[3], 4.0, tolerance);
  }
  {
    const su2double* solution = flow_sol->node[5]->GetSolution();
    BOOST_CHECK_SMALL(solution[2], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[3], 5.0, tolerance);
  }

  /*--- Points 3, 4, 5 have the same x, r as points 0, 1, 2, but
   * they have a different theta.  So the momentum should be rotated ---*/

  {
    const su2double* solution = flow_sol->node[6]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[2], -3.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[7]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[2], -4.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = flow_sol->node[8]->GetSolution();
    BOOST_CHECK_CLOSE_FRACTION(solution[2], -5.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
}

BOOST_FIXTURE_TEST_CASE(ReadTurbVarsInCartesian, SliceFixture) {
  CFileReader_Cartesian file_reader;

  file_reader.Read_SliceFile_ASCII(config, slice_filename);
  file_reader.LoadSolutionFromSlice(slice_filename, config, geometry, TURB_SOL, turb_sol->node);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  {
    const su2double* solution = turb_sol->node[0]->GetSolution();
    BOOST_CHECK_SMALL(solution[0], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 3.35e-05, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
}

BOOST_FIXTURE_TEST_CASE(ReadTurbVarsInPolar, SliceFixture) {
  CFileReader_Cylindrical file_reader;

  /*--- Instead of y = 0, 1, 2, 0, 1, 2, 0, 1, 2
   *               z = 0, 0, 0, 1, 1, 1, 2, 2, 2
   *    Set        r = 0, 1, 2, 0, 1, 2, 0, 1, 2
   *               theta = 0, 0, 0, pi/2, pi/2, pi/2, pi, pi, pi
   *---*/
  for (unsigned short i = 0; i < 3; i++) {
    for (unsigned short j = 0; j < 3; j++) {
      const su2double r = 1 + j;
      const su2double theta = M_PI*i/2.0;
      const su2double y = r*cos(theta);
      const su2double z = r*sin(theta);
      su2double* coord = geometry->node[i*3+j]->GetCoord();
      coord[0] = -3+j; coord[1] = y; coord[2] = z;
    }
  }

  file_reader.Read_SliceFile_ASCII(config, slice_filename);
  file_reader.LoadSolutionFromSlice(slice_filename, config, geometry, TURB_SOL, turb_sol->node);

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  {
    const su2double* solution = turb_sol->node[0]->GetSolution();
    BOOST_CHECK_SMALL(solution[0], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 3.35e-05, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = turb_sol->node[1]->GetSolution();
    BOOST_CHECK_SMALL(solution[0], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 3.35e-05, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 2.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  /*--- points 3, 4, 5 have the same x, r as points 0, 1, 2, but
   * they have a different theta.  so the solution should be identical ---*/
  {
    const su2double* solution = turb_sol->node[3]->GetSolution();
    BOOST_CHECK_SMALL(solution[0], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 3.35e-05, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 1.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
  {
    const su2double* solution = turb_sol->node[4]->GetSolution();
    BOOST_CHECK_SMALL(solution[0], tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[1], 3.35e-05, tolerance);
    BOOST_CHECK_CLOSE_FRACTION(solution[2], 2.0, tolerance);
    BOOST_CHECK_SMALL(solution[3], tolerance);
  }
}

BOOST_AUTO_TEST_SUITE_END();

} // end namespace

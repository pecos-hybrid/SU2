/*!
 * \file fluctuating_stress_test.cpp
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

#define BOOST_TEST_MODULE ConvectiveBlending
#include "boost/test/included/unit_test.hpp"

#include <cstdio> // std::remove
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/numerics_structure.hpp"

#include "../../Common/test/MPI_global_fixture.hpp"

#include <cstdio> // std::remove
#include <fstream>

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

void WriteCfgFile(unsigned short nDim, const char* filename,
                  std::string blending) {
  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "ROE_LOW_DISSIPATION= " << blending << std::endl;

  cfg_file.close();
}

BOOST_AUTO_TEST_CASE(BadSensorsAllowedForNTS) {

  /*--- Setup ---*/

  const unsigned short nDim = 3;

  char cfg_filename[100] = "convective_blending_test.cfg";
  WriteCfgFile(nDim, cfg_filename, "NTS");
  CConfig* config = new CConfig(cfg_filename, SU2_CFD, 0, 1, 2, VERB_NONE);
  std::remove(cfg_filename);

  const su2double dissipation_i = 0.4;
  const su2double dissipation_j = 0.6;
  su2double dissipation;
  // Intentionally unphysical:
  const su2double sensor_i = NAN;
  const su2double sensor_j = NAN;

  /*--- Test ---*/

  CNumerics numerics;
  numerics.SetRoe_Dissipation(dissipation_i, dissipation_j,
                              sensor_i, sensor_j,
                              dissipation, config);

  const su2double tolerance = std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(dissipation, 0.5, tolerance);

  /*--- Teardown ---*/
  delete config;
}

BOOST_AUTO_TEST_CASE(NTSAndFDGiveSameAverageDissipation) {

  /*--- Setup ---*/

  const unsigned short nDim = 3;

  char cfg_filename[100] = "convective_blending_test.cfg";
  WriteCfgFile(nDim, cfg_filename, "NTS");
  CConfig* nts_config = new CConfig(cfg_filename, SU2_CFD, 0, 1, 2, VERB_NONE);
  std::remove(cfg_filename);

  WriteCfgFile(nDim, cfg_filename, "FD");
  CConfig* fd_config = new CConfig(cfg_filename, SU2_CFD, 0, 1, 2, VERB_NONE);
  std::remove(cfg_filename);

  const su2double dissipation_i = 0.2;
  const su2double dissipation_j = 0.4;
  const su2double sensor_i = 0;
  const su2double sensor_j = 0;
  su2double nts_dissipation, fd_dissipation;

  /*--- Test ---*/

  CNumerics numerics;
  numerics.SetRoe_Dissipation(dissipation_i, dissipation_j,
                              sensor_i, sensor_j,
                              nts_dissipation, nts_config);
  numerics.SetRoe_Dissipation(dissipation_i, dissipation_j,
                              sensor_i, sensor_j,
                              fd_dissipation, fd_config);

  const su2double tolerance = std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(fd_dissipation, nts_dissipation, tolerance);

  /*--- Teardown ---*/
  delete fd_config;
  delete nts_config;
}

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

#define BOOST_TEST_MODULE hybrid_rdelta
#include "MPI_global_fixture.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../include/solver_structure.hpp"
#include "../include/variable_structure.hpp"
#include "../include/hybrid_RANS_LES_model.hpp"

void WriteCfgFile(const unsigned short& nDim) {
  std::ofstream cfg_file;

  cfg_file.open("hybrid_rdelta_test.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
  cfg_file << "HYBRID_TURB_MODEL= YES" << std::endl;
  cfg_file << "HYBRID_BLENDING_SCHEME= CONVECTIVE" << std::endl;
  cfg_file << "HYBRID_RESOLUTION_INDICATOR= RDELTA" << std::endl;
  cfg_file << "HYBRID_ANISOTROPY_MODEL= ISOTROPIC" << std::endl;
  cfg_file.close();

}


/* ----------------------------------------------------------------------------
 *  Test Fixtures
 * --------------------------------------------------------------------------*/

struct HybridRdeltaFixture {
  HybridRdeltaFixture()
    : machine_eps(std::numeric_limits<su2double>::epsilon())
  {

    WriteCfgFile(3);
    mock_config = new CConfig("hybrid_rdelta_test.cfg", SU2_CFD, 0, 1, 3, VERB_NONE);

    su2double flow_soln[5], hybrid_soln;

    flow_soln[0] = 1.0;
    flow_soln[1] = 0.0;
    flow_soln[2] = 0.0;
    flow_soln[3] = 0.0;
    flow_soln[4] = 2.5;

    hybrid_soln = 1.0;

    // Alloc mock solver objects
    mock_solver_array = new CSolver* [MAX_SOLS];
    for (unsigned short ii=0; ii<MAX_SOLS; ii++) {
      mock_solver_array[ii] = NULL;
    }
    mock_solver_array[FLOW_SOL]   = new CNSSolver();
    mock_solver_array[TURB_SOL]   = new CTurbKESolver();
    mock_solver_array[HYBRID_SOL] = new CHybridConvSolver();

    // Alloc mock node within the solver
    mock_solver_array[FLOW_SOL]->node    = new CVariable*[1];
    mock_solver_array[FLOW_SOL]->node[0] =
      new CNSVariable(flow_soln, 3, 5, mock_config);

    mock_solver_array[TURB_SOL]->node    = new CVariable*[1];
    mock_solver_array[TURB_SOL]->node[0] =
      new CTurbKEVariable(0.0, 0.0, 0.0, 0.0, // k, eps, v2, f
                          0.0, 0.0, 0.0,      // muT, Tm, Lm
                          3, 4,               // ndim, nvar
                          NULL, mock_config); // constants, config

    mock_solver_array[HYBRID_SOL]->node    = new CVariable*[1];
    mock_solver_array[HYBRID_SOL]->node[0] =
      new CHybridConvVariable(hybrid_soln, 3, 1, mock_config);

    mock_mediator = new CHybrid_Mediator(3, mock_config);
  }

  ~HybridRdeltaFixture() {
    delete mock_mediator;

    delete mock_solver_array[HYBRID_SOL]->node[0];
    delete mock_solver_array[HYBRID_SOL]->node;
    mock_solver_array[HYBRID_SOL]->node = NULL; // avoid double free

    delete mock_solver_array[TURB_SOL]->node[0];
    delete mock_solver_array[TURB_SOL]->node;
    mock_solver_array[TURB_SOL]->node = NULL; // avoid double free

    delete mock_solver_array[FLOW_SOL]->node[0];
    delete mock_solver_array[FLOW_SOL]->node;
    mock_solver_array[FLOW_SOL]->node = NULL; // avoid double free

    delete mock_solver_array[HYBRID_SOL];
    delete mock_solver_array[TURB_SOL];
    delete mock_solver_array[FLOW_SOL];

    delete mock_solver_array;
    delete mock_config;
  }

  const su2double machine_eps;

  CConfig *mock_config;
  CSolver **mock_solver_array;
  CHybrid_Mediator *mock_mediator;
};

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

#ifdef BUILD_TESTS

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

BOOST_FIXTURE_TEST_CASE(ZeroGradientTrivial, HybridRdeltaFixture) {

  BOOST_TEST_MESSAGE( "Running ZeroGradientTrival test..." );

  mock_solver_array[FLOW_SOL]->node[0]->SetGradient_PrimitiveZero(7);
  mock_solver_array[FLOW_SOL]->node[0]->SetEddyViscosity(0.0);
  mock_mediator->ComputeProdLengthTensor(NULL, mock_solver_array, 0);

  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,2),0.0);

}

#endif

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

#include "../include/solver_structure_v2f.hpp"
#include "../include/variable_structure_v2f.hpp"
#include "../include/hybrid_RANS_LES_model.hpp"

void WriteCfgFile(const unsigned short& nDim) {
  std::ofstream cfg_file;

  cfg_file.open("hybrid_rdelta_test.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
  cfg_file << "HYBRID_RANSLES= DYNAMIC_HYBRID" << std::endl;
  cfg_file << "HYBRID_RESOLUTION_INDICATOR= RDELTA_STRAIN_ONLY" << std::endl;
  cfg_file << "HYBRID_ANISOTROPY_MODEL= ISOTROPIC" << std::endl;
  // This option is deprecated
  // cfg_file << "HYBRID_MODEL_CONSTANT= 1.0" << std::endl;
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
    mock_var_array = new CVariable* [3];
    for (unsigned short ii=0; ii<3; ii++) {
      mock_var_array[ii] = NULL;
    }

    // Alloc mock node within the solver
    mock_var_array[0] = new CNSVariable(flow_soln, 3, 5, mock_config);

    mock_var_array[1] = new CTurbKEVariable(0.0, 0.0, 0.0, 0.0, // k, eps, v2, f
                                            0.0, 0.0, 0.0,      // muT, Tm, Lm
                                            3, 4,               // ndim, nvar
                                            NULL, mock_config); // constants, config

    mock_var_array[2] = new CHybridConvVariable(hybrid_soln, 3, 1, mock_config);

    mock_mediator = new CHybrid_Mediator(3, mock_config);
  }

  ~HybridRdeltaFixture() {
    delete mock_mediator;
    delete mock_var_array[2];
    delete mock_var_array[1];
    delete mock_var_array[0];
    delete mock_var_array;
    delete mock_config;
  }

  const su2double machine_eps;

  CConfig *mock_config;
  CVariable **mock_var_array;
  CHybrid_Mediator *mock_mediator;
};

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

#ifdef BUILD_TESTS

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

BOOST_FIXTURE_TEST_CASE(ZeroGradientTrivial, HybridRdeltaFixture) {

  mock_var_array[0]->SetGradient_PrimitiveZero(7);
  mock_var_array[0]->SetEddyViscosity(0.0);
  mock_mediator->ComputeInvLengthTensor(mock_var_array[0],
                                        mock_var_array[1],
                                        mock_var_array[2],
                                        mock_config->GetKind_Hybrid_Resolution_Indicator());

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

BOOST_FIXTURE_TEST_CASE(Shear_dudy, HybridRdeltaFixture) {

  //------------------------------------------------------
  // Simplest possible case:
  // du/dy = 1
  // mut = 1
  // v2 = 1
  // alpha = 1
  //------------------------------------------------------

  mock_var_array[0]->SetGradient_PrimitiveZero(7);
  mock_var_array[0]->SetGradient_Primitive(1, 1, 1.0);
  mock_var_array[0]->SetEddyViscosity(1.0);
  mock_var_array[1]->SetSolution(2, 1.0);  // v2
  mock_var_array[2]->SetSolution(0, 1.0);  // alpha

  mock_mediator->ComputeInvLengthTensor(mock_var_array[0],
                                        mock_var_array[1],
                                        mock_var_array[2],
                                        mock_config->GetKind_Hybrid_Resolution_Indicator());

  // 0,0 and 1,1 entries should be 0.5
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(0,0),0.5,machine_eps);
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(1,1),0.5,machine_eps);

  // everybody else is 0
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,2),0.0);


  //------------------------------------------------------
  // Check responds to mut correctly
  // du/dy = 1
  // mut = 2
  // v2 = 1
  // alpha = 1
  //------------------------------------------------------

  mock_var_array[0]->SetGradient_PrimitiveZero(7);
  mock_var_array[0]->SetGradient_Primitive(1, 1, 1.0);
  mock_var_array[0]->SetEddyViscosity(2.0);
  mock_var_array[1]->SetSolution(2, 1.0);  // v2
  mock_var_array[2]->SetSolution(0, 1.0);  // alpha

  mock_mediator->ComputeInvLengthTensor(mock_var_array[0],
                                        mock_var_array[1],
                                        mock_var_array[2],
                                        mock_config->GetKind_Hybrid_Resolution_Indicator());

  // 0,0 and 1,1 entries should be 1.0
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(0,0),1.0,machine_eps);
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(1,1),1.0,machine_eps);

  // everybody else is 0
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,2),0.0);

  //------------------------------------------------------
  // Check responds to v2 correctly
  // du/dy = 1
  // mut = 1
  // v2 = 4
  // alpha = 1
  //------------------------------------------------------

  mock_var_array[0]->SetGradient_PrimitiveZero(7);
  mock_var_array[0]->SetGradient_Primitive(1, 1, 1.0);
  mock_var_array[0]->SetEddyViscosity(1.0);
  mock_var_array[1]->SetSolution(2, 4.0);  // v2
  mock_var_array[2]->SetSolution(0, 1.0);  // alpha

  mock_mediator->ComputeInvLengthTensor(mock_var_array[0],
                                        mock_var_array[1],
                                        mock_var_array[2],
                                        mock_config->GetKind_Hybrid_Resolution_Indicator());

  // 0,0 and 1,1 entries should be 0.0625
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(0,0),0.0625,machine_eps);
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(1,1),0.0625,machine_eps);

  // everybody else is 0
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,2),0.0);

  //------------------------------------------------------
  // Check responds to alpha correctly
  // du/dy = 1
  // mut = 1
  // v2 = 1
  // alpha = 0.5
  //------------------------------------------------------

  mock_var_array[0]->SetGradient_PrimitiveZero(7);
  mock_var_array[0]->SetGradient_Primitive(1, 1, 1.0);
  mock_var_array[0]->SetEddyViscosity(1.0);
  mock_var_array[1]->SetSolution(2, 1.0);  // v2
  mock_var_array[2]->SetSolution(0, 0.25);  // alpha

  mock_mediator->ComputeInvLengthTensor(mock_var_array[0],
                                        mock_var_array[1],
                                        mock_var_array[2],
                                        mock_config->GetKind_Hybrid_Resolution_Indicator());

  // 0,0 and 1,1 entries should be 0.0625
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(0,0),0.25,machine_eps);
  BOOST_CHECK_CLOSE(mock_mediator->GetInvLengthScale(1,1),0.25,machine_eps);

  // everybody else is 0
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(0,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(1,2),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,0),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,1),0.0);
  BOOST_CHECK_EQUAL(mock_mediator->GetInvLengthScale(2,2),0.0);

}

BOOST_FIXTURE_TEST_CASE(PureRotation, HybridRdeltaFixture) {

  mock_var_array[0]->SetGradient_PrimitiveZero(7);


  // grad(u)
  mock_var_array[0]->SetGradient_Primitive(1, 0, 0.0 );
  mock_var_array[0]->SetGradient_Primitive(1, 1, 3.14);
  mock_var_array[0]->SetGradient_Primitive(1, 2, 1.59);

  // grad(v)
  mock_var_array[0]->SetGradient_Primitive(2, 0, -3.14);
  mock_var_array[0]->SetGradient_Primitive(2, 1,  0.0 );
  mock_var_array[0]->SetGradient_Primitive(2, 2, 42.0 );

  // grad(w)
  mock_var_array[0]->SetGradient_Primitive(3, 0,  -1.59);
  mock_var_array[0]->SetGradient_Primitive(3, 1, -42.0 );
  mock_var_array[0]->SetGradient_Primitive(3, 2,   0.0 );

  // mut = 1.0
  mock_var_array[0]->SetEddyViscosity(1.0);

  // v2 = 2.0
  mock_var_array[1]->SetSolution(2, 1.0);

  mock_mediator->ComputeInvLengthTensor(mock_var_array[0],
                                        mock_var_array[1],
                                        mock_var_array[2],
                                        mock_config->GetKind_Hybrid_Resolution_Indicator());

  // Everything should be zero
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

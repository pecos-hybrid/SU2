/*!
 * \file MPI_global_fixture.hpp
 * \brief This header defines a global fixture that allows use of MPI
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

#ifdef BUILD_TESTS
#include "boost/test/included/unit_test.hpp"

#ifdef HAVE_MPI
#include <mpi.h>
#endif

/*
 * After including this header, include the line:
 *
 * BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );
 *
 * to use this fixture in addition to any other fixtures.  This macro will
 * have file scope, affecting all tests and test suites in the current file.
 */

// A fixture to automate setup/teardown of MPI.
struct MPIGlobalFixture {

  /**
   * Setup MPI
   *
   * When this fixture is used, MPI will be setup whenever a test starts.
   */
  MPIGlobalFixture() {
#ifdef HAVE_MPI
  MPI_Init(NULL,NULL);
  BOOST_TEST_MESSAGE("Calling MPI_Init...");
#endif
  }

  /**
   * Tear down MPI
   *
   * When this fixture is used, MPI will be finalize whenever a test starts.
   * This will happen whether or not a particular test fails or throws an
   * exception.  If the test is expected to call MPI_Abort or MPI_Finalize
   * itself, you probably shouldn't use this fixture.
   */
  ~MPIGlobalFixture() {
#ifdef HAVE_MPI
  MPI_Finalize();
  BOOST_TEST_MESSAGE("Calling MPI_Finalize...");
#endif
  }

};

#endif

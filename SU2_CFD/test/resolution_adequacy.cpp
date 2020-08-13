/*!
 * \file resolution_adequacy.cpp
 * \brief Tests the resolution adequacy calculation
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

#include <random>
#include <limits> // used to find machine epsilon

#include "../include/solver_structure.hpp"
#include "../include/solver_structure_v2f.hpp"
#include "../include/variable_structure_v2f.hpp"
#include "../include/hybrid_RANS_LES_model.hpp"

namespace resolution_adequacy_test {

const unsigned short nDim = 3;
const unsigned short nVar = nDim+2;
const unsigned short nPrimVar = nDim+9;
const unsigned short nSecVar = 4;

/**
 * Write a cfg file to be used in initializing the CConfig object.
 */
static void WriteCfgFile(const char* filename) {

  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "HYBRID_RANSLES= MODEL_SPLIT" << std::endl;
  cfg_file << "RUNTIME_AVERAGING= POINTWISE" << std::endl;
  cfg_file << "UNSTEADY_SIMULATION= TIME_STEPPING" << std::endl;
  cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
  cfg_file.close();

}

class CTestGeometry : public CGeometry {
 public:
  CTestGeometry(CConfig* config) : CGeometry() {
    nPoint = 1;
    /*--- Local scope overrides file scope ---*/
    nDim = resolution_adequacy_test::nDim;

    node = new CPoint*[nPoint];
    for (unsigned short iPoint = 0; iPoint < nPoint; iPoint++) {
      node[iPoint] = new CPoint(0, 0, 0, iPoint, config);
      node[iPoint]->SetVolume(1);
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        for (unsigned short jDim = 0; jDim < nDim; jDim++) {
          node[iPoint]->SetResolutionTensor(iDim, jDim, (iDim == jDim));
        }
      }
    }
  }
};

class CTestFlowSolver : public CNSSolver {
 public:
  CTestFlowSolver(unsigned short val_nDim, unsigned short val_nVar,
                  CConfig* config)
       : CNSSolver() {

    nDim = val_nDim;
    nVar = val_nVar;
    nPoint = 1;
    nPointDomain = 1;

    Solution = new su2double[val_nVar];

    node = new CVariable*[nPoint];
    average_node = new CVariable*[nPoint];

    const su2double density = 1.0;
    const su2double laminar_viscosity = 1.0;
    su2double* velocity = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      velocity[iDim] = 0;
    }

    node[0] = new CNSVariable(density, velocity, 0, nDim, nVar, config);
    average_node[0] =
        new CNSVariable(density, velocity, 0, nDim, nVar, config);

    node[0]->SetLaminarViscosity(laminar_viscosity);
    average_node[0]->SetLaminarViscosity(laminar_viscosity);
    average_node[0]->SetKineticEnergyRatio(1.0);

    delete [] velocity;
  };
};

class CTestTurbSolver : public CTurbKESolver {
 public:
  CTestTurbSolver(su2double kine, su2double epsi, su2double v2,
                  unsigned short val_nDim, CConfig* config)
       : CTurbKESolver() {

    nDim = val_nDim;
    nVar = 4;
    nPoint = 1;
    nPointDomain = 1;

    Solution = new su2double[nVar];

    node = new CVariable*[nPoint];

    const su2double f = 1;
    const su2double muT = 1;
    const su2double Tm = 1;
    const su2double Lm = 1;

    node[0] = new CTurbKEVariable(kine, epsi, v2, f, muT, Tm, Lm,
                                  nDim, nVar, config);
  };
};

BOOST_AUTO_TEST_SUITE(ResolutionAdequacy);

struct RkFixture {
  RkFixture() {
    /*--- Setup ---*/

    char cfg_filename[100] = "resolution_adequacy.cfg";
    WriteCfgFile(cfg_filename);
    const unsigned short iZone = 0;
    const unsigned short nZone = 1;
    config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, 2, VERB_NONE);
    std::remove(cfg_filename);

    geometry = new CTestGeometry(config);

    solver_container = new CSolver*[TURB_SOL+1];
    solver_container[FLOW_SOL] =
          new CTestFlowSolver(nDim, nVar, config);

    flow_vars = solver_container[FLOW_SOL]->node[0];
    flow_avgs = solver_container[FLOW_SOL]->average_node[0];

    const su2double kine = 3.0;
    const su2double epsi = 1.0;
    const su2double v2 = 2.0;
    solver_container[TURB_SOL] =
          new CTestTurbSolver(kine, epsi, v2, nDim, config);
  };

  ~RkFixture() {
    delete solver_container[FLOW_SOL];
    delete solver_container[TURB_SOL];
    delete [] solver_container;
    delete geometry;
    delete config;
  };

  CConfig* config;
  CGeometry* geometry;
  CSolver** solver_container;
  CVariable* flow_vars;
  CVariable* flow_avgs;
};

/**
 */
BOOST_FIXTURE_TEST_CASE(MaxIsOneWhenAlphaGreaterThanOne, RkFixture) {

  /*--- Setup ---*/

  /*--- To force rk to be very large, set dUdy to be small and
   * M to be large ---*/
  const su2double very_large = 1E3;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    const unsigned short u_index = iDim % 3;
    const unsigned short x_index = (iDim+1) % 3;
    flow_vars->AddGradient_Primitive(u_index+1, x_index, very_large);
    flow_avgs->AddGradient_Primitive(u_index+1, x_index, very_large);
  }

  const unsigned short iPoint = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      geometry->node[iPoint]->SetResolutionTensor(iDim, jDim, very_large*(iDim == jDim));
    }
  }

  /*--- Set alpha to be an unrealistic value ---*/
  flow_avgs->SetKineticEnergyRatio(pow(1.1, 1.0/1.7));

  /*--- Test --*/

  CHybrid_Mediator mediator(nDim, config);
  mediator.ComputeResolutionAdequacy(geometry, solver_container, iPoint);
  const su2double r_k = flow_vars->GetResolutionAdequacy();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(r_k, su2double(1.0), tolerance);
}

/**
 */
BOOST_FIXTURE_TEST_CASE(MaxIs30WhenAlphaLessThanOne, RkFixture) {

  /*--- Setup ---*/

  /*--- To force rk to be very large, set dUdy to be large and
   * M to be large ---*/
  const su2double very_large = 1E3;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    const unsigned short u_index = iDim % 3;
    const unsigned short x_index = (iDim+1) % 3;
    flow_vars->AddGradient_Primitive(u_index+1, x_index, very_large);
    flow_avgs->AddGradient_Primitive(u_index+1, x_index, very_large);
  }

  const unsigned short iPoint = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      geometry->node[iPoint]->SetResolutionTensor(iDim, jDim, very_large*(iDim == jDim));
    }
  }

  /*--- Set alpha to be a realistic value ---*/
  flow_avgs->SetKineticEnergyRatio(pow(0.5, 1.0/1.7));

  /*--- Test --*/

  CHybrid_Mediator mediator(nDim, config);
  mediator.ComputeResolutionAdequacy(geometry, solver_container, iPoint);
  const su2double r_k = flow_vars->GetResolutionAdequacy();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(r_k, su2double(30.0), tolerance);
}

BOOST_FIXTURE_TEST_CASE(MinEnforcedWhenRkIsSmall, RkFixture) {

  /*--- Setup ---*/

  /*--- To force rk to be very small, leave dUdy as 0 ---*/
  /*--- Ensure isotropic contribution is zero --*/
  const unsigned short iPoint = 0;
  solver_container[TURB_SOL]->node[iPoint]->SetSolution(0, 0);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      geometry->node[iPoint]->SetResolutionTensor(iDim, jDim, (iDim == jDim));
    }
  }

  /*--- Set alpha to be a realistic value ---*/
  flow_avgs->SetKineticEnergyRatio(pow(0.5, 1.0/1.7));

  /*--- Test --*/

  CHybrid_Mediator mediator(nDim, config);
  mediator.ComputeResolutionAdequacy(geometry, solver_container, iPoint);
  const su2double r_k = flow_vars->GetResolutionAdequacy();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(r_k, su2double(1E-8), tolerance);
}

BOOST_FIXTURE_TEST_CASE(SimpleRkTest, RkFixture) {

  /*--- Setup ---*/

  const unsigned short iPoint = 0;
  flow_vars->AddGradient_Primitive(1, 1, 1.0);
  flow_avgs->AddGradient_Primitive(1, 1, 1.0);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      geometry->node[iPoint]->SetResolutionTensor(iDim, jDim, (iDim == jDim));
    }
  }

  /*--- Set v2 to 4/3 make (3/2*alpha*v2)^(1.5)=1.0 ---*/
  solver_container[TURB_SOL]->node[iPoint]->SetSolution(0, 0);
  solver_container[TURB_SOL]->node[iPoint]->SetSolution(2, 4.0/3);

  /*--- Set alpha to be a realistic value ---*/
  flow_avgs->SetKineticEnergyRatio(pow(0.5, 1.0/1.7));

  /*--- Test --*/

  CHybrid_Mediator mediator(nDim, config);
  mediator.ComputeResolutionAdequacy(geometry, solver_container, iPoint);
  const su2double r_k = flow_vars->GetResolutionAdequacy();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(r_k, su2double(0.75), tolerance);
}

BOOST_AUTO_TEST_SUITE_END();

} // end namespace

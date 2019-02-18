/*!
 * \file forcing_test.cpp
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

#define BOOST_TEST_MODULE Hybrid_Model
#include "boost/test/included/unit_test.hpp"
#include "MPI_global_fixture.hpp"

#include <limits>

#include "../include/hybrid_RANS_LES_forcing.hpp"
#include "../include/solver_structure.hpp"
#include "../include/variable_structure.hpp"

// Setup MPI
BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

void WriteCfgFile(const char* filename) {
  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
  cfg_file << "HYBRID_RANSLES= MODEL_SPLIT" << std::endl;
  cfg_file << "RUNTIME_AVERAGING= POINTWISE" << std::endl;
  cfg_file << "UNSTEADY_SIMULATION= TIME_STEPPING" << std::endl;
  cfg_file << "HYBRID_FORCING_PERIODIC_LENGTH= (";
  cfg_file << std::setprecision(std::numeric_limits<su2double>::digits10 + 1);
  cfg_file << M_PI << ", 0.25, 1.1780972451)" << std::endl;
  cfg_file.close();
}

/*--- Test Doubles ---*/

class TestGeometry : public CGeometry {
 public:
  TestGeometry(CConfig* config) : CGeometry() {
    nPoint = 1;
    nPointDomain = nPoint;
    node = new CPoint*[nPoint];
    unsigned short iPoint = 0;
    const su2double coord[3] = {0, 0, 0};
    node[iPoint] = new CPoint(coord[0], coord[1], coord[2], iPoint, config);
  }
};

class TestFlowVariable : public CVariable {
 private:
  su2double LaminarViscosity;
  su2double* Primitive;
 public:
  TestFlowVariable(unsigned short val_nDim, CConfig *config)
      : CVariable(val_nDim, 5, config) {

    Solution[0] = 1;
    Solution[1] = 23.0867;
    Solution[2] = -0.0728996;
    Solution[3] = -0.311579;

    nPrimVar = nVar;
    Primitive = new su2double[nPrimVar];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Primitive[iDim+1] = Solution[iDim+1]/Solution[0];
    }

    LaminarViscosity = 0.000192831;
    Delta_Time = 1E-3;
  }
  ~TestFlowVariable() {
    delete [] Primitive;
  }

  su2double GetLaminarViscosity() { return LaminarViscosity; }
  su2double* GetPrimitive() { return Primitive; }
};

class TestAvgVariable : public CVariable {
 private:
  su2double* Primitive;
  su2double KineticEnergyRatio, ResolutionAdequacy;
 public:
  TestAvgVariable(unsigned short val_nDim, CConfig *config)
      : CVariable(val_nDim, 5, config) {

    nPrimVar = nVar;
    Primitive = new su2double[nPrimVar];
    Primitive[1] = 22.8796;
    Primitive[2] = -0.266714;
    Primitive[3] = 0.00357933;
    ResolutionAdequacy = 0.8;
    KineticEnergyRatio = 0.598898;
  }
  ~TestAvgVariable() {
    delete [] Primitive;
  }

  su2double* GetPrimitive() { return Primitive; }
  su2double GetKineticEnergyRatio() const { return KineticEnergyRatio; }
  su2double GetResolutionAdequacy() const { return ResolutionAdequacy; }
};

class TestTurbVariable : public CVariable {
 private:
  su2double TurbTimescale, TurbLengthscale;
 public:
  TestTurbVariable(unsigned short val_nDim,
                   CConfig *config)
      : CVariable(val_nDim, 4, config) {

    Solution[0] = 2.28081;
    Solution[1] = 9.82732;
    Solution[2] = 1.22631;
  }
};

class TestSolver : public CSolver {
 public:
  TestSolver(unsigned short val_nPoint) : CSolver() {
    nPoint = val_nPoint;
    node = new CVariable*[nPoint];
    average_node = new CVariable*[nPoint];
  }
};

struct ForcingFixture {
  ForcingFixture() {
    nDim = 3;
    char cfg_filename[100] = "forcing_test.cfg";
    WriteCfgFile(cfg_filename);
    const unsigned short iZone = 0;
    const unsigned short nZone = 1;
    config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, nDim, VERB_NONE);
    std::remove(cfg_filename);

    config->SetCurrent_UnstTime(1.0);
    // Check that periodic length was set correctly in *.cfg file:
    su2double* periodic_length = config->GetHybrid_Forcing_Periodic_Length();
    assert(abs(periodic_length[0] - M_PI) < 1E-8);

    geometry = new TestGeometry(config);
    geometry->SetnDim(nDim);
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      for (unsigned short jDim = 0; jDim < nDim; jDim++) {
        geometry->node[0]->SetResolutionTensor(iDim, jDim, 0);
      }
    }
    geometry->node[0]->SetResolutionTensor(0, 0, 0.0634665);
    geometry->node[0]->SetResolutionTensor(1, 1, 0.0146755);
    geometry->node[0]->SetResolutionTensor(2, 2, 0.0785398);
    geometry->node[0]->SetWall_Distance(0.171989);

    solver = new CSolver*[MAX_SOLS];
    solver[FLOW_SOL] = new TestSolver(1);
    solver[TURB_SOL] = new TestSolver(1);
    solver[FLOW_SOL]->node[0] = new TestFlowVariable(nDim, config);
    solver[FLOW_SOL]->average_node[0] = new TestAvgVariable(nDim, config);
    solver[TURB_SOL]->node[0] = new TestTurbVariable(nDim, config);
    solver[TURB_SOL]->average_node[0] = new TestTurbVariable(nDim, config);
  }

  ~ForcingFixture() {
    delete geometry;
    delete solver[FLOW_SOL];
    delete solver[TURB_SOL];
    delete [] solver;
    delete config;
  }

  unsigned short nDim;
  CGeometry* geometry;
  CSolver** solver;
  CConfig* config;
};

/*--- PINCH POINT 1 ---*/
BOOST_FIXTURE_TEST_CASE(ForcingTest, ForcingFixture) {

  CHybridForcingTG0 forcing(geometry, config);
  forcing.ComputeForcingField(solver, geometry, config);
  const su2double* F = forcing.GetForcingVector(0);
  // Test values taken from CDP
  const su2double dt = 1E-3;
  const su2double true_F[3] = {-2.5810564961917728E-005/dt,
                               -1.0608578335743336E-004/dt,
                               -2.2531677802402187E-004/dt};

  const su2double tolerance = 1E-4;
  BOOST_CHECK_CLOSE_FRACTION(F[0], true_F[0], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(F[1], true_F[1], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(F[2], true_F[2], tolerance);

}

/*--- PINCH POINT 2 ---*/
BOOST_FIXTURE_TEST_CASE(ScalingCoefficient, ForcingFixture) {

  const su2double alpha = 0.598898;
  const su2double Ftar = 49.324157965434289;
  const su2double resolution_adequacy = 0.8;
  const su2double alpha_kol = 1.9086085598731258E-002;
  const su2double PFtest = 3.8390135906329317E-004;

  CHybridForcingTG0 forcing(geometry, config);
  forcing.ComputeForcingField(solver, geometry, config);
  const su2double eta = forcing.ComputeScalingFactor(Ftar, resolution_adequacy,
                                                     alpha, alpha_kol, PFtest);

  const su2double true_eta = 5.7950398607929117;
  const su2double tolerance = 1E-4;
  BOOST_CHECK_CLOSE_FRACTION(eta, true_eta, tolerance);
}

/*--- PINCH POINT 3 ---*/
BOOST_FIXTURE_TEST_CASE(TaylorGreenFields, ForcingFixture) {

  su2double x[3] = {22.879600000000000, -1.2667139999999999, 3.5793299999999999E-003};
  su2double Lsgs = 0.64981218416406406;
  su2double Lmesh[3] = {0.0634665, 0.0146755, 0.0785398};
  su2double D[3] = {M_PI, 0.25, 1.1780972451};
  su2double dwall = 0.171989;
  su2double h[3];

  CHybridForcingTG0 forcing(geometry, config);
  forcing.SetTGField(x, Lsgs, Lmesh, D, dwall, h);

  su2double true_h[3] = {-4.4539063719893331E-003,  -1.8306307791801466E-002,  -3.8880971216165662E-002};
  const su2double tolerance = 1E-4;
  BOOST_CHECK_CLOSE_FRACTION(h[0], true_h[0], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(h[1], true_h[1], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(h[2], true_h[2], tolerance);
}

/*--- PINCH POINT 4 ---*/
BOOST_FIXTURE_TEST_CASE(TargetForcing, ForcingFixture) {

  su2double v2 = 1.22631;
  su2double Ttot = 0.23208870780640092;
  su2double alpha = 0.598898;
  su2double Tsgs = Ttot*alpha;

  CHybridForcingTG0 forcing(geometry, config);
  const su2double F_target = forcing.GetTargetProduction(v2, Tsgs, alpha);

  const su2double true_F_target = 49.324157965434289;
  const su2double tolerance = 1E-4;
  BOOST_CHECK_CLOSE_FRACTION(F_target, true_F_target, tolerance);

}

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

#include "../include/numerics_structure.hpp"
#include "../include/numerics_structure_v2f.hpp"
#include "../include/solver_structure.hpp"
#include "../include/solver_structure_v2f.hpp"

// Setup MPI
BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

void WriteCfgFile(const char* filename) {
  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
  cfg_file.close();
}

/*--- Test Doubles ---*/

class TestGeometry : public CGeometry {
 public:
  TestGeometry(CConfig* config) : CGeometry() {
    nPoint = 1;
    nPointDomain = nPoint;
    nPointNode = nPoint;
    node = new CPoint*[nPoint];
    unsigned short iPoint = 0;
    const su2double coord[3] = {0, 0, 0};
    node[iPoint] = new CPoint(coord[0], coord[1], coord[2], iPoint, config);
    node[iPoint]->SetVolume(1.0);
  }
};

class TestFlowVariable : public CVariable {
 private:
  su2double LaminarViscosity;
  su2double* Primitive;
  su2double** PrimVarGrad;
  su2double StrainMag;
  su2double Density;
 public:
  TestFlowVariable(unsigned short val_nDim, CConfig *config)
      : CVariable(val_nDim, 5, config) {


    Delta_Time = 1E-3;
    StrainMag = 1.0;
    LaminarViscosity = 1E-5;
    Density = 1.0;

    Solution[0] = Density;
    Solution[1] = 20;
    Solution[2] = 0;
    Solution[3] = 0;

    nPrimVar = nDim+7;
    Primitive = new su2double[nPrimVar];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      Primitive[iDim+1] = Solution[iDim+1]/Density;
    }
    Primitive[nDim+2] = Density; // density
    Primitive[nDim+5] = LaminarViscosity; // laminar_viscosity
    Primitive[nDim+6] = 1E-3; // eddy_viscosity;

    PrimVarGrad = new su2double*[nPrimVar];
    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
      PrimVarGrad[iVar] = new su2double[nDim];
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        PrimVarGrad[iVar][iDim] = 0;
      }
    }
    PrimVarGrad[1][1] = -1; // dU/dy
  }

  ~TestFlowVariable() {
    delete [] Primitive;
    for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
      delete [] PrimVarGrad[iVar];
    }
    delete [] PrimVarGrad;
  }

  su2double* GetPrimitive() { return Primitive; }
  su2double** GetGradient_Primitive() { return PrimVarGrad; }
  su2double* GetVorticity() { return NULL; }
  su2double GetStrainMag() { return StrainMag; }
  su2double GetDensity() { return Density; }
  su2double GetLaminarViscosity() { return LaminarViscosity; }
};

class TestSolver : public CSolver {
 public:
  TestSolver(unsigned short val_nPoint) : CSolver() {
    nPoint = val_nPoint;
    node = new CVariable*[nPoint];
    average_node = NULL;
  }
};

BOOST_AUTO_TEST_CASE(ForcingTest) {

  /*--- ARRANGE --*/

  unsigned short nDim = 3;

  char cfg_filename[100] = "turb_scales_test.cfg";
  WriteCfgFile(cfg_filename);
  const unsigned short iZone = 0;
  const unsigned short nZone = 1;
  CConfig* config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, nDim, VERB_NONE);
  std::remove(cfg_filename);

  config->SetDensity_FreeStreamND(1.0);
  config->SetVelocity_FreeStreamND(20, 0);
  config->SetVelocity_FreeStreamND(0, 1);
  config->SetVelocity_FreeStreamND(0, 2);
  config->SetViscosity_FreeStreamND(1E-5);

  CGeometry* geometry = new TestGeometry(config);
  geometry->SetnDim(nDim);

  const unsigned short iMesh = 0;
  CSolver** solver = new CSolver*[MAX_SOLS];
  solver[FLOW_SOL] = new TestSolver(1);
  solver[FLOW_SOL]->node[0] = new TestFlowVariable(nDim, config);
  solver[TURB_SOL] = new CTurbKESolver(geometry, config, iMesh);

  su2double* constants = solver[TURB_SOL]->GetConstants();
  CNumerics *numerics = new CSourcePieceWise_TurbKE(nDim, 4, constants, config);
  CNumerics *second_numerics = NULL;

  /*--- ACT ---*/
  const unsigned short iRKStep = 0;
  solver[TURB_SOL]->Source_Residual(geometry, solver, numerics,
                                    second_numerics, config, iMesh, iRKStep);

  /*--- ASSERT ---*/
  const su2double tol=1E-8;
  // Values calculated before making changes:
  const su2double true_residual[4] = {1742.3990000000001,
                                    3845544.7801321908,
                                    1161.5993333333327,
                                    21488354022.044903};
  BOOST_CHECK_CLOSE_FRACTION(solver[TURB_SOL]->LinSysRes[0], true_residual[0], tol);
  BOOST_CHECK_CLOSE_FRACTION(solver[TURB_SOL]->LinSysRes[1], true_residual[1], tol);
  BOOST_CHECK_CLOSE_FRACTION(solver[TURB_SOL]->LinSysRes[2], true_residual[2], tol);
  BOOST_CHECK_CLOSE_FRACTION(solver[TURB_SOL]->LinSysRes[3], true_residual[3], tol);

  /*--- TEARDOWN ---*/

  delete numerics;
  delete geometry;
  delete solver[FLOW_SOL];
  delete solver[TURB_SOL];
  delete [] solver;
  delete config;
}

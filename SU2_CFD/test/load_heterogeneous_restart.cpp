/*!
 * \file
 * \brief
 * \author C. Pederson
 * \version
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

#define BOOST_TEST_MODULE LoadRestart
#include "MPI_global_fixture.hpp"

#include <iostream>
#include <string>

#include "../include/solver_structure.hpp"
#include "../include/solver_structure_v2f.hpp"

// Setup MPI
BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

/*--- Constants used in all of the tests ---*/

const unsigned short nDim = 3;
const unsigned short nVar = nDim + 2;

const su2double density = 1;
const su2double momentum[] = {1.0, 0, 0};
const su2double energy = 290930;

class TestGeometry : public CGeometry {
 public:
  unsigned long GetGlobal_nPointDomain() { return 1; };
  long GetGlobal_to_Local_Point(unsigned long val_ipoint) { return 0; };
};

class TestSolver : public CNSSolver {
 public:
  TestSolver(unsigned short val_nDim, unsigned short val_nVar, CConfig* config)
       : CNSSolver() {

    nDim = val_nDim;
    nVar = val_nVar;
    nPoint = 1;
    nPointDomain = 1;

    Solution = new su2double[val_nVar];

    node = new CVariable*[nPoint];
    average_node = new CVariable*[nPoint];

    Density_Inf = 0;
    Velocity_Inf = new su2double[nDim];
    for (unsigned short iDim = 0; iDim < nDim; iDim++) Velocity_Inf[iDim] = 0;
    Energy_Inf = 0;

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++)
      average_node[iPoint] = new CNSVariable(Density_Inf, Velocity_Inf, Energy_Inf, nDim, nVar, config);
  };

  void SetRestart_Data(passivedouble* restart_data) { Restart_Data = restart_data; };
  void SetRestart_Vars(int* restart_vars) { Restart_Vars = restart_vars; };
};

struct HetergeneousRestartFixture {
  HetergeneousRestartFixture() {

    geometry = new CGeometry*[1];
    geometry[0] = new TestGeometry();
    geometry[0]->SetnDim(nDim);

    char cfg_filename[100] = "load_heterogeneous_restart.cfg";
    WriteCfgFile(nDim, cfg_filename);
    const unsigned short iZone = 0;
    const unsigned short nZone = 1;
    config = new CConfig(cfg_filename, SU2_CFD, iZone, nZone, nDim, VERB_NONE);
    std::remove(cfg_filename);

    solver = new TestSolver(nDim, nVar, config);
  }

  void WriteCfgFile(unsigned short nDim, const char* filename) {
    std::ofstream cfg_file;

    cfg_file.open(filename, ios::out);
    cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
    cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
    cfg_file << "HYBRID_RANSLES= MODEL_SPLIT" << std::endl;
    cfg_file << "RUNTIME_AVERAGING= POINTWISE" << std::endl;
    cfg_file << "UNSTEADY_SIMULATION= TIME_STEPPING" << std::endl;
    cfg_file.close();
  }

  void SetupRestart(char** arr, size_t count, passivedouble* restart_data) {
    config->fields.clear();
    config->fields.insert(config->fields.begin(), arr, arr+count);
    const int nFields = (int)config->fields.size() - 1;
    /*--- Restart vars isn't actually a 2 element array, but we only need [1] ---*/
    restart_vars[1] = nFields;

    solver->SetRestart_Data(restart_data);
    solver->SetRestart_Vars(restart_vars);
  }

  ~HetergeneousRestartFixture() {
    solver->SetRestart_Data(NULL);
    solver->SetRestart_Vars(NULL);

//    delete geometry[1];
//    delete [] geometry;
    delete solver;
    delete config;
  }

  int restart_vars[2];
  CConfig* config;
  CGeometry** geometry;
  TestSolver* solver;
};

BOOST_FIXTURE_TEST_CASE(ResolvedSolutionLoadsFromResolvedRestart,
                        HetergeneousRestartFixture) {

  char* name_array[] = {"Point_ID",
      "\"x\"", "\"y\"", "\"z\"",
      "\"Density\"", "\"X-Momentum\"", "\"Y-Momentum\"", "\"Z-Momentum\"", "\"Energy\"",
      "\"TKE\"", "\"Dissipation\"", "\"v2\"", "\"f\"",
      "\"Pressure\"", "\"Temperature\"", "\"Mach\"",
      "\"Pressure_Coefficient\"", "\"Laminar_Viscosity\"",
      "\"Skin_Friction_Coefficient_X\"", "\"Skin_Friction_Coefficient_Y\"",
      "\"Skin_Friction_Coefficient_Z\"", "\"Heat_Flux\"", "\"Y_Plus\"",
      "\"Eddy_Viscosity\"", "\"L_m\"", "\"T_m\""};

  passivedouble restart_data[] =
        {/* x */ 0, /* y */ 0, /* z */ 0,
         density, momentum[0], momentum[1], momentum[2], energy,
         /* TKE */ 0,
         /* Dissipation */ 0.17797,
         /* v2 */ 0,
         /* f */ 0,
         /* Other variables */ 116372, 303.772, 0, 0.598286, 1.86371e-05,
         1.1907e-06, 0, -9.30875e-15, -0.286451, 1276.33, 3.76891e-31, 0.0246201,
         0.0531493};

  const size_t nFields = sizeof(name_array)/sizeof(name_array[0]);
  SetupRestart(name_array, nFields, restart_data);
  solver->LoadSolution(false, "dummy_string", config, geometry);
  su2double* solution = solver->node[0]->GetSolution();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(solution[0], density, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(solution[1], momentum[0], tolerance);
  BOOST_CHECK_SMALL(solution[2], tolerance);
  BOOST_CHECK_SMALL(solution[3], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(solution[4], energy, tolerance);
}

BOOST_FIXTURE_TEST_CASE(AverageSolutionLoadsFromAverageRestart,
                        HetergeneousRestartFixture) {

  char* name_array[] = {"Point_ID",
      "\"x\"", "\"y\"", "\"z\"",
      "\"Density\"", "\"X-Momentum\"", "\"Y-Momentum\"", "\"Z-Momentum\"", "\"Energy\"",
      "\"TKE\"", "\"Dissipation\"", "\"v2\"", "\"f\"",
      "\"Average_Density\"","\"Average_X-Momentum\"", "\"Average_Y-Momentum\"", "\"Average_Z-Momentum\"", "\"Average_Energy\"",
      "\"Pressure\"", "\"Temperature\"", "\"Mach\"",
      "\"Pressure_Coefficient\"", "\"Laminar_Viscosity\"",
      "\"Skin_Friction_Coefficient_X\"", "\"Skin_Friction_Coefficient_Y\"",
      "\"Skin_Friction_Coefficient_Z\"", "\"Heat_Flux\"", "\"Y_Plus\"",
      "\"Eddy_Viscosity\"", "\"L_m\"", "\"T_m\"",
      "\"alpha\"", "\"k_res\"",
      "\"tau_res_11\"", "\"tau_res_12\"", "\"tau_res_13\"", "\"tau_res_21\"",
      "\"tau_res_22\"", "\"tau_res_23\"", "\"tau_res_31\"", "\"tau_res_32\"",
      "\"tau_res_33\"", "\"mu_SGET_11\"", "\"mu_SGET_12\"", "mu_SGET_13\"",
      "\"mu_SGET_21\"", "\"mu_SGET_22\"", "\"mu_SGET_23\"", "\"mu_SGET_31\"",
      "\"mu_SGET_32\"", "\"mu_SGET_33\""};

  passivedouble restart_data[] =
        {/* x */ 0, /* y */ 0, /* z */ 0,
         /* resolved state */ density, momentum[0], momentum[1], momentum[2], energy,
         /* TKE */ 0,
         /* Dissipation */ 0.17797,
         /* v2 */ 0,
         /* f */ 0,
         /* average state */ density, momentum[0], momentum[1], momentum[2], energy,
         /* Other variables */ 116372, 303.772, 0, 0.598286, 1.86371e-05,
         1.1907e-06, 0, -9.30875e-15, -0.286451, 1276.33, 3.76891e-31, 0.0246201,
         0.0531493,
         /* alpha */ 1,
         /* k_res */ 0,
         /* tau_res */ 0, 0, 0, 0, 0, 0, 0, 0, 0,
         /* mu_SGET */ 0, 0, 0, 0, 0, 0, 0, 0, 0};

  const size_t nFields = sizeof(name_array)/sizeof(name_array[0]);
  SetupRestart(name_array, nFields, restart_data);
  solver->LoadSolution(false, "dummy_string", config, geometry);
  su2double* solution = solver->average_node[0]->GetSolution();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(solution[0], density, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(solution[1], momentum[0], tolerance);
  BOOST_CHECK_SMALL(solution[2], tolerance);
  BOOST_CHECK_SMALL(solution[3], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(solution[4], energy, tolerance);
}

BOOST_FIXTURE_TEST_CASE(AverageSolutionLoadsFromResolvedRestart,
                        HetergeneousRestartFixture) {

  char* name_array[] = {"Point_ID",
      "\"x\"", "\"y\"", "\"z\"",
      "\"Density\"", "\"X-Momentum\"", "\"Y-Momentum\"", "\"Z-Momentum\"", "\"Energy\"",
      "\"TKE\"", "\"Dissipation\"", "\"v2\"", "\"f\"",
      "\"Pressure\"", "\"Temperature\"", "\"Mach\"",
      "\"Pressure_Coefficient\"", "\"Laminar_Viscosity\"",
      "\"Skin_Friction_Coefficient_X\"", "\"Skin_Friction_Coefficient_Y\"",
      "\"Skin_Friction_Coefficient_Z\"", "\"Heat_Flux\"", "\"Y_Plus\"",
      "\"Eddy_Viscosity\"", "\"L_m\"", "\"T_m\""};

  passivedouble restart_data[] =
        {/* x */ 0, /* y */ 0, /* z */ 0,
         density, momentum[0], momentum[1], momentum[2], energy,
         /* TKE */ 0,
         /* Dissipation */ 0.17797,
         /* v2 */ 0,
         /* f */ 0,
         /* Other variables */ 116372, 303.772, 0, 0.598286, 1.86371e-05,
         1.1907e-06, 0, -9.30875e-15, -0.286451, 1276.33, 3.76891e-31, 0.0246201,
         0.0531493};

  const size_t nFields = sizeof(name_array)/sizeof(name_array[0]);
  SetupRestart(name_array, nFields, restart_data);
  solver->LoadSolution(false, "dummy_string", config, geometry);
  su2double* solution = solver->average_node[0]->GetSolution();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(solution[0], density, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(solution[1], momentum[0], tolerance);
  BOOST_CHECK_SMALL(solution[2], tolerance);
  BOOST_CHECK_SMALL(solution[3], tolerance);
  BOOST_CHECK_CLOSE_FRACTION(solution[4], energy, tolerance);
}

BOOST_FIXTURE_TEST_CASE(HybridSolutionLoadsProperlyFromRANSRestart,
                        HetergeneousRestartFixture) {

  char* name_array[] = {"Point_ID",
      "\"x\"", "\"y\"", "\"z\"",
      "\"Density\"", "\"X-Momentum\"", "\"Y-Momentum\"", "\"Z-Momentum\"", "\"Energy\"",
      "\"TKE\"", "\"Dissipation\"", "\"v2\"", "\"f\"",
      "\"Pressure\"", "\"Temperature\"", "\"Mach\"",
      "\"Pressure_Coefficient\"", "\"Laminar_Viscosity\"",
      "\"Skin_Friction_Coefficient_X\"", "\"Skin_Friction_Coefficient_Y\"",
      "\"Skin_Friction_Coefficient_Z\"", "\"Heat_Flux\"", "\"Y_Plus\"",
      "\"Eddy_Viscosity\"", "\"L_m\"", "\"T_m\""};

  passivedouble restart_data[] =
        {/* x */ 0, /* y */ 0, /* z */ 0,
         density, momentum[0], momentum[1], momentum[2], energy,
         /* TKE */ 0,
         /* Dissipation */ 0.17797,
         /* v2 */ 0,
         /* f */ 0,
         /* Other variables */ 116372, 303.772, 0, 0.598286, 1.86371e-05,
         1.1907e-06, 0, -9.30875e-15, -0.286451, 1276.33, 3.76891e-31, 0.0246201,
         0.0531493};

  const size_t nFields = sizeof(name_array)/sizeof(name_array[0]);
  SetupRestart(name_array, nFields, restart_data);
  solver->LoadSolution(false, "dummy_string", config, geometry);
  solver->InitAverages();

  const su2double k_resolved = solver->average_node[0]->GetResolvedKineticEnergy();
  su2double** resolved_stress = solver->average_node[0]->GetResolvedTurbStress();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_SMALL(k_resolved, tolerance);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      BOOST_CHECK_SMALL(resolved_stress[iDim][jDim], tolerance);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(HybridSolutionLoadsFromHybridRestart,
                        HetergeneousRestartFixture) {

  char* name_array[] = {"Point_ID",
      "\"x\"", "\"y\"", "\"z\"",
      "\"Density\"", "\"X-Momentum\"", "\"Y-Momentum\"", "\"Z-Momentum\"", "\"Energy\"",
      "\"TKE\"", "\"Dissipation\"", "\"v2\"", "\"f\"",
      "\"Average_Density\"","\"Average_X-Momentum\"", "\"Average_Y-Momentum\"", "\"Average_Z-Momentum\"", "\"Average_Energy\"",
      "\"Pressure\"", "\"Temperature\"", "\"Mach\"",
      "\"Pressure_Coefficient\"", "\"Laminar_Viscosity\"",
      "\"Skin_Friction_Coefficient_X\"", "\"Skin_Friction_Coefficient_Y\"",
      "\"Skin_Friction_Coefficient_Z\"", "\"Heat_Flux\"", "\"Y_Plus\"",
      "\"Eddy_Viscosity\"", "\"L_m\"", "\"T_m\"",
      "\"alpha\"", "\"k_res\"",
      "\"tau_res_11\"", "\"tau_res_12\"", "\"tau_res_13\"", "\"tau_res_21\"",
      "\"tau_res_22\"", "\"tau_res_23\"", "\"tau_res_31\"", "\"tau_res_32\"",
      "\"tau_res_33\"", "\"mu_SGET_11\"", "\"mu_SGET_12\"", "mu_SGET_13\"",
      "\"mu_SGET_21\"", "\"mu_SGET_22\"", "\"mu_SGET_23\"", "\"mu_SGET_31\"",
      "\"mu_SGET_32\"", "\"mu_SGET_33\""};

  const su2double resolved_stress = -4;
  const su2double k_resolved = -0.5*(3*resolved_stress)/density;
  passivedouble restart_data[] =
        {/* x */ 0, /* y */ 0, /* z */ 0,
         /* resolved state */ density, momentum[0], momentum[1], momentum[2], energy,
         /* TKE */ 0,
         /* Dissipation */ 0.17797,
         /* v2 */ 0,
         /* f */ 0,
         /* average state */ density, momentum[0], momentum[1], momentum[2], energy,
         /* Other variables */ 116372, 303.772, 0, 0.598286, 1.86371e-05,
         1.1907e-06, 0, -9.30875e-15, -0.286451, 1276.33, 3.76891e-31, 0.0246201,
         0.0531493,
         /* alpha */ 1,
         /* k_res */ k_resolved,
         /* tau_res */ resolved_stress, resolved_stress, resolved_stress,
         resolved_stress, resolved_stress, resolved_stress,
         resolved_stress, resolved_stress, resolved_stress,
         /* mu_SGET */ 0, 0, 0, 0, 0, 0, 0, 0, 0};


  const size_t nFields = sizeof(name_array)/sizeof(name_array[0]);
  SetupRestart(name_array, nFields, restart_data);
  solver->LoadSolution(false, "dummy_string", config, geometry);

  const su2double k_resolved_actual =
      solver->average_node[0]->GetResolvedKineticEnergy();
  su2double** resolved_stress_actual =
      solver->average_node[0]->GetResolvedTurbStress();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE(k_resolved_actual, k_resolved, tolerance);
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      BOOST_CHECK_CLOSE_FRACTION(resolved_stress_actual[iDim][jDim],
                                 resolved_stress, tolerance);
    }
  }
}

BOOST_FIXTURE_TEST_CASE(HybridSolutionLoadsProductionFromHybridRestart,
                       HetergeneousRestartFixture) {

  char* name_array[] = {"Point_ID",
      "\"x\"", "\"y\"", "\"z\"",
      "\"Density\"", "\"X-Momentum\"", "\"Y-Momentum\"", "\"Z-Momentum\"", "\"Energy\"",
      "\"TKE\"", "\"Dissipation\"", "\"v2\"", "\"f\"",
      "\"Average_Density\"","\"Average_X-Momentum\"", "\"Average_Y-Momentum\"", "\"Average_Z-Momentum\"", "\"Average_Energy\"",
      "\"Pressure\"", "\"Temperature\"", "\"Mach\"",
      "\"Pressure_Coefficient\"", "\"Laminar_Viscosity\"",
      "\"Skin_Friction_Coefficient_X\"", "\"Skin_Friction_Coefficient_Y\"",
      "\"Skin_Friction_Coefficient_Z\"", "\"Heat_Flux\"", "\"Y_Plus\"",
      "\"Eddy_Viscosity\"", "\"L_m\"", "\"T_m\"",
      "\"alpha\"", "\"k_res\"", "\"Production\"", "\"mu_SGET_11\"", "\"mu_SGET_12\"", "mu_SGET_13\"",
      "\"mu_SGET_21\"", "\"mu_SGET_22\"", "\"mu_SGET_23\"", "\"mu_SGET_31\"",
      "\"mu_SGET_32\"", "\"mu_SGET_33\""};

  const su2double production = 4;
  const su2double k_resolved = 0.5;
  passivedouble restart_data[] =
        {/* x */ 0, /* y */ 0, /* z */ 0,
         /* resolved state */ density, momentum[0], momentum[1], momentum[2], energy,
         /* TKE */ 0,
         /* Dissipation */ 0.17797,
         /* v2 */ 0,
         /* f */ 0,
         /* average state */ density, momentum[0], momentum[1], momentum[2], energy,
         /* Other variables */ 116372, 303.772, 0, 0.598286, 1.86371e-05,
         1.1907e-06, 0, -9.30875e-15, -0.286451, 1276.33, 3.76891e-31, 0.0246201,
         0.0531493,
         /* alpha */ 1,
         /* k_res */ k_resolved,
         /* improved production */ production,
         /* mu_SGET */ 0, 0, 0, 0, 0, 0, 0, 0, 0};


  const size_t nFields = sizeof(name_array)/sizeof(name_array[0]);
  SetupRestart(name_array, nFields, restart_data);
  solver->LoadSolution(false, "dummy_string", config, geometry);

  const su2double actual_production = solver->average_node[0]->GetProduction();
  const su2double actual_k_resolved = solver->average_node[0]->GetResolvedKineticEnergy();

  const su2double tolerance = 10*std::numeric_limits<su2double>::epsilon();
  BOOST_CHECK_CLOSE_FRACTION(actual_production, production, tolerance);
  BOOST_CHECK_CLOSE_FRACTION(actual_k_resolved, k_resolved, tolerance);
}

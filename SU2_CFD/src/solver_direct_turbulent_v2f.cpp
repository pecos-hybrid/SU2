/*!
 * \file solution_direct_turbulent.cpp
 * \brief Main subrotuines for solving direct problems
 * \author F. Palacios, A. Bueno
 * \version 6.0.1 "Falcon"
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

#include "../include/solver_structure.hpp"
#include "../include/numerics_structure_v2f.hpp"
#include "../include/solver_structure_v2f.hpp"
#include "../include/variable_structure_v2f.hpp"

CTurbKESolver::CTurbKESolver(void) : CTurbSolver() {

  /*--- Array initialization ---*/
  constants = NULL;

}

CTurbKESolver::CTurbKESolver(CGeometry *geometry, CConfig *config,
                             unsigned short iMesh)
    : CTurbSolver(geometry, config) {

  unsigned short iVar, iDim, nLineLets;
  unsigned long iPoint;
  ifstream restart_file;
  string text_line;
  const bool runtime_averaging = (config->GetKind_Averaging() != NO_AVERAGING);

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  /*--- Array initialization ---*/
  constants = NULL;

  Gamma = config->GetGamma();
  Gamma_Minus_One = Gamma - 1.0;

  /*--- Dimension of the problem --> dependent of the turbulent model ---*/
  nVar = 4;

  nPoint = geometry->GetnPoint();
  nPointDomain = geometry->GetnPointDomain();

  /*--- Initialize nVarGrad for deallocation ---*/

  nVarGrad = nVar;

  /*--- Define geometry constants in the solver structure ---*/
  nDim = geometry->GetnDim();
  node = new CVariable*[nPoint];
  if (runtime_averaging) {
    average_node = new CVariable*[nPoint];
  }

  /*--- Single grid simulation ---*/
  if (iMesh == MESH_0) {

    /*--- Define some auxiliary vector related with the residual ---*/
    Residual     = new su2double[nVar];
    Residual_RMS = new su2double[nVar];
    Residual_i   = new su2double[nVar];
    Residual_j   = new su2double[nVar];
    Residual_Max = new su2double[nVar];

    for (iVar=0; iVar<nVar; iVar++) {
       Residual[iVar]     = 0.0;
       Residual_RMS[iVar] = 0.0;
       Residual_i[iVar]   = 0.0;
       Residual_j[iVar]   = 0.0;
       Residual_Max[iVar] = 0.0;
    }

    /*--- Define some structures for locating max residuals ---*/
    Point_Max = new unsigned long[nVar];
    for (iVar = 0; iVar < nVar; iVar++) Point_Max[iVar] = 0;

    Point_Max_Coord = new su2double*[nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Point_Max_Coord[iVar] = new su2double[nDim];
      for (iDim = 0; iDim < nDim; iDim++) Point_Max_Coord[iVar][iDim] = 0.0;
    }

    /*--- Define some auxiliary vector related with the solution ---*/
    Solution = new su2double[nVar];
    Solution_i = new su2double[nVar]; Solution_j = new su2double[nVar];

    /*--- Define some auxiliary vector related with the geometry ---*/
    Vector_i = new su2double[nDim]; Vector_j = new su2double[nDim];

    /*--- Define some auxiliary vector related with the flow solution ---*/
    FlowPrimVar_i = new su2double [nDim+9]; FlowPrimVar_j = new su2double [nDim+9];

    /*--- Jacobians and vector structures for implicit computations ---*/
    Jacobian_i = new su2double* [nVar];
    Jacobian_j = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Jacobian_i[iVar] = new su2double [nVar];
      Jacobian_j[iVar] = new su2double [nVar];
    }

    /*--- Initialization of the structure of the whole Jacobian ---*/
    if (rank == MASTER_NODE) {
      cout << "Initialize Jacobian structure (KE model)." << endl;
    }

    Jacobian.Initialize(nPoint, nPointDomain, nVar, nVar,
                        true, geometry, config);

    if (rank == MASTER_NODE) cout << "Finished." << endl;

    if ((config->GetKind_Linear_Solver_Prec() == LINELET) ||
        (config->GetKind_Linear_Solver() == SMOOTHER_LINELET)) {
      nLineLets = Jacobian.BuildLineletPreconditioner(geometry, config);
      if (rank == MASTER_NODE) {
        cout << "Compute linelet structure. " << nLineLets
             << " elements in each line (average)." << endl;
      }
    }

    LinSysSol.Initialize(nPoint, nPointDomain, nVar, 0.0);
    LinSysRes.Initialize(nPoint, nPointDomain, nVar, 0.0);
  }

  /*--- Computation of gradients by least squares ---*/
  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/
    Smatrix = new su2double* [nDim];
    for (iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];

    /*--- c vector := transpose(WA)*(Wb) ---*/
    Cvector = new su2double* [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
    Cvector[iVar] = new su2double [nDim];
  }

  /*--- Initialize value for model constants ---*/
  constants = new su2double[11];

  /* v2-f */
  constants[0]  = 0.22;   //C_mu
  constants[1]  = 1.0/1.0;    //1/sigma_k
  constants[2]  = 1.0/1.3;    //1/sigma_e
  constants[3]  = 1.0/1.0;    //1/sigma_z
  constants[4]  = 1.4;    //C_e1^o
  constants[5]  = 1.9;    //C_e2
  constants[6]  = 1.4;    //C_1
  constants[7]  = 0.3;   //C_2p
  constants[8]  = 6.0;    //C_T
  constants[9]  = 0.23;   //C_L
  constants[10] = 70.0;   //C_eta

  /*--- Initialize lower and upper limits---*/
  lowerlimit = new su2double[nVar];
  upperlimit = new su2double[nVar];


  // These are used in the AddConservativeSolution calls
  // k
  lowerlimit[0] = -1.0e10;
  upperlimit[0] =  1.0e10;

  // epsi
  lowerlimit[1] = -1.0e10;
  upperlimit[1] =  1.0e10;

  // v2
  lowerlimit[2] = -1.0e10;
  upperlimit[2] =  1.0e10;

  // f
  lowerlimit[3] = -1.0e10;
  upperlimit[3] =  1.0e10;


  /*--- Flow infinity initialization stuff ---*/
  su2double rhoInf, *VelInf, muLamInf, Intensity, viscRatio, muT_Inf, Tm_Inf, Lm_Inf;
  rhoInf    = config->GetDensity_FreeStreamND();
  VelInf    = config->GetVelocity_FreeStreamND();
  muLamInf  = config->GetViscosity_FreeStreamND();
  Intensity = config->GetTurbulenceIntensity_FreeStream();
  viscRatio = config->GetTurb2LamViscRatio_FreeStream();
  
  su2double VelMag = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    VelMag += VelInf[iDim]*VelInf[iDim];
  }

  VelMag = sqrt(VelMag);

  su2double L_Inf = config->GetLength_Reynolds();
  su2double scale = 1.0e-14;
  su2double scalar_min = scale/(VelMag*VelMag);
  su2double tke_min = scalar_min*VelMag*VelMag;
  su2double tdr_min = scalar_min*pow(VelMag,3.0)/L_Inf;


  // Freestream eddy visc
  muT_Inf = muLamInf*viscRatio;

  // Convenience: freestream kinematic viscosities
  const su2double nuInf  = muLamInf/rhoInf;
  const su2double nutInf = muT_Inf /rhoInf;

  // Freestream TKE
  kine_Inf = 1.5*(VelMag*VelMag*Intensity*Intensity);

  // Freestream dissipation
  epsi_Inf = (2.0/3.0)*constants[0]*(kine_Inf*kine_Inf)/nutInf;
  const su2double ktmp = 2.0/3.0*constants[0]*constants[8]*kine_Inf/viscRatio;
  const su2double epsi_Inf_alt = ktmp*ktmp/nuInf;
  epsi_Inf = min( epsi_Inf, epsi_Inf_alt );

  // Fresstream v2
  zeta_Inf = 2.0/3.0*kine_Inf;


  // Freestream time scale
  Tm_Inf = kine_Inf/max(epsi_Inf,tdr_min);
  su2double Tkol_inf = constants[8]*sqrt(nuInf/max(epsi_Inf,tdr_min));
  Tm_Inf = max( Tm_Inf, Tkol_inf );

  // Freestream length scale
  Lm_Inf = pow(kine_Inf,1.5)/max(epsi_Inf,tdr_min);
  const su2double nu3 = nuInf*nuInf*nuInf;
  const su2double Lkol_Inf = constants[10]*pow(nu3/max(epsi_Inf,tdr_min),0.25);
  Lm_Inf = constants[9] * max( Lm_Inf, Lkol_Inf);

  // Freestream f
  f_Inf = (10.0/3.0+0.3)*epsi_Inf/max(kine_Inf,tke_min);


  /*--- Initialize the solution to the far-field state everywhere. ---*/
  CTurbKEVariable::SetConstants(constants);

  for (iPoint = 0; iPoint < nPoint; iPoint++)
    node[iPoint] = new CTurbKEVariable(kine_Inf, epsi_Inf, zeta_Inf, f_Inf,
                                       muT_Inf, Tm_Inf, Lm_Inf,
                                       nDim, nVar, config);

  if (runtime_averaging) {
    for (iPoint = 0; iPoint < nPoint; iPoint++) {
      average_node[iPoint] =
          new CTurbKEVariable(kine_Inf, epsi_Inf, zeta_Inf, f_Inf, muT_Inf,
                              Tm_Inf, Lm_Inf, nDim, nVar, config);
    }
  }

  /*--- MPI solution ---*/

  //TODO fix order of comunication the periodic should be first otherwise you have wrong values on the halo cell after restart
  Set_MPI_Solution(geometry, config);
  Set_MPI_Solution(geometry, config);
  if (runtime_averaging) {
    Set_MPI_Average_Solution(geometry, config);
    Set_MPI_Average_Solution(geometry, config);
  }

  /*--- Initializate quantities for SlidingMesh Interface ---*/

  unsigned long iMarker;

  SlidingState       = new su2double*** [nMarker];
  SlidingStateNodes  = new int*         [nMarker];

  for (iMarker = 0; iMarker < nMarker; iMarker++){

    SlidingState[iMarker]      = NULL;
    SlidingStateNodes[iMarker] = NULL;

    if (config->GetMarker_All_KindBC(iMarker) == FLUID_INTERFACE){

      SlidingState[iMarker]       = new su2double**[geometry->GetnVertex(iMarker)];
      SlidingStateNodes[iMarker]  = new int        [geometry->GetnVertex(iMarker)];

      for (iPoint = 0; iPoint < geometry->GetnVertex(iMarker); iPoint++){
        SlidingState[iMarker][iPoint] = new su2double*[nPrimVar+1];

        SlidingStateNodes[iMarker][iPoint] = 0;
        for (iVar = 0; iVar < nPrimVar+1; iVar++)
          SlidingState[iMarker][iPoint][iVar] = NULL;
      }

    }
  }

  /*-- Allocation of inlets has to happen in derived classes (not CTurbSolver),
   * due to arbitrary number of turbulence variables ---*/

  Inlet_TurbVars = new su2double**[nMarker];
  for (unsigned long iMarker = 0; iMarker < nMarker; iMarker++) {
    if (config->GetMarker_All_KindBC(iMarker) == INLET_FLOW ||
        config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_INLET) {
      Inlet_TurbVars[iMarker] = new su2double*[nVertex[iMarker]];
      for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
        Inlet_TurbVars[iMarker][iVertex] = new su2double[nVar];
        Inlet_TurbVars[iMarker][iVertex][0] = kine_Inf;
        Inlet_TurbVars[iMarker][iVertex][1] = epsi_Inf;
        Inlet_TurbVars[iMarker][iVertex][2] = zeta_Inf;
        Inlet_TurbVars[iMarker][iVertex][3] = f_Inf;
      }
    } else {
      Inlet_TurbVars[iMarker] = NULL;
    }
  }

}

CTurbKESolver::~CTurbKESolver(void) {

  if (constants != NULL) delete [] constants;

  unsigned long iMarker, iVertex;
  unsigned short iVar;

  if ( SlidingState != NULL ) {
    for (iMarker = 0; iMarker < nMarker; iMarker++) {
      if ( SlidingState[iMarker] != NULL ) {
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++)
          if ( SlidingState[iMarker][iVertex] != NULL ){
            for (iVar = 0; iVar < nPrimVar+1; iVar++)
              delete [] SlidingState[iMarker][iVertex][iVar];
            delete [] SlidingState[iMarker][iVertex];
          }
        delete [] SlidingState[iMarker];
      }
    }
    delete [] SlidingState;
  }
  
  if ( SlidingStateNodes != NULL ){
    for (iMarker = 0; iMarker < nMarker; iMarker++){
        if (SlidingStateNodes[iMarker] != NULL)
            delete [] SlidingStateNodes[iMarker];  
    }
    delete [] SlidingStateNodes;
  }

}

void CTurbKESolver::Preprocessing(CGeometry *geometry,
       CSolver **solver_container, CConfig *config, unsigned short iMesh,
       unsigned short iRKStep, unsigned short RunTime_EqSystem, bool Output) {

  unsigned long iPoint;

  unsigned long ExtIter = config->GetExtIter();
  bool disc_adjoint     = config->GetDiscrete_Adjoint();
  bool limiter_flow     = ((config->GetKind_SlopeLimit_Flow() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));
  bool limiter_turb     = ((config->GetKind_SlopeLimit_Turb() != NO_LIMITER) && (ExtIter <= config->GetLimiterIter()) && !(disc_adjoint && config->GetFrozen_Limiter_Disc()));


  for (iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Initialize the residual vector ---*/
    LinSysRes.SetBlock_Zero(iPoint);

  }

  /*--- Initialize the Jacobian matrices ---*/
  Jacobian.SetValZero();

  /*--- Upwind second order reconstruction ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config);
  }

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config);
  }

  /*--- Upwind second order reconstruction ---*/

  if (limiter_turb) SetSolution_Limiter(geometry, config);

  if (limiter_flow) solver_container[FLOW_SOL]->SetPrimitive_Limiter(geometry, config);

}

void CTurbKESolver::Postprocessing(CGeometry *geometry,
                                   CSolver **solver_container,
                                   CConfig *config, unsigned short iMesh) {

  su2double rho, nu, S;
  su2double tke, v2, zeta, muT;
  const bool model_split = (config->GetKind_HybridRANSLES() == MODEL_SPLIT);

  /*--- Compute mean flow and turbulence gradients ---*/
  if (config->GetKind_Gradient_Method() == GREEN_GAUSS) {
    SetSolution_Gradient_GG(geometry, config);
  }

  if (config->GetKind_Gradient_Method() == WEIGHTED_LEAST_SQUARES) {
    SetSolution_Gradient_LS(geometry, config);
  }

  su2double* VelInf = config->GetVelocity_FreeStreamND();
  su2double VelMag = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);
  const su2double L_inf = config->GetLength_Reynolds();

  CVariable** flow_node;
  if (model_split) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  const unsigned short realizability_limit = config->GetKind_v2f_Limit();
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint ++) {

    /*--- Compute turbulence scales ---*/

    rho = flow_node[iPoint]->GetDensity();
    nu  = flow_node[iPoint]->GetLaminarViscosity() / rho;
    S   = flow_node[iPoint]->GetStrainMag();

    /*--- Scalars ---*/

    tke = max(node[iPoint]->GetSolution(0), 0.0);
    v2  = max(node[iPoint]->GetSolution(2), 0.0);

    /*--- T & L ---*/

    node[iPoint]->SetTurbScales(nu, S, VelMag, L_inf, config);
    node[iPoint]->SetKolKineticEnergyRatio(nu);

    const su2double Tm = node[iPoint]->GetTurbTimescale();

    /*--- Compute the eddy viscosity ---*/

    muT = constants[0]*rho*v2*Tm;
    if (realizability_limit == EDDY_VISC_LIMIT) {
      const su2double C_lim = config->Getv2f_Realizability_Constant();
      muT = min(muT, C_lim*rho*tke/(sqrt(3)*max(S, 1E-16)));
    }

    node[iPoint]->SetmuT(muT);

    /* Compute resolution adequacy */
    if (model_split) {
      HybridMediator->ComputeResolutionAdequacy(geometry, solver_container,
                                                iPoint);
    }
  }
}

void CTurbKESolver::CalculateTurbScales(CSolver **solver_container,
                                        CConfig *config) {

  CVariable** flow_node;
  if (config->GetKind_HybridRANSLES() == MODEL_SPLIT) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  su2double* VelInf = config->GetVelocity_FreeStreamND();
  su2double VelMag = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);
  const su2double L_inf = config->GetLength_Reynolds();

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    const su2double nu  = flow_node[iPoint]->GetLaminarViscosity() /
                          flow_node[iPoint]->GetDensity();
    const su2double S   = flow_node[iPoint]->GetStrainMag();

    node[iPoint]->SetTurbScales(nu, S, VelMag, L_inf, config);
  }
}

void CTurbKESolver::Source_Residual(CGeometry *geometry,
        CSolver **solver_container, CNumerics *numerics,
        CNumerics *second_numerics, CConfig *config, unsigned short iMesh,
        unsigned short iRKStep) {

  const bool model_split = (config->GetKind_HybridRANSLES() == MODEL_SPLIT);

  CVariable** flow_node;
  if (model_split) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  /*--- Update turb scales ---*/
  // TODO: Evaluate if this should be left in.
  CalculateTurbScales(solver_container, config);

  for (unsigned long iPoint = 0; iPoint < nPointDomain; iPoint++) {

    numerics->SetTimeStep(
      solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time());

    numerics->SetCoord(geometry->node[iPoint]->GetCoord(), NULL);

    /*--- Conservative variables w/o reconstruction ---*/
    numerics->SetPrimitive(flow_node[iPoint]->GetPrimitive(), NULL);

    /*--- Gradient of the primitive and conservative variables ---*/
    numerics->SetPrimVarGradient(flow_node[iPoint]->GetGradient_Primitive(), NULL);

    /*--- Turbulent variables w/o reconstruction, and its gradient ---*/
    numerics->SetTurbVar(node[iPoint]->GetSolution(), NULL);
    numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), NULL);

    /*--- Set volume ---*/
    numerics->SetVolume(geometry->node[iPoint]->GetVolume());

    /*--- Pass T, L to the lengthscale ---*/
    numerics->SetTurbLengthscale(node[iPoint]->GetTurbLengthscale());
    numerics->SetTurbTimescale(node[iPoint]->GetTurbTimescale());

    /*--- Set vorticity and strain rate magnitude ---*/
    numerics->SetVorticity(flow_node[iPoint]->GetVorticity(), NULL);
    numerics->SetStrainMag(flow_node[iPoint]->GetStrainMag(), 0.0);

    /*--- Set the resolved Reynolds stress ---*/
    if (model_split) {
      HybridMediator->SetupRANSNumerics(solver_container, numerics, iPoint);
    }

    /*--- Compute the source term ---*/
    numerics->ComputeResidual(Residual, Jacobian_i, NULL, config);

    if (model_split) {
      CSourcePieceWise_TurbKE* v2f_numerics = dynamic_cast<CSourcePieceWise_TurbKE*>(numerics);
      const su2double sgs_production = v2f_numerics->GetSGSProduction();
      solver_container[FLOW_SOL]->average_node[iPoint]->SetSGSProduction(sgs_production);
      if (config->GetUse_Resolved_Turb_Stress()) {
        /*--- Production is not passed in, so output it instead ---*/
        const su2double production = numerics->GetProduction();
        solver_container[FLOW_SOL]->average_node[iPoint]->SetProduction(production);
      }
    }
    node[iPoint]->SetProduction(numerics->GetProduction());
    assert(node[iPoint]->GetProduction() == numerics->GetProduction());

    /*--- Subtract residual and the Jacobian ---*/
    LinSysRes.SubtractBlock(iPoint, Residual);
    Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

  }
}

void CTurbKESolver::Source_Template(CGeometry *geometry,
       CSolver **solver_container, CNumerics *numerics,
       CConfig *config, unsigned short iMesh) {
  // No-op
}

void CTurbKESolver::BC_HeatFlux_Wall(CGeometry *geometry,
       CSolver **solver_container, CNumerics *conv_numerics,
       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker,
       unsigned short iRKStep) {

  unsigned long iPoint, jPoint, iVertex, total_index;
  unsigned short iDim, iVar, jVar;
  su2double distance, wall_k;
  su2double density = 0.0, laminar_viscosity = 0.0;

  bool compressible = (config->GetKind_Regime() == COMPRESSIBLE);
  bool incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);
  CVariable** flow_node;
  if (config->GetKind_HybridRANSLES() == MODEL_SPLIT) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- distance to closest neighbor ---*/
      jPoint = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();
      distance = 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        const su2double dx = (geometry->node[iPoint]->GetCoord(iDim) -
                              geometry->node[jPoint]->GetCoord(iDim));
        distance += dx*dx;
      }
      distance = sqrt(distance);

      /*--- Set wall values ---*/
      if (compressible) {
        density = flow_node[iPoint]->GetDensity();
        laminar_viscosity = flow_node[iPoint]->GetLaminarViscosity();
      }
      if (incompressible) {
        density = flow_node[iPoint]->GetDensity();
        laminar_viscosity = flow_node[iPoint]->GetLaminarViscosity();
      }

      if (config->GetBoolUse_v2f_Explicit_WallBC()) {

        wall_k = node[jPoint]->GetSolution(0);

        // wall boundary conditions (https://turbmodels.larc.nasa.gov/k-e-zeta-f.html)
        Solution[0] = 0.0;
        Solution[1] = 2.0*laminar_viscosity*wall_k/(density*distance*distance);
        Solution[2] = 0.0;
        Solution[3] = 0.0;

        /*--- Set the solution values and zero the residual ---*/
        node[iPoint]->SetSolution_Old(Solution);
        node[iPoint]->SetSolution(Solution);
        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Change rows of the Jacobian (to be just 1 on the diagonal) ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }

      } else {
        /*--- The conditions below are a mess.  Here is why...

          Epsilon at the wall depends on k at the first node off the
          wall.  It is observed that, if this is enforced following the
          'standard' SU2 procedure (see e.g., the wall conditions for SA
          or SST), it is difficult to converge epsilon and that the
          large residuals tend to appear very close to the wall.

          Thus, it is preferred to enforce the BC by including the
          equation in the system being solved directly.  This is easy
          enough to do in the residual by setting

          Res[iPoint*nVal+1] = (epsilon - desired epsilon based on k).

          It is also simple to set the appropriate Jacobian.  However,
          Vol/dt is added to the diagonal of the Jacobian in
          CTurbSolver::ImplicitEuler_Iteration.  So, to ensure we get
          the right Jacobian, we subtract Vol/dt here.  But, obviously
          this is very brittle b/c it assumes we're using backward
          Euler. ---*/

        wall_k = node[jPoint]->GetSolution(0);

        const su2double Vol = geometry->node[iPoint]->GetVolume();
        const su2double dt = solver_container[FLOW_SOL]->node[iPoint]->GetDelta_Time();

        // Set k, v2, and f as for explicit method,
        // but keep epsilon as is b/c it is dealt with below
        Solution[0] = 0.0;
        Solution[1] = node[iPoint]->GetSolution(1);
        Solution[2] = 0.0;
        Solution[3] = 0.0;

        /*--- Set the solution values and zero the residual ---*/
        node[iPoint]->SetSolution_Old(Solution);
        node[iPoint]->SetSolution(Solution);
        LinSysRes.SetBlock_Zero(iPoint);

        /*--- Change rows of the Jacobian (to be just 1 on the diagonal) ---*/
        for (iVar = 0; iVar < nVar; iVar++) {
          total_index = iPoint*nVar+iVar;
          Jacobian.DeleteValsRowi(total_index);
        }


        for (iVar=0; iVar<nVar; iVar++) {
          for (jVar=0; jVar<nVar; jVar++) {
            Jacobian_j[iVar][jVar] = 0.0;
          }
        }

        const su2double rho_epsi_wall = 2.0*laminar_viscosity*wall_k/(distance*distance);
        const su2double rho_epsi_res = density*node[iPoint]->GetSolution(1) - rho_epsi_wall;
        LinSysRes.SetBlock(iPoint, 1, rho_epsi_res);
        Jacobian.DeleteValsRowi(iPoint*nVar+1); // zeros this row and puts 1 on diagonal

        // WARNING: Hackery
        // Subtract Vol/dt from jacobian to offset addition of this later on
        Jacobian_j[1][1] = Vol/dt;
        Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_j);

        // Jacobian wrt conserved state at node j
        Jacobian_j[1][0] = -2.0*laminar_viscosity/(density*distance*distance);
        Jacobian_j[1][1] = 0.0;
        Jacobian_j[1][2] = 0.0;
        Jacobian_j[1][3] = 0.0;

        Jacobian.AddBlock(iPoint, jPoint, Jacobian_j);
      }

    }
  }
}

void CTurbKESolver::BC_Isothermal_Wall(CGeometry *geometry,
       CSolver **solver_container, CNumerics *conv_numerics,
       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker,
       unsigned short iRKStep) {

  // Turbulence model BC doesn't care whether wall is isothermal or
  // heat flux condition, so just call the heat flux BC function to
  // minimize code duplication
  BC_HeatFlux_Wall(geometry, solver_container, conv_numerics,
                   visc_numerics, config, val_marker, iRKStep);

}

void CTurbKESolver::BC_Far_Field(CGeometry *geometry,
       CSolver **solver_container, CNumerics *conv_numerics,
       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker,
       unsigned short iRKStep) {

  unsigned long iPoint, iVertex;
  su2double *Normal, *V_infty, *V_domain;
  unsigned short iVar, iDim;

  bool grid_movement = config->GetGrid_Movement();

  CVariable** flow_node;
  if (config->GetKind_HybridRANSLES() == MODEL_SPLIT) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  Normal = new su2double[nDim];

  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e, not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Allocate the value at the infinity ---*/
      V_infty = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = flow_node[iPoint]->GetPrimitive();
      conv_numerics->SetPrimitive(V_domain, V_infty);

      /*--- Set turbulent variable at infinity ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      Solution_j[0] = node[iPoint]->GetSolution(0);
      Solution_j[1] = node[iPoint]->GetSolution(1);
      Solution_j[2] = node[iPoint]->GetSolution(2);
      Solution_j[3] = node[iPoint]->GetSolution(3);

      conv_numerics->SetTurbVar(Solution_i, Solution_j);
      /*--- Set Normal (it is necessary to change the sign) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      /*--- Grid Movement ---*/
      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute residuals and Jacobians ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);

      /*--- Add residuals and Jacobians ---*/
      LinSysRes.AddBlock(iPoint, Residual);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

      // FIXME: Jacobian should also depend on Jacobian_j b/c
      // Solution_j = Solution_i
    }
  }

  delete [] Normal;
}

void CTurbKESolver::BC_Inlet(CGeometry *geometry, CSolver **solver_container,
                             CNumerics *conv_numerics, CNumerics *visc_numerics,
                             CConfig *config, unsigned short val_marker,
                             unsigned short iRKStep) {

  unsigned short iVar, iDim;
  unsigned long iVertex, iPoint, Point_Normal;
  su2double *V_inlet, *V_domain, *Normal;

  /*--- Check that the inlet has been set up correctly ---*/
  assert(Inlet_TurbVars[val_marker] != NULL);

  Normal = new su2double[nDim];

  bool grid_movement  = config->GetGrid_Movement();

  CVariable** flow_node;
  if (config->GetKind_HybridRANSLES() == MODEL_SPLIT) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  string Marker_Tag = config->GetMarker_All_TagBound(val_marker);

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {

    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Normal vector for this vertex (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++) Normal[iDim] = -Normal[iDim];

      /*--- Allocate the value at the inlet ---*/
      V_inlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = flow_node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_inlet);

      /*--- Set the turbulent variable states. Use free-stream KE
       values for the turbulent state at the inflow. ---*/
      for (iVar = 0; iVar < nVar; iVar++)
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);

      Solution_j[0] = Inlet_TurbVars[val_marker][iVertex][0];
      Solution_j[1] = Inlet_TurbVars[val_marker][iVertex][1];
      Solution_j[2] = Inlet_TurbVars[val_marker][iVertex][2];
      Solution_j[3] = Inlet_TurbVars[val_marker][iVertex][3];

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set various other quantities in the solver class ---*/
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);

//      /*--- Viscous contribution, commented out because serious convergence problems  ---*/
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
//                              geometry->node[Point_Normal]->GetCoord());
//
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//      visc_numerics->SetPrimitive(V_domain, V_inlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(), node[iPoint]->GetGradient());
//
//      /*--- Compute residual, and Jacobians ---*/
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);

      // Set f residual correctly to get right update
      LinSysRes.SetBlock(iPoint, 3, Solution_i[3] - f_Inf);
      Jacobian.DeleteValsRowi(iPoint*nVar+3);

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

void CTurbKESolver::BC_Outlet(CGeometry *geometry,
       CSolver **solver_container, CNumerics *conv_numerics,
       CNumerics *visc_numerics, CConfig *config, unsigned short val_marker,
       unsigned short iRKStep) {

  unsigned long iPoint, iVertex, Point_Normal;
  unsigned short iVar, iDim;
  su2double *V_outlet, *V_domain, *Normal;

  bool grid_movement  = config->GetGrid_Movement();

  CVariable** flow_node;
  if (config->GetKind_HybridRANSLES() == MODEL_SPLIT) {
    /*--- Use explicit average values instead of fluctuating values ---*/
    flow_node = solver_container[FLOW_SOL]->average_node;
  } else {
    flow_node = solver_container[FLOW_SOL]->node;
  }

  Normal = new su2double[nDim];

  /*--- Loop over all the vertices on this boundary marker ---*/
  for (iVertex = 0; iVertex < geometry->nVertex[val_marker]; iVertex++) {
    iPoint = geometry->vertex[val_marker][iVertex]->GetNode();

    /*--- Check if the node belongs to the domain (i.e., not a halo node) ---*/
    if (geometry->node[iPoint]->GetDomain()) {

      /*--- Index of the closest interior node ---*/
      Point_Normal = geometry->vertex[val_marker][iVertex]->GetNormal_Neighbor();

      /*--- Allocate the value at the outlet ---*/
      V_outlet = solver_container[FLOW_SOL]->GetCharacPrimVar(val_marker, iVertex);

      /*--- Retrieve solution at the farfield boundary node ---*/
      V_domain = flow_node[iPoint]->GetPrimitive();

      /*--- Set various quantities in the solver class ---*/
      conv_numerics->SetPrimitive(V_domain, V_outlet);

      /*--- Set the turbulent variables. Here we use a Neumann BC such
       that the turbulent variable is copied from the interior of the
       domain to the outlet before computing the residual.
       Solution_i --> TurbVar_internal,
       Solution_j --> TurbVar_outlet ---*/
      for (iVar = 0; iVar < nVar; iVar++) {
        Solution_i[iVar] = node[iPoint]->GetSolution(iVar);
        Solution_j[iVar] = node[iPoint]->GetSolution(iVar);
      }

      conv_numerics->SetTurbVar(Solution_i, Solution_j);

      /*--- Set Normal (negate for outward convention) ---*/
      geometry->vertex[val_marker][iVertex]->GetNormal(Normal);
      for (iDim = 0; iDim < nDim; iDim++)
      Normal[iDim] = -Normal[iDim];
      conv_numerics->SetNormal(Normal);

      if (grid_movement)
        conv_numerics->SetGridVel(geometry->node[iPoint]->GetGridVel(),
                                  geometry->node[iPoint]->GetGridVel());

      /*--- Compute the residual using an upwind scheme ---*/
      conv_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
      LinSysRes.AddBlock(iPoint, Residual);

      /*--- Jacobian contribution for implicit integration ---*/
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_i);
      Jacobian.AddBlock(iPoint, iPoint, Jacobian_j); // since soln_j = soln_i

//      /*--- Viscous contribution, commented out because serious convergence problems  ---*/
//      visc_numerics->SetCoord(geometry->node[iPoint]->GetCoord(),
//                              geometry->node[Point_Normal]->GetCoord());
//      visc_numerics->SetNormal(Normal);
//
//      /*--- Conservative variables w/o reconstruction ---*/
//      visc_numerics->SetPrimitive(V_domain, V_outlet);
//
//      /*--- Turbulent variables w/o reconstruction, and its gradients ---*/
//      visc_numerics->SetTurbVar(Solution_i, Solution_j);
//      visc_numerics->SetTurbVarGradient(node[iPoint]->GetGradient(),
//                                        node[iPoint]->GetGradient());
//
//      /*--- Compute residual, and Jacobians ---*/
//      visc_numerics->ComputeResidual(Residual, Jacobian_i, Jacobian_j, config);
//
//      /*--- Subtract residual, and update Jacobians ---*/
//      LinSysRes.SubtractBlock(iPoint, Residual);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_i);
//      Jacobian.SubtractBlock(iPoint, iPoint, Jacobian_j); // since soln_j = soln_i

    }
  }

  /*--- Free locally allocated memory ---*/
  delete[] Normal;

}

su2double* CTurbKESolver::GetConstants() {
  return constants;
}

void CTurbKESolver::SetInletAtVertex(su2double *val_inlet,
                                     unsigned short iMarker,
                                     unsigned long iVertex,
                                     CConfig* config) {
  for (unsigned short iVar = 0; iVar < nVar; iVar++) {
    Inlet_TurbVars[iMarker][iVertex][iVar] = val_inlet[nDim+2+nDim + iVar];
  }
}

su2double CTurbKESolver::GetInletAtVertex(su2double *val_inlet,
                                           unsigned long val_inlet_point,
                                           unsigned short val_kind_marker,
                                           string val_marker,
                                           CGeometry *geometry,
                                           CConfig *config) {

  /*--- Local variables ---*/

  unsigned short iMarker, iDim;
  unsigned long iPoint, iVertex;
  su2double Area = 0.0;
  su2double Normal[3] = {0.0,0.0,0.0};

  /*--- Alias positions within inlet file for readability ---*/

  if (val_kind_marker == INLET_FLOW || val_kind_marker == SUPERSONIC_INLET) {

    // Offset for the coordinates + flow variables
    const unsigned short offset = nDim + 2 + nDim;

    for (iMarker = 0; iMarker < config->GetnMarker_All(); iMarker++) {
      if ((config->GetMarker_All_KindBC(iMarker) == INLET_FLOW ||
           config->GetMarker_All_KindBC(iMarker) == SUPERSONIC_INLET) &&
          (config->GetMarker_All_TagBound(iMarker) == val_marker)) {
        
        for (iVertex = 0; iVertex < nVertex[iMarker]; iVertex++){

          iPoint = geometry->vertex[iMarker][iVertex]->GetNode();

          if (iPoint == val_inlet_point) {

            /*-- Compute boundary face area for this vertex. ---*/

            geometry->vertex[iMarker][iVertex]->GetNormal(Normal);
            Area = 0.0;
            for (iDim = 0; iDim < nDim; iDim++)
              Area += Normal[iDim]*Normal[iDim];
            Area = sqrt(Area);

            /*--- Access and store the inlet variables for this vertex. ---*/

            for (unsigned short iVar = 0; iVar < nVar; iVar++) {
              val_inlet[offset + iVar] = Inlet_TurbVars[iMarker][iVertex][iVar];
            }

            /*--- Exit once we find the point. ---*/

            return Area;

          }
        }
      }
    }

  }

  /*--- If we don't find a match, then the child point is not on the
   current inlet boundary marker. Return zero area so this point does
   not contribute to the restriction operator and continue. ---*/
  
  return Area;
  
}

void CTurbKESolver::SetUniformInlet(CConfig* config, unsigned short iMarker) {

  assert(nVertex != NULL);
  for(unsigned long iVertex=0; iVertex < nVertex[iMarker]; iVertex++){
    Inlet_TurbVars[iMarker][iVertex][0] = kine_Inf;
    Inlet_TurbVars[iMarker][iVertex][1] = epsi_Inf;
    Inlet_TurbVars[iMarker][iVertex][2] = zeta_Inf;
    Inlet_TurbVars[iMarker][iVertex][3] = f_Inf;
  }
}

void CTurbKESolver::BC_Supersonic_Inlet(CGeometry *geometry,
                                       CSolver **solver_container,
                                       CNumerics *conv_numerics,
                                       CNumerics *visc_numerics,
                                       CConfig *config,
                                       unsigned short val_marker,
                                       unsigned short iRKStep) {
  BC_Inlet(geometry, solver_container, conv_numerics, visc_numerics,
           config, val_marker, iRKStep);
}

void CTurbKESolver::BC_Supersonic_Outlet(CGeometry *geometry,
                                        CSolver **solver_container,
                                        CNumerics *conv_numerics,
                                        CNumerics *visc_numerics,
                                        CConfig *config,
                                        unsigned short val_marker,
                                        unsigned short iRKStep) {
  BC_Outlet(geometry, solver_container, conv_numerics, visc_numerics,
            config, val_marker, iRKStep);
}

/*!
 * \file hybrid_RANS_LES_forcing.cpp
 * \brief Implementation of the hybrid RANS/LES forcing terms
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

#include "../include/hybrid_RANS_LES_forcing.hpp"
#include "../include/solver_structure.hpp"

CHybridForcingAbstractBase::CHybridForcingAbstractBase(
                              const unsigned short nDim,
                              const unsigned long nPoint,
                              const unsigned long nPointDomain)
  : nDim(nDim), nVar(nDim), nVarGrad(nDim),
    nPoint(nPoint), nPointDomain(nPointDomain) {
}

CHybridForcingAbstractBase::~CHybridForcingAbstractBase() {
}

CHybridForcingTG0::CHybridForcingTG0(const unsigned short nDim,
                               const unsigned long nPoint,
                               const unsigned long nPointDomain)
    : CHybridForcingAbstractBase(nDim, nPoint, nPointDomain),
      forcing_scale(4), C_F(8) {

  node = new su2double*[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new su2double[nVar];
  }

}


CHybridForcingTG0::CHybridForcingTG0(CGeometry* geometry, const CConfig* config)
    : CHybridForcingAbstractBase(geometry->GetnDim(),
                                 geometry->GetnPoint(),
                                 geometry->GetnPointDomain()),
      forcing_scale(config->GetHybrid_Forcing_Vortex_Length()),
      C_F(config->GetHybrid_Forcing_Strength()) {

  node = new su2double*[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new su2double[nVar];
  }

}

CHybridForcingTG0::~CHybridForcingTG0() {

  if (node != NULL) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      delete [] node[iPoint];
    }
    delete [] node;
  }
}

void CHybridForcingTG0::ComputeForcingField(CSolver** solver, CGeometry *geometry,
                                            CConfig *config) {

  const su2double TKE_MIN = 1.0E-8;
  const su2double V2_MIN = 2.0/3*TKE_MIN;
  const su2double TDR_MIN = 1.0E-8;

  const unsigned short kind_time_marching = config->GetUnsteady_Simulation();

  assert(kind_time_marching == TIME_STEPPING   ||
         kind_time_marching == DT_STEPPING_1ST ||
         kind_time_marching == DT_STEPPING_2ND );

  const su2double time = config->GetCurrent_UnstTimeND();
  const su2double dt = config->GetDelta_UnstTimeND();
  assert(time >= 0);
  assert(dt > 0);

  /*--- Allocate some scratch arrays to avoid continual reallocation ---*/
  su2double h[nDim]; // Initial TG vortex field.
  su2double Lsgs; // SGS length scale
  su2double x[nDim]; // Position
  su2double Lmesh[nDim]; // Mesh length scales in coord direction (computed from res tensor)
  su2double dwall; // distance to the wall
  su2double uf[nDim];

  // Domain lengths for periodic directions
  su2double *D = config->GetHybrid_Forcing_Periodic_Length();

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- This velocity retrieval works only for compressible ---*/
    assert(config->GetKind_Regime() == COMPRESSIBLE);
    const su2double* prim_vars = solver[FLOW_SOL]->node[iPoint]->GetPrimitive();
    const su2double* prim_mean =
      solver[FLOW_SOL]->average_node[iPoint]->GetPrimitive();

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      /*--- Copy the values (and not the pointer), since we're changing them. ---*/
      x[iDim] = geometry->node[iPoint]->GetCoord(iDim);
      x[iDim] = TransformCoords(x[iDim], prim_mean[iDim+1], time);
    }

    /*--- Setup the TG vortex ---*/

    // total k, epsilon, and v2 from model
    // TODO: Generalize beyond v2-f
    if (config->GetKind_Turb_Model() != KE) {
      SU2_MPI::Error("Hybrid forcing is only compatible with the v2-f model.", CURRENT_FUNCTION);
    }

    const su2double ktot = max(solver[TURB_SOL]->node[iPoint]->GetSolution(0), TKE_MIN);
    const su2double tdr  = max(solver[TURB_SOL]->node[iPoint]->GetSolution(1), TDR_MIN);
    const su2double v2 = max(solver[TURB_SOL]->node[iPoint]->GetSolution(2), V2_MIN);

    // ratio of modeled to total TKE
    const su2double alpha = solver[FLOW_SOL]->average_node[iPoint]->GetKineticEnergyRatio();

    // average of r_M
    const su2double resolution_adequacy =
      solver[FLOW_SOL]->average_node[iPoint]->GetResolutionAdequacy();

    const su2double density = solver[FLOW_SOL]->average_node[iPoint]->GetSolution(0);
    const su2double nu =
      solver[FLOW_SOL]->average_node[iPoint]->GetLaminarViscosity() / density;

    // TODO: Move length, timescale calculations to CTurbKEVariable class,
    // and keep track of all three timescales, plus alpha_kol
    const su2double V2F_CETA = 70;
    Lsgs = forcing_scale * pow(alpha * ktot, 1.5) / tdr;
    const su2double Lkol = V2F_CETA*pow(nu, 0.75)/pow(tdr, 0.25);
    Lsgs = max(Lsgs, Lkol);

    // FIXME: I think this is equivalent to repo version of CDP,but
    // not consistent with paper description, except for orthogonal
    // grids aligned with coordinate axes.  Check with Sigfried.
    const su2double* const* ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
    Lmesh[0] = ResolutionTensor[0][0];
    Lmesh[1] = ResolutionTensor[1][1];
    Lmesh[2] = ResolutionTensor[2][2];

    // Get dwall
    dwall = geometry->node[iPoint]->GetWall_Distance();

    // Compute TG velocity at this point
    this->SetTGField(x, Lsgs, Lmesh, D, dwall, h);

    const su2double V2F_CT = 6.0;
    const su2double T1 = alpha*ktot/tdr;
    const su2double Tkol = V2F_CT*sqrt(nu/tdr);
    const su2double Tsgs = max(T1, Tkol);

    const su2double Ftar = this->GetTargetProduction(v2, Tsgs, alpha);

    // Compute PFtest
    su2double PFtest = 0.0;
    for (unsigned short iDim=0; iDim<nDim; iDim++) {
      uf[iDim] = prim_vars[iDim+1] - prim_mean[iDim+1];
      PFtest += uf[iDim]*h[iDim];
    }
    PFtest *= Ftar*dt;

    const su2double Cnu = 1.0;
    const su2double alpha_kol = min(Cnu*std::sqrt(nu*tdr)/ktot, 1.0);

    const su2double eta = this->ComputeScalingFactor(Ftar, resolution_adequacy,
                                                     alpha, alpha_kol, PFtest);

    /*--- Store eta*h so we can compute the derivatives ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      node[iPoint][iDim] = eta*h[iDim];
    }

  } // end loop over points
}

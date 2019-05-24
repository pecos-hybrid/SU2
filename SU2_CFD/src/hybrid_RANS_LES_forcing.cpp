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
#include "../include/variable_structure_v2f.hpp"

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

  const unsigned short kind_time_marching = config->GetUnsteady_Simulation();

  assert(kind_time_marching == TIME_STEPPING   ||
         kind_time_marching == DT_STEPPING_1ST ||
         kind_time_marching == DT_STEPPING_2ND );

  const su2double time = config->GetCurrent_UnstTimeND();
  assert(time >= 0);

  /*--- Allocate some scratch arrays to avoid continual reallocation ---*/
  su2double h[nDim]; // Initial TG vortex field.
  su2double Lsgs; // SGS length scale
  su2double Tsgs; // SGS time scale
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
    const su2double v2 = max(solver[TURB_SOL]->node[iPoint]->GetSolution(2), V2_MIN);

    // ratio of modeled to total TKE
    su2double alpha = solver[FLOW_SOL]->average_node[iPoint]->GetKineticEnergyRatio();
    alpha = max(alpha, 1e-8);

    // average of r_M
    const su2double resolution_adequacy =
      solver[FLOW_SOL]->average_node[iPoint]->GetResolutionAdequacy();

    const su2double L_typical = solver[TURB_SOL]->node[iPoint]->GetTypicalLengthscale();
    const su2double L_kol = solver[TURB_SOL]->node[iPoint]->GetKolLengthscale();
    const su2double T_typical = solver[TURB_SOL]->node[iPoint]->GetTypicalTimescale();
    const su2double T_kol = solver[TURB_SOL]->node[iPoint]->GetKolTimescale();
    assert(L_typical >= 0);
    assert(L_kol > 0);
    assert(T_typical >= 0);
    assert(T_kol > 0);

    Lsgs = max(forcing_scale * pow(alpha, 1.5) * L_typical, L_kol);
    Lsgs = max(Lsgs, L_kol);
    Tsgs = alpha * T_typical;
    Tsgs = max(Tsgs, T_kol);

    // FIXME: I think this is equivalent to repo version of CDP,but
    // not consistent with paper description, except for orthogonal
    // grids aligned with coordinate axes.  Check with Sigfried.
    const su2double* const* ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
    Lmesh[0] = ResolutionTensor[0][0];
    Lmesh[1] = ResolutionTensor[1][1];
    Lmesh[2] = ResolutionTensor[2][2];

    // Get dwall
    dwall = geometry->node[iPoint]->GetWall_Distance();

    if (config->isHybrid_Forced_Axi()) {
      // Angular periodic version 
      this->SetAxiTGField(x, Lsgs, Lmesh, D, dwall, h);
    } else {
      // Compute TG velocity at this point
      this->SetTGField(x, Lsgs, Lmesh, D, dwall, h);
    }

    const su2double Ftar = this->GetTargetProduction(v2, Tsgs, alpha);

    // Compute PFtest
    su2double PFtest = 0.0;
    for (unsigned short iDim=0; iDim<nDim; iDim++) {
      uf[iDim] = prim_vars[iDim+1] - prim_mean[iDim+1];
      PFtest += uf[iDim]*h[iDim];
    }
    PFtest *= Ftar;

    const su2double alpha_kol = solver[TURB_SOL]->node[iPoint]->GetKolKineticEnergyRatio();

    const su2double eta = this->ComputeScalingFactor(Ftar, resolution_adequacy,
                                                     alpha, alpha_kol, PFtest);

    /*--- Store eta*h so we can compute the derivatives ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      node[iPoint][iDim] = eta*h[iDim];
    }

  } // end loop over points
}

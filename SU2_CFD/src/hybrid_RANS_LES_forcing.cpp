/*!
 * \file hybrid_RANS_LES_forcing.cpp
 * \brief Implementation of the hybrid RANS/LES forcing terms
 * \author C. Pederson
 * \version 6.2.0 "Falcon"
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

    // Get dwall
    dwall = geometry->node[iPoint]->GetWall_Distance();

    if (config->isHybrid_Forced_Axi()) {
      // Angular periodic version 
      this->SetAxiTGField(x, Lsgs, D, dwall, h);
    } else {
      // Compute TG velocity at this point
      this->SetTGField(x, Lsgs, D, dwall, h);
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

    /*---
     * Ensure that forcing doesn't push the flow into an unphysical state.
     *
     * Forcing is supposed to just exchange energy between modeled
     * turbulence and the resolved flow.  So we shouldn't have to worry about
     * the internal energy (or temperature) here.  If the resolved flow gains
     * more energy, alpha*k_tot should go down to match.
     *
     * But there's two conditions where a mismatch can occur:
     *   1. Since k_tot is an average quantity and the resolved
     *      kinetic energy is a fluctuating quantity, a local hotspot in
     *      resolved kinetic energy can exist, independent of k_tot.
     *   2. We expect alpha to lag behind the actual value as the amount of
     *      resolved turbulence increases.  So the resolved state will
     *      instantly feel the effects of forcing, but alpha will not
     *      adjust automatically.  This can lead to a double-accounting
     *      of energy and therefore a negative internal energy
     *      (or temperature).
     *
     * This correction clips the forcing when such a state may occur.
     * It is an estimate, since we don't have the actual values of
     * d(rho e)/dt.
     * ---*/
    su2double clip = 1.0;
    if (eta > 0) {
      const su2double k_tot = max(solver[TURB_SOL]->node[iPoint]->GetSolution(0), su2double(0.0));
      const su2double total_energy = solver[FLOW_SOL]->node[iPoint]->GetEnergy();
      const su2double resolved_ke = 0.5*solver[FLOW_SOL]->node[iPoint]->GetVelocity2();
      const su2double static_energy = total_energy - resolved_ke - alpha*k_tot;
      su2double resolved_u[3];
      for (unsigned short iDim=0; iDim<nDim; iDim++) {
        resolved_u[iDim] = prim_vars[iDim+1];
      }
      if (static_energy < 0) {
        ostringstream error_msg;
        error_msg << "Detected a negative static energy!\n";
        error_msg << "  total energy: " << total_energy << "\n";
        error_msg << "  resolved KE:  " << resolved_ke << "\n";
        error_msg << "  alpha:        " << alpha << "\n";
        error_msg << "  total TKE:    " << k_tot << "\n";
        SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
      }


      // By assumption, we are using 'TIME_STEPPING' for these
      // simulations, such that dt is the same for each node and we can
      // just grab it off the first node.
      assert(config->GetUnsteady_Simulation() == TIME_STEPPING);
      const su2double delta_t = solver[FLOW_SOL]->node[0]->GetDelta_Time();

      /*--- Calculate proposed forcing vector ---*/
      su2double current_f[3]; // instantaneous value
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        current_f[iDim] = eta*h[iDim];
      }
      const su2double* avg_f = solver[FLOW_SOL]->average_node[iPoint]->GetForcingVector();

      clip = CheckRealizability(current_f, avg_f, resolved_u, delta_t, static_energy);
    }

    /*--- Store eta*h so we can compute the derivatives ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      node[iPoint][iDim] = clip*eta*h[iDim];
    }

  } // end loop over points
}

su2double CHybridForcingTG0::ComputeScalingFactor(
                     const su2double Ftar,
                     const su2double resolution_adequacy,
                     const su2double alpha,
                     const su2double alpha_kol,
                     const su2double PFtest) const {

  su2double eta = 0.0;

  // TODO: Compare this with Sigfried's improved version once channel
  // validation is successful.
  if ( (PFtest >= 0.0) && (resolution_adequacy < 1.0) ) {
    const su2double Sr = tanh(1.0 - 1.0/sqrt(resolution_adequacy));
    if (alpha <= alpha_kol) {
      eta = -Ftar * Sr * (alpha - alpha_kol);
    } else {
      eta = -Ftar * Sr;
    }
  }

  return eta;
}

su2double CHybridForcingTG0::CheckRealizability(const su2double* current_f,
                                                const su2double* avg_f,
                                                const su2double* resolved_u,
                                                const su2double delta_t,
                                                const su2double static_energy) const {

  // This is the added *specific* energy.  We omit density from
  // both the static energy and the added resolved kinetic energy.
  su2double added_energy = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    added_energy += delta_t*resolved_u[iDim]*(current_f[iDim] - avg_f[iDim]);
  }

  const su2double small_correction = 0.1;
  return max(0.0, min(1.0, static_energy / added_energy));
}

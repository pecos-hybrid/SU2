/*!
 * \file hybrid_RANS_LES_model.cpp
 * \brief Describes the hybrid RANS/LES models
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
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/hybrid_RANS_LES_model.hpp"
#include "../include/numerics_structure.hpp"
#include "../include/solver_structure.hpp"
#ifdef HAVE_LAPACK
#include "mkl.h"
#include "mkl_lapacke.h"
#endif


CHybrid_Visc_Anisotropy::CHybrid_Visc_Anisotropy(unsigned short nDim)
    : nDim(nDim) {
  eddy_visc_anisotropy = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    eddy_visc_anisotropy[iDim] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++)
      eddy_visc_anisotropy[iDim][jDim] = 0.0;
  }
}

CHybrid_Visc_Anisotropy::~CHybrid_Visc_Anisotropy() {
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] eddy_visc_anisotropy[iDim];
  delete [] eddy_visc_anisotropy;
}

su2double** CHybrid_Visc_Anisotropy::GetViscAnisotropy() {
  return eddy_visc_anisotropy;
}

CHybrid_Isotropic_Visc::CHybrid_Isotropic_Visc(unsigned short nDim)
: CHybrid_Visc_Anisotropy(nDim) {
  for (unsigned short iDim=0; iDim < nDim; iDim++) {
    for (unsigned short jDim=0; jDim < nDim; jDim++) {
      eddy_visc_anisotropy[iDim][jDim] = (su2double)(iDim == jDim);
    }
  }
}

void CHybrid_Isotropic_Visc::CalculateViscAnisotropy() {
}

CHybrid_Aniso_Q::CHybrid_Aniso_Q(unsigned short nDim)
  : CHybrid_Visc_Anisotropy(nDim), rans_weight(1.0) {
}

void CHybrid_Aniso_Q::CalculateViscAnisotropy() {
#ifndef NDEBUG
  if (rans_weight < 0 || rans_weight > 1) {
    ostringstream error_msg;
    error_msg << "Isotropic weighting in hybrid RANS/LES was not in range [0,1]" << endl;
    error_msg << "    weight = " << rans_weight << endl;
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }
#endif
  su2double LES_weight = (1.0 - rans_weight);

  unsigned short iDim, jDim;

  // Qstar is already normalized.
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eddy_visc_anisotropy[iDim][jDim] = rans_weight*double(iDim == jDim);
      eddy_visc_anisotropy[iDim][jDim] += LES_weight*sqrt(nDim) *
                                          Qstar[iDim][jDim];
    }
  }
}


inline void CHybrid_Aniso_Q::SetTensor(su2double** val_approx_struct_func) {
  Qstar = val_approx_struct_func;
}

inline void CHybrid_Aniso_Q::SetScalars(vector<su2double> val_scalars) {
#ifndef NDEBUG
  if (val_scalars.size() != 1) {
    ostringstream error_msg;
    error_msg << "Improper number of scalars passed to anisotropy model!" << endl;
    error_msg << "    Expected: 1" << endl;
    error_msg << "    Found:    " << val_scalars.size() << endl;
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }
#endif
  rans_weight = val_scalars[0];
}



CHybrid_Mediator::CHybrid_Mediator(unsigned short nDim, CConfig* config,
                                   string filename)
   : nDim(nDim), C_sf(0.367), config(config) {

  hybrid_forcing = new CHybridForcing(nDim);

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  Q = NULL;
  Qapprox = NULL;
  invLengthTensor = NULL;

  if (config->GetKind_Hybrid_Resolution_Indicator() == RK_INDICATOR) {
    // TODO: Evaluate if this is all necessary

    /*--- Allocate the approximate structure function (used in calcs) ---*/

    Q = new su2double*[nDim];
    Qapprox = new su2double*[nDim];
    invLengthTensor = new su2double*[nDim];
    for (unsigned int iDim = 0; iDim < nDim; iDim++) {
      Q[iDim] = new su2double[nDim];
      Qapprox[iDim] = new su2double[nDim];
      invLengthTensor[iDim] = new su2double[nDim];
    }

    /*--- Load the constants for mapping M to M-tilde ---*/

    if (filename == "") {
      if (rank == MASTER_NODE) {
        cout << "WARNING: No file given for hybrid RANS/LES constants." << endl;
        cout << "         Default (hardcoded) values used." << endl;
      }
      constants.resize(3, vector<su2double>(15));

      constants[0][0] =  2.4253168624203;
      constants[0][1] =  0.3377273122202;
      constants[0][2] =  0.2454150824949;
      constants[0][3] = -0.4015732570841;
      constants[0][4] = -0.3023468205458;
      constants[0][5] = -0.1386196773678;
      constants[0][6] =  0.0451966752099;
      constants[0][7] =  0.0325650620151;
      constants[0][8] =  0.0093904940685;
      constants[0][9] = -0.0052447144608;
      constants[0][10] = -0.0019607919411;
      constants[0][11] = -0.0005522138218;
      constants[0][12] =  0.0013947282467;
      constants[0][13] =  0.0012723863199;
      constants[0][14] = -0.0000420559137;

      constants[1][0] =  0.6999425502058;
      constants[1][1] =  0.3056790968854;
      constants[1][2] = -0.1914576501370;
      constants[1][3] =  0.0713376305722;
      constants[1][4] =  0.2874057660774;
      constants[1][5] =  0.1107104307784;
      constants[1][6] = -0.0215754933753;
      constants[1][7] = -0.0652953391552;
      constants[1][8] = -0.0460413983614;
      constants[1][9] = -0.0131511446213;
      constants[1][10] =  0.0015258919631;
      constants[1][11] =  0.0046851430319;
      constants[1][12] =  0.0046149483796;
      constants[1][13] =  0.0020781858721;
      constants[1][14] =  0.0001722924891;

      constants[2][0] = -0.1451211648913;
      constants[2][1] = -0.0419089159238;
      constants[2][2] = -0.0090912831194;
      constants[2][3] =  0.0120968852318;
      constants[2][4] = -0.0318033690621;
      constants[2][5] = -0.0157539031345;
      constants[2][6] = -0.0007323909092;
      constants[2][7] =  0.0105452780759;
      constants[2][8] =  0.0089366657596;
      constants[2][9] =  0.0030581437094;
      constants[2][10] = -0.0000170956796;
      constants[2][11] = -0.0009297436006;
      constants[2][12] = -0.0010752469431;
      constants[2][13] = -0.0005650127892;
      constants[2][14] = -0.0000591358738;

    } else {
      constants = LoadConstants(filename);
    }

    /*--- Calculate scaling constants so that zeta -> kroneckor delta for
     * isotropic cells ---*/
    vector<su2double> temp_values = GetEigValues_Q(vector<su2double>(3, 1.0));
    C_zeta = pow(temp_values[0], 0.5);

  } else {
    // Resolution indicator is a variant of RDELTA

    invLengthTensor = new su2double*[nDim];
    for (unsigned int iDim = 0; iDim < nDim; iDim++) {
      invLengthTensor[iDim] = new su2double[nDim];
    }
  }

}

CHybrid_Mediator::~CHybrid_Mediator() {

  delete hybrid_forcing;

  if (Q != NULL) {
    for (unsigned int iDim = 0; iDim < nDim; iDim++)
      if (Q[iDim] != NULL) delete [] Q[iDim];
    delete [] Q;
  }

  if (Qapprox != NULL) {
    for (unsigned int iDim = 0; iDim < nDim; iDim++)
      if (Qapprox[iDim] != NULL) delete [] Qapprox[iDim];
    delete [] Qapprox;
  }

  if (invLengthTensor != NULL) {
    for (unsigned int iDim = 0; iDim < nDim; iDim++)
      if (invLengthTensor[iDim] != NULL) delete [] invLengthTensor[iDim];
    delete [] invLengthTensor;
  }

#ifdef HAVE_LAPACK
  mkl_free_buffers();
#endif
}

void CHybrid_Mediator::SetupRANSNumerics(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         unsigned short iPoint,
                                         unsigned short jPoint) {
  su2double* alpha =
      solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  /*--- Since this is a source term, we don't need a second point ---*/
  rans_numerics->SetHybridParameter(alpha, NULL);

  su2double forcing = hybrid_forcing->GetProduction();
  rans_numerics->SetForcingProduction(forcing);

}

void CHybrid_Mediator::SetupHybridParamSolver(CGeometry* geometry,
                                              CSolver **solver_container,
                                              unsigned short iPoint) {

  unsigned short iDim, jDim, kDim, lDim;
  // XXX: This floor is arbitrary.
  const su2double TKE_MIN = EPS;
  su2double r_k, w_rans;


  /*--- Find eigenvalues and eigenvecs for grid-based resolution tensor ---*/
  su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
  su2double* ResolutionValues = geometry->node[iPoint]->GetResolutionValues();
  su2double** ResolutionVectors = geometry->node[iPoint]->GetResolutionVectors();

  if (config->GetKind_Hybrid_Resolution_Indicator() != RK_INDICATOR) {
    // Compute inverse length scale tensor
    ComputeInvLengthTensor(solver_container[FLOW_SOL]->node[iPoint],
                           solver_container[TURB_SOL]->node[iPoint],
                           solver_container[HYBRID_SOL]->node[iPoint],
                           config->GetKind_Hybrid_Resolution_Indicator());

    vector<su2double> eigvalues_iLM;
    vector<vector<su2double> > eigvectors_iLM;
    SolveGeneralizedEigen(invLengthTensor, ResolutionTensor,
                          eigvalues_iLM, eigvectors_iLM);
    std::vector<su2double>::iterator iter;
    iter = max_element(eigvalues_iLM.begin(), eigvalues_iLM.end());
    unsigned short max_index = distance(eigvalues_iLM.begin(), iter);

    r_k = C_sf*eigvalues_iLM[max_index];

    // NB: with the r_{\Delta} variants of the model, we should never
    // use w_rans.  So, I set it to nan s.t. if it ever gets used
    // inadvertantly, the soln will be obviously wrong.
    // FIXME: RANS weight still necessary?

    w_rans = std::numeric_limits<su2double>::quiet_NaN();
  }
  else if (config->GetKind_Hybrid_Resolution_Indicator() == RK_INDICATOR) {

    /*--- Transform the approximate structure function ---*/
    su2double** PrimVar_Grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
    CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Qapprox);
    vector<vector<su2double> > zeta = BuildZeta(ResolutionValues, ResolutionVectors);
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        Q[iDim][jDim] = 0.0;
        for (kDim = 0; kDim < nDim; kDim++) {
          for (lDim = 0; lDim < nDim; lDim++) {
            // Now zeta*Q*zeta, not Q
            Q[iDim][jDim] += zeta[iDim][kDim] * Qapprox[kDim][lDim] *
                             zeta[lDim][jDim];
          }
        }
      }
    }

    /*--- Find eigenvalues and eigenvectors ---*/
    su2double total_vel_differences = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        total_vel_differences += abs(Q[iDim][jDim]);
      }
    }

    su2double k = solver_container[TURB_SOL]->node[iPoint]->GetSolution(0);
    const su2double MIN_VEL_DIFF = fmax(EPS, EPS*k);
    if (total_vel_differences > MIN_VEL_DIFF) {
      /*--- Only calculate r_k, w_rans if there are resolved velocity differences
       * at resolution scale.  Otherwise, eigenvector calculation is arbitrary */

      /*--- Calculate eigenvectors and eigenvalues of zeta*Q*zeta ---*/
      vector<su2double> eigvalues_zQz;
      vector<vector<su2double> > eigvectors_zQz;
      SolveEigen(Q, eigvalues_zQz, eigvectors_zQz);
      std::vector<su2double>::iterator iter;
      iter = max_element(eigvalues_zQz.begin(), eigvalues_zQz.end());
      unsigned short max_index = distance(eigvalues_zQz.begin(), iter);
      vector<su2double> max_eigenvalue_direction = eigvectors_zQz[max_index];

      /*---Find the largest product of resolved fluctuations at the cutoff---*/
      su2double aniso_ratio = solver_container[TURB_SOL]->node[iPoint]->GetAnisoRatio(); 
      su2double C_kQ = 16.0;
      su2double max_resolved = aniso_ratio*C_kQ*C_sf*TWO3*eigvalues_zQz[max_index];

      /*--- Find the smallest product of unresolved fluctuations at the cutoff ---*/
      su2double min_unresolved;
      switch (config->GetKind_Turb_Model()) {
        case SST: {
          su2double C_mu = 0.22;
          su2double TurbT = solver_container[TURB_SOL]->node[iPoint]->GetTurbTimescale();
          su2double omega = solver_container[TURB_SOL]->node[iPoint]->GetSolution(1);
          min_unresolved = TurbT*k*omega/C_mu;
          break;
        }
        case KE: {
          min_unresolved = solver_container[TURB_SOL]->node[iPoint]->GetSolution(2);
          break;
        }
        default: {
          SU2_MPI::Error("The hybrid mediator is not set up for your turb. model!", CURRENT_FUNCTION);
        }
      }

      /*--- Calculate the resolution adequacy parameter ---*/
      r_k = max(max_resolved / fmax(min_unresolved, TKE_MIN), EPS);


      /*--- Find the dissipation ratio ---*/
      su2double C_eps = 0.25;
      su2double TurbL = solver_container[TURB_SOL]->node[iPoint]->GetTurbLengthscale();
      su2double d_max = GetProjResolution(ResolutionTensor,
                                          eigvectors_zQz[max_index]);
      su2double r_eps = C_eps * pow(r_k, 1.5) * TurbL / d_max;

      /*--- Calculate the RANS weight ---*/
      w_rans = tanh(0.5*pow(fmax(r_eps - 1, 0), 0.25));
    } else {
      r_k = 0.0;
      w_rans = 0.0;
    }
  }
  else {
    SU2_MPI::Error("Unrecognized HYBRID_RESOLUTION_INDICATOR value!", CURRENT_FUNCTION);
  }

  solver_container[HYBRID_SOL]->node[iPoint]->SetResolutionAdequacy(r_k);
  // FIXME: RANS weight still necessary?
  solver_container[HYBRID_SOL]->node[iPoint]->SetRANSWeight(w_rans);

}

void CHybrid_Mediator::SetupHybridParamNumerics(CGeometry* geometry,
                                                CSolver **solver_container,
                                                CNumerics *hybrid_param_numerics,
                                                unsigned short iPoint,
                                                unsigned short jPoint) {

  /*--- Find and store turbulent length and timescales ---*/
  // FIXME: Keep both types of timescales?
  su2double turb_T =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbTimescale();
  su2double turb_L =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbLengthscale();

#ifndef NDEBUG
  if (turb_T <= 0) {
    SU2_MPI::Error("Turbulent timescale was non-positive!", CURRENT_FUNCTION);
  }
  if (turb_L <= 0) {
    SU2_MPI::Error("Turbulent lengthscale was non-positive!", CURRENT_FUNCTION);
  }
#endif

  hybrid_param_numerics->SetTurbTimescale(turb_T);
  hybrid_param_numerics->SetTurbLengthscale(turb_L);

  /*--- Pass resolution adequacy into the numerics object ---*/

  su2double r_k = solver_container[HYBRID_SOL]->node[iPoint]->GetResolutionAdequacy();
  hybrid_param_numerics->SetResolutionAdequacy(r_k);

  /*--- Pass RANS weight into the numerics object ---*/
  // FIXME: RANS weight still necessary?
  su2double w_rans = solver_container[HYBRID_SOL]->node[iPoint]->GetRANSWeight();
  hybrid_param_numerics->SetRANSWeight(w_rans);

  /*--- Pass ratio of P_F_tilde to P_lim into numerics object ---*/
  su2double production_ratio =
      solver_container[HYBRID_SOL]->node[iPoint]->GetForcingRatio();
  hybrid_param_numerics->SetForcingRatio(production_ratio);
}

void CHybrid_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                                             unsigned short iPoint) {
  if (dynamic_cast<CHybrid_Aniso_Q*>(hybrid_anisotropy)) {
    /*--- XXX: Come up with a beter way to do this, so we only pass what's
     * necessary to the anisotropy model. Multiple dispatch? */
    unsigned int iDim, jDim, kDim, lDim;
      /*--- Use the grid-based resolution tensor, not the anisotropy-correcting
       * resolution tensor ---*/
      su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
      su2double** PrimVar_Grad =
            solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
      CalculateApproxStructFunc(ResolutionTensor, PrimVar_Grad, Q);

      su2double total= 0.0;
      for (iDim = 0; iDim < nDim; iDim++) {
        for (jDim = 0; jDim < nDim; jDim++) {
          total += abs(Q[iDim][jDim]);
        }
      }
      if (total < 9*EPS) {
        /*--- If there are negligible velocity differences at the resolution scale,
         * set the tensor to the Kroneckor delta ---*/
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            Q[iDim][jDim] = (iDim == jDim);
          }
        }
      } else {
        /*--- Normalize the approximate structure function tensor ---*/
        su2double norm = 0.0;
        for (iDim = 0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            norm += Q[iDim][jDim]*Q[jDim][iDim];
          }
        }
        norm = sqrt(norm);
        for (iDim=0; iDim < nDim; iDim++) {
          for (jDim = 0; jDim < nDim; jDim++) {
            Q[iDim][jDim] /= norm;
          }
        }
      }

#ifndef NDEBUG
  su2double trace = Q[0][0] + Q[1][1] + Q[2][2];
  if (trace > 1e7) {
    ostringstream error_msg;
    error_msg << "Trace of the stress anisotropy was unusually large." << endl;
    error_msg << "       Trace: " << trace << endl;
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  } else if (trace < 1e-7) {
    ostringstream error_msg;
    error_msg << "Trace of the stress anisotropy was unusually small." << endl;
    error_msg << "       Trace: " << trace << endl;
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }
#endif

      hybrid_anisotropy->SetTensor(Q);

      /*--- Retrieve and pass along the resolution adequacy parameter ---*/

      vector<su2double> scalars;
      scalars.push_back(solver_container[HYBRID_SOL]->node[iPoint]->GetRANSWeight());
      hybrid_anisotropy->SetScalars(scalars);
  }
}


/**
 * \brief The hybrid solver needs the resolution adequacy parameter, which
 *        is dependent on RANS results.
 *
 * \param[in] geometry - A pointer to the geometry
 * \param[in] solver_container - An array of solvers
 * \param[in] iPoint - The node being evaluated
 */
void CHybrid_Mediator::SetupResolvedFlowSolver(CGeometry* geometry,
                                               CSolver **solver_container,
                                               unsigned short iPoint) {
  // TODO: Set up viscous anisotropy here instead

  /*--- Set up forcing ---*/
  SetupForcing(geometry, solver_container, iPoint);
}

void CHybrid_Mediator::SetupForcing(CGeometry* geometry,
                                    CSolver **solver_container,
                                    unsigned short iPoint) {

  /*--- Pull in all the relevant parameters ---*/
  su2double* prim_vars =
      solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
  su2double** prim_var_grad =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  su2double* turb_vars =
      solver_container[TURB_SOL]->node[iPoint]->GetSolution();
  su2double turb_L =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbLengthscale();
  su2double turb_T =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbTimescale();
  su2double* alpha =
      solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  su2double* S_terms =
      solver_container[HYBRID_SOL]->node[iPoint]->GetSourceTerms();
  hybrid_forcing->SetCoord(geometry->node[iPoint]->GetCoord());
  hybrid_forcing->SetPrimitive(prim_vars);
  hybrid_forcing->SetPrimVarGradient(prim_var_grad);
  hybrid_forcing->SetTurbVar(turb_vars);
  hybrid_forcing->SetTurbLengthscale(turb_L);
  hybrid_forcing->SetTurbTimescale(turb_L);
  hybrid_forcing->SetHybridParameter(alpha);
  hybrid_forcing->SetHybridSourceTerms(S_terms);

  /*--- Run calculation now that we've initialized everything ---*/
  hybrid_forcing->CalculateForcing();

  /*--- Pass results on to respective solvers
   * We take over some of the work that would normally be in
   * `SetupRANSSolver` or `SetupHybridSolver` so that the same forcing
   * values are used for all equations.  This has the advantage of only
   * requiring one forcing calculation per iteration. ---*/
  su2double** tau_F = hybrid_forcing->GetStress();
  solver_container[FLOW_SOL]->node[iPoint]->SetForcingStress(tau_F);
  su2double P_F = hybrid_forcing->GetProduction();
  solver_container[TURB_SOL]->node[iPoint]->SetForcingProduction(P_F);
  su2double P_F_ratio = hybrid_forcing->GetForcingRatio();
  solver_container[HYBRID_SOL]->node[iPoint]->SetForcingRatio(P_F_ratio);
}


void CHybrid_Mediator::SetupResolvedFlowNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {

  /*--- Pass alpha to the resolved flow ---*/

  su2double* alpha_i = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  su2double* alpha_j = solver_container[HYBRID_SOL]->node[jPoint]->GetSolution();
  visc_numerics->SetHybridParameter(alpha_i, alpha_j);

  /*--- Pass the stress anisotropy tensor to the resolved flow ---*/

  su2double** aniso_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscAnisotropy();
  su2double** aniso_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscAnisotropy();
  visc_numerics->SetEddyViscAnisotropy(aniso_i, aniso_j);

  /*--- Pass the forcing stress tensor to the resolved flow ---*/

  su2double** tau_F_i = solver_container[FLOW_SOL]->node[iPoint]->GetForcingStress();
  su2double** tau_F_j = solver_container[FLOW_SOL]->node[iPoint]->GetForcingStress();
  visc_numerics->SetForcingStress(tau_F_i, tau_F_j);

}


void CHybrid_Mediator::ComputeInvLengthTensor(CVariable* flow_vars,
                                              CVariable* turb_vars,
                                              CVariable* hybr_vars,
                                              int short hybrid_res_ind) {

  unsigned short iDim, jDim, kDim;
  su2double Sd[3][3], Om[3][3], delta[3][3], Pij[3][3];
  su2double div_vel;

  // Not intended for use in 2D
  if (nDim == 2) {
    SU2_MPI::Error("The RDELTA resolution adequacy option is not implemented for 2D!", CURRENT_FUNCTION);
  }

  // Get primative variables and gradients
  su2double** val_gradprimvar =  flow_vars->GetGradient_Primitive();

  // Get eddy viscosity
  su2double eddy_viscosity = flow_vars->GetEddyViscosity();

  // Get hybrid blend (i.e., alpha)
  su2double* blend = hybr_vars->GetSolution();

  // Bound alpha away from zero
  su2double alpha = max(blend[0],1e-8);

  // delta tensor (for convenience)
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (iDim==jDim) {
	delta[iDim][jDim] = 1.0;
      } else {
	delta[iDim][jDim] = 0.0;
      }
    }
  }

  // Compute divergence
  div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  // The deviatoric part of the strain rate tensor
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Sd[iDim][jDim] = 0.5*(val_gradprimvar[iDim+1][jDim] + val_gradprimvar[jDim+1][iDim]) -
                       delta[iDim][jDim]*div_vel/3.0;
    }
  }

  // The rotation rate tensor
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Om[iDim][jDim] = 0.5*(val_gradprimvar[iDim+1][jDim] - val_gradprimvar[jDim+1][iDim]);
    }
  }

  // v2 here is *subgrid*, so must multiply by alpha
  su2double v2;
  if (config->GetKind_Turb_Model() == KE) {
    v2 = alpha*max(turb_vars->GetSolution(2),1e-8);
  } else if (config->GetKind_Turb_Model() == SST) {
    v2 = alpha*2.0*max(turb_vars->GetSolution(0),1e-8)/3.0;
  } else {
    SU2_MPI::Error("The RDELTA resolution adequacy option is only implemented for KE and SST turbulence models!", CURRENT_FUNCTION);
  }

  // NB: Assumes isotropic eddy viscosity
  // TODO: If/when we go back to anisotropic eddy viscosity, make sure
  // change propagates to here

  // Strain only part
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Pij[iDim][jDim] = 0.0;
      for (kDim = 0; kDim < nDim; kDim++) {
        // FIXME: Not true...Check that proper eddy viscosity is used.
        // NB: Assumes that eddy_viscosity is *full* eddy
        // viscosity---i.e., not multiplied by alpha^2---so that
        // alpha^2 is necessary here
        Pij[iDim][jDim] += 2.0*alpha*alpha*eddy_viscosity*Sd[iDim][kDim]*Sd[kDim][jDim];
      }
    }
  }

  switch (hybrid_res_ind) {
  case RDELTA_INDICATOR_STRAIN_ONLY:
    // nothing to do here
    break;
  case RDELTA_INDICATOR_FULLP:
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        for (kDim = 0; kDim < nDim; kDim++) {
          Pij[iDim][jDim] += 2.0*alpha*alpha*eddy_viscosity*Sd[iDim][kDim]*Om[jDim][kDim];
        }
      }
    }
    break;
  case RDELTA_INDICATOR_FULLP_VELCON:
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        for (kDim = 0; kDim < nDim; kDim++) {
          Pij[iDim][jDim] += 2.0*alpha*alpha*eddy_viscosity*Sd[iDim][kDim]*Om[kDim][jDim];
        }
      }
    }
    break;
  default:
    SU2_MPI::Error("Unrecognized RDELTA option.", CURRENT_FUNCTION);
    break;
  }


  // set inverse length scale tensor
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      invLengthTensor[iDim][jDim] = 0.5*(Pij[iDim][jDim] + Pij[jDim][iDim]) / (v2*sqrt(v2)) ;
    }
  }

#ifndef NDEBUG
  // check for nans
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (invLengthTensor[iDim][jDim] != invLengthTensor[iDim][jDim]) {
        SU2_MPI::Error("invLengthTensor has NaN!", CURRENT_FUNCTION);
      }
    }
  }
#endif

}


su2double CHybrid_Mediator::GetProjResolution(su2double** resolution_tensor,
                                              vector<su2double> direction) {

#ifndef NDEBUG
  su2double magnitude_squared = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    magnitude_squared += direction[iDim]*direction[iDim];
  if (abs(magnitude_squared - 1.0) > 1e-7) {
    ostringstream error_msg;
    error_msg << "The unit vector for the projected resolution calc had a ";
    error_msg << "magnitude greater than 1!" << endl;
    error_msg << "    Magnitude: " << sqrt(magnitude_squared);
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }
#endif

  su2double temp, result = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    temp = 0;
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      temp += resolution_tensor[iDim][jDim]*direction[jDim];
    }
    result += temp*temp;
  }
  return sqrt(result);
}

vector<vector<su2double> > CHybrid_Mediator::LoadConstants(string filename) {
  vector<vector<su2double> > output;
  output.resize(nDim);
  ifstream file;
  for (int i=0; i<nDim; i++) {
    stringstream ss;
    ss << filename << i << ".dat";
    string fullname = ss.str();
    file.open(fullname.c_str());
    if (file.is_open()) {
      su2double value;
      while (file >> value) {
        output[i].push_back(value);
      }
      file.close();
    } else {
      ostringstream error_msg;
      error_msg << "Could not open the hybrid constants file." << endl;
      error_msg << "       Tried reading file " << fullname;
      SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
    }
  }
  return output;
};

vector<su2double> CHybrid_Mediator::GetEigValues_Q(vector<su2double> eigvalues_M) {
  su2double dnorm = *min_element(eigvalues_M.begin(), eigvalues_M.end());

  /*--- Normalize eigenvalues ---*/
  su2double a, b;
  if (eigvalues_M[0] == dnorm) {
    a = eigvalues_M[1]/dnorm;
    b = eigvalues_M[2]/dnorm;
  } else if (eigvalues_M[1] == dnorm) {
    a = eigvalues_M[0]/dnorm;
    b = eigvalues_M[2]/dnorm;
  } else {
    a = eigvalues_M[0]/dnorm;
    b = eigvalues_M[1]/dnorm;
  }
#ifndef NDEBUG
  if (a < 1 || b < 1) {
    SU2_MPI::Error("Normalization in the zeta transformation failed!", CURRENT_FUNCTION);
  }
#endif

  /*--- Convert to cylindrical coordinates ---*/
  su2double r = sqrt(a*a + b*b);
  su2double theta = acos(max(a,b)/r);

  /*--- Convert to more convenient log coordinates ---*/
  su2double y = log(sin(2*theta));
  su2double x = log(r);

  vector<su2double> g(3);
  for (int iDim = 0; iDim < nDim; iDim++) {
    g[iDim] = constants[iDim][0];
    g[iDim] += constants[iDim][1]*x;
    g[iDim] += constants[iDim][2]*y;
    g[iDim] += constants[iDim][3]*x*x;
    g[iDim] += constants[iDim][4]*x*y;
    g[iDim] += constants[iDim][5]*y*y;
    g[iDim] += constants[iDim][6]*x*x*x;
    g[iDim] += constants[iDim][7]*x*x*y;
    g[iDim] += constants[iDim][8]*x*y*y;
    g[iDim] += constants[iDim][9]*y*y*y;
    g[iDim] += constants[iDim][10]*x*x*x*x;
    g[iDim] += constants[iDim][11]*x*x*x*y;
    g[iDim] += constants[iDim][12]*x*x*y*y;
    g[iDim] += constants[iDim][13]*x*y*y*y;
    g[iDim] += constants[iDim][14]*y*y*y*y;
  }

  vector<su2double> eigvalues_Q(3);
  for (int iDim = 0; iDim < nDim; iDim++) {
    su2double d_in = log(eigvalues_M[iDim]/dnorm);
    eigvalues_Q[iDim] = g[0] + g[1]*d_in + g[2]*pow(d_in, 2);
    eigvalues_Q[iDim] = exp(eigvalues_Q[iDim]);
  }

  return eigvalues_Q;
}

vector<vector<su2double> > CHybrid_Mediator::BuildZeta(su2double* values_M,
                                                       su2double** vectors_M) {

  vector<vector<su2double> > zeta(3, vector<su2double>(3,0));

#ifdef HAVE_LAPACK
  unsigned short iDim, jDim, kDim, lDim;

  vector<su2double> eigvalues_M;
  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues_M.push_back(values_M[iDim]);

  /*--- Solve for the modified resolution tensor  ---*/
  vector<su2double> eigvalues_Zeta = GetEigValues_Zeta(eigvalues_M);
  vector<vector<su2double> > temp(3, vector<su2double>(3));
  for (iDim = 0; iDim < nDim; iDim++) {
    temp[iDim][iDim] = eigvalues_Zeta[iDim];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      zeta[iDim][jDim] = 0.0;
      for (kDim = 0; kDim < nDim; kDim++) {
        for (lDim = 0; lDim < nDim; lDim++) {
          zeta[iDim][jDim] += vectors_M[kDim][iDim] * temp[kDim][lDim] *
                              vectors_M[lDim][jDim];
        }
      }
    }
  }
#else
  SU2_MPI::Error("Eigensolver without LAPACK not implemented; please use LAPACK.", CURRENT_FUNCTION);
#endif

  /*--- This return is placed here to avoid "no return" compiler warnings.---*/
  return zeta;

}

vector<su2double> CHybrid_Mediator::GetEigValues_Zeta(vector<su2double> eigvalues_M) {
  /*--- Find the minimum eigenvalue ---*/

  su2double dnorm = *min_element(eigvalues_M.begin(), eigvalues_M.end());

  /*--- Get eigenvalues to the normalized gradient-gradien tensor ---*/

  vector<su2double> eigvalues_Q = GetEigValues_Q(eigvalues_M);
#ifndef NDEBUG
  for (int iDim = 0; iDim < nDim; iDim++) {
    if (eigvalues_Q[iDim] != eigvalues_Q[iDim]) {
      SU2_MPI::Error("At least one computed eigenvalue was NaN!", CURRENT_FUNCTION);
    }
  }
#endif

  /*--- Use eigenvalues from M and G to find eigenvalues for modified M ---*/

  vector<su2double> eigvalues_zeta(3);
  for (int iDim = 0; iDim < nDim; iDim++) {
    eigvalues_zeta[iDim] = C_zeta*pow((eigvalues_M[iDim]/dnorm),1.0/3)*
                           pow(eigvalues_Q[iDim],-0.5);

    // XXX: Numerical fit for anisotropic resolution doesn't include > 256
    if (eigvalues_M[iDim]/dnorm > 256) {
      eigvalues_zeta[iDim] = max(eigvalues_zeta[iDim], 0.90);
    }
  }

  return eigvalues_zeta;
}

void CHybrid_Mediator::SolveEigen(su2double** M,
                                  vector<su2double> &eigvalues,
                                  vector<vector<su2double> > &eigvectors) {
unsigned short iDim, jDim;

#ifndef NDEBUG
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (M[iDim][jDim] != M[iDim][jDim]) {
        SU2_MPI::Error("SolveEigen received a matrix with NaN!", CURRENT_FUNCTION);
      }
    }
  }
  su2double sum = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      sum += abs(M[iDim][jDim]);
    }
  }
  if (sum < EPS) {
    SU2_MPI::Error("SolveEigen received an empty matrix!", CURRENT_FUNCTION);
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (iDim == jDim) continue;
      if (abs(M[iDim][jDim] - M[jDim][iDim]) > 1e-10 &&
          abs(M[iDim][jDim] - M[jDim][iDim])/abs(M[iDim][jDim]) > 1e-6) {
        ostringstream error_msg;
        error_msg << "SolveEigen received a non-symmetric matrix!" << endl;
        error_msg << "    The difference between elements" << endl;
        error_msg << "      [" << iDim << ", " << jDim << "] and [" << jDim << ", " << iDim << "]" << endl;
        error_msg << "      was: " << M[iDim][jDim] - M[jDim][iDim] << endl;
        error_msg << "    Matrix:" << endl;
        error_msg << "      [[" << M[0][0] << ", " << M[0][1] << ", " << M[0][2] << "]" << endl;
        error_msg << "       [" << M[1][0] << ", " << M[1][1] << ", " << M[1][2] << "]" << endl;
        error_msg << "       [" << M[2][0] << ", " << M[2][1] << ", " << M[2][2] << "]]";
        SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
      }
    }
  }
#endif

  eigvalues.resize(nDim);
  eigvectors.resize(nDim, std::vector<su2double>(nDim));

#ifdef HAVE_LAPACK
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim< nDim; jDim++) {
      mat[iDim*nDim+jDim] = M[iDim][jDim];
    }
  }

  /*--- Call LAPACK routines ---*/

  info = LAPACKE_dsyev(LAPACK_ROW_MAJOR, 'V', 'U', nDim, mat, nDim, eigval);
  if (info != 0) {
    ostringstream error_msg;
    error_msg << "The solver failed to compute eigenvalues." << endl;
    error_msg << "The info variable was set to: " << info;
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }

  /*--- Rewrite arrays to eigenvalues output ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues[iDim] = eigval[iDim];

  /*--- Check the values ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    if (eigval[iDim] < 0.0 && eigval[iDim] > -1e-4) {
        eigvalues[iDim] = 0.0;
    } else if (eigval[iDim] < -1e-4) {
      ostringstream error_msg;
      error_msg << "The solver returned a large negative eigenvalue!" << endl;
      error_msg << "    Eigenvalues: [";
      error_msg << eigval[0] << ", ";
      error_msg << eigval[1] << ", ";
      error_msg << eigval[2] << "]" << endl;
      SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
    }
  }

//  su2double max_val = max(eigval[0], max(eigval[1], eigval[2]));
//  for (iDim = 0; iDim < nDim; iDim++) {
//    if (eigvalues[iDim] > 0 && log10(max_val/eigvalues[iDim]) > 10) {
//      // If the condition number is bad, give up
//      eigvalues[iDim] = 0.0;
//    }
//  }

  /*--- Normalize the eigenvectors by the L2 norm of each vector ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eigvectors[iDim][jDim] = mat[jDim*nDim+iDim];
    }
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    su2double norm = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      norm += eigvectors[iDim][jDim]*eigvectors[iDim][jDim];
    }
    norm = sqrt(norm);
    for (jDim = 0; jDim < nDim; jDim++) {
      eigvectors[iDim][jDim] /= norm;
    }
  }
#else
  SU2_MPI::Error("Eigensolver without LAPACK not implemented; please use LAPACK.", CURRENT_FUNCTION);
#endif
}

void CHybrid_Mediator::SolveGeneralizedEigen(su2double** A, su2double** B,
					     vector<su2double> &eigvalues,
					     vector<vector<su2double> > &eigvectors) {

  eigvalues.resize(nDim);
  eigvectors.resize(nDim, std::vector<su2double>(nDim));

#ifdef HAVE_LAPACK
  unsigned short iDim, jDim;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim< nDim; jDim++) {
      mat[iDim*nDim+jDim]  = A[iDim][jDim];
      matb[iDim*nDim+jDim] = B[iDim][jDim];
    }
  }

  /*--- Call LAPACK routines ---*/

  // itype = 2 ==> we're solving A*B*v = \lambda v
  // when A = L^{-1} and B = M, this is what we want
  info = LAPACKE_dsygv(LAPACK_ROW_MAJOR, 2, 'V', 'U',
                       nDim, mat, nDim, matb, nDim, eigval);
  if (info != 0) {
    ostringstream error_msg;
    error_msg << "The generalized eigensolver failed with info = "
         << info << ".";
    SU2_MPI::Error(error_msg.str(), CURRENT_FUNCTION);
  }

  /*--- Rewrite arrays to eigenvalues output ---*/
  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues[iDim] = eigval[iDim];

  /*--- Check the values ---*/
  // NB: Here, B is SPD and A is assumed symmetric but not necessarily
  // positive definite.  Thus, don't check positivity of eigenvalues.

  /*--- Normalize the eigenvectors by the L2 norm of each vector ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      eigvectors[iDim][jDim] = mat[jDim*nDim+iDim];
    }
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    su2double norm = 0.0;
    for (jDim = 0; jDim < nDim; jDim++) {
      norm += eigvectors[iDim][jDim]*eigvectors[iDim][jDim];
    }
    norm = sqrt(norm);
    for (jDim = 0; jDim < nDim; jDim++) {
      eigvectors[iDim][jDim] /= norm;
    }
  }
#else
  SU2_MPI::Error("Eigensolver without LAPACK not implemented; please use LAPACK.", CURRENT_FUNCTION);
#endif

}


vector<vector<su2double> > CHybrid_Mediator::GetConstants() {
  return constants;
}



void CHybrid_Mediator::AdjustEddyViscosity(CSolver **solver_container,
                                           unsigned short iPoint,
                                           su2double *eddy_viscosity) {
  su2double* alpha = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  *eddy_viscosity *= alpha[0]*alpha[0];
}



CHybrid_Dummy_Mediator::CHybrid_Dummy_Mediator(unsigned short nDim,
                                               CConfig* config)
    : nDim(nDim) {

  /*--- Set the default value of the hybrid parameter ---*/
  dummy_alpha = new su2double[1];
  dummy_alpha[0] = 1.0;
}

CHybrid_Dummy_Mediator::~CHybrid_Dummy_Mediator() {
  delete [] dummy_alpha;
}

void CHybrid_Dummy_Mediator::SetupResolvedFlowSolver(CGeometry* geometry,
                                                     CSolver **solver_container,
                                                     unsigned short iPoint) {
}

void CHybrid_Dummy_Mediator::SetupRANSNumerics(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         unsigned short iPoint,
                                         unsigned short jPoint) {
  // This next line is just here for testing purposes.
  su2double* alpha = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  assert(alpha != NULL);
  rans_numerics->SetHybridParameter(dummy_alpha, dummy_alpha);
}

void CHybrid_Dummy_Mediator::SetupHybridParamSolver(CGeometry* geometry,
                                           CSolver **solver_container,
                                           unsigned short iPoint) {
}

void CHybrid_Dummy_Mediator::SetupHybridParamNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics *hybrid_param_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {
}

void CHybrid_Dummy_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                                             unsigned short iPoint) {
}

void CHybrid_Dummy_Mediator::SetupResolvedFlowNumerics(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             unsigned short iPoint,
                                             unsigned short jPoint) {

  /*--- Pass alpha to the resolved flow ---*/

  // This next two lines are just here for testing purposes.
  su2double* alpha_i = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  su2double* alpha_j = solver_container[HYBRID_SOL]->node[jPoint]->GetSolution();
  assert(alpha_i != NULL);
  assert(alpha_j != NULL);
  visc_numerics->SetHybridParameter(dummy_alpha, dummy_alpha);

  /*--- Pass the stress anisotropy tensor to the resolved flow ---*/

  su2double** aniso_i = solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscAnisotropy();
  su2double** aniso_j = solver_container[FLOW_SOL]->node[jPoint]->GetEddyViscAnisotropy();
  visc_numerics->SetEddyViscAnisotropy(aniso_i, aniso_j);
}

void CHybrid_Dummy_Mediator::AdjustEddyViscosity(CSolver **solver_container,
                                                 unsigned short iPoint,
                                                 su2double *eddy_viscosity) {
  /*--- Don't adjust the eddy viscosity for a no-hybrid model.
   * Technically, if alpha=1 everywhere, we could still multiply eddy
   * viscosity by alpha.  But there's no point in checking multiplication ---*/
}

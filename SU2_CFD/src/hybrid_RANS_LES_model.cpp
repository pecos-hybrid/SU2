/*!
 * \file hybrid_RANS_LES_model.cpp
 * \brief Describes the hybrid RANS/LES models
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
#include "../include/numerics_direct_mean_hybrid.hpp"
#include "../include/solver_structure.hpp"
#ifdef HAVE_LAPACK
#include "mkl.h"
#include "mkl_lapacke.h"
#endif

CHybrid_Mediator::CHybrid_Mediator(unsigned short nDim, CConfig* config,
                                   const string& filename)
    : nDim(nDim), fluct_stress_model(NULL), forcing_model(NULL), config(config) {

  int rank = MASTER_NODE;
#ifdef HAVE_MPI
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
#endif

  aniso_eddy_viscosity = new su2double*[nDim];
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    aniso_eddy_viscosity[iDim] = new su2double[nDim];
  }

  /*--- Allocate the inverse length tensor (used in calcs) ---*/

  invLengthTensor = new su2double*[nDim];
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    invLengthTensor[iDim] = new su2double[nDim];
  }
}

CHybrid_Mediator::~CHybrid_Mediator() {
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    delete [] invLengthTensor[iDim];
    delete [] aniso_eddy_viscosity[iDim];
  }
  delete [] invLengthTensor;
  delete [] aniso_eddy_viscosity;

  /*--- Through dependency injection, we've taken over responsibility for
   * the fluctuating stress and forcing models.
   */
  if (fluct_stress_model != NULL) {
    delete fluct_stress_model;
  }

  if (forcing_model != NULL) {
    delete forcing_model;
  }

#ifdef HAVE_LAPACK
  mkl_free_buffers();
#endif
}

void CHybrid_Mediator::SetupRANSNumerics(CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         const unsigned long iPoint) {

  CVariable** flow_node = solver_container[FLOW_SOL]->average_node;
  rans_numerics->SetKineticEnergyRatio(flow_node[iPoint]->GetKineticEnergyRatio());
  if (config->GetUse_Resolved_Turb_Stress()) {
    rans_numerics->SetResolvedTurbStress(flow_node[iPoint]->GetResolvedTurbStress());
  } else {
    rans_numerics->SetProduction(flow_node[iPoint]->GetProduction());
  }
}

void CHybrid_Mediator::ComputeResolutionAdequacy(const CGeometry* geometry,
                                                 CSolver **solver_container,
                                                 unsigned long iPoint) {

  unsigned short iDim, jDim, kDim, lDim;
  // XXX: This floor is arbitrary.
  const su2double TKE_MIN = EPS;

  /*--- Find eigenvalues and eigenvecs for grid-based resolution tensor ---*/
  const su2double* const* ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();

  // Compute inverse length scale tensor
  const su2double alpha =
      solver_container[FLOW_SOL]->average_node[iPoint]->GetKineticEnergyRatio();
  ComputeInvLengthTensor(solver_container[FLOW_SOL]->node[iPoint],
                         solver_container[FLOW_SOL]->average_node[iPoint],
                         solver_container[TURB_SOL]->node[iPoint],
                         alpha,
                         config->GetKind_Hybrid_Resolution_Indicator());

  vector<su2double> eigvalues_iLM;
  vector<vector<su2double> > eigvectors_iLM;
  SolveGeneralizedEigen(invLengthTensor, ResolutionTensor,
                        eigvalues_iLM, eigvectors_iLM);

  su2double max_eigval = 0.0;
  for (iDim=0; iDim<3; iDim++) {
    if( abs(eigvalues_iLM[iDim]) > max_eigval) {
      max_eigval = abs(eigvalues_iLM[iDim]);
    }
  }

  const su2double C_r = 1.0;
  const su2double r_k_min = 1.0E-8;
  const su2double r_k_max = (alpha > 1) ? 1.0 : 30;
  const su2double r_k = max(min(C_r*max_eigval, r_k_max), r_k_min);

  // Set resolution adequacy in the CNSVariables class
  // Use the resolved (e.g. instantaneous) variables, instead of the
  // average variables.  Although r_k uses some average variables, such
  // as the RANS turbulent stress, the instantaneous resolution adequacy
  // is still different than the time-averaged resolution adequacy.
  solver_container[FLOW_SOL]->node[iPoint]->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::ComputeForcingField(CSolver** solver,
                                           CGeometry *geometry,
                                           CConfig *config) {
  forcing_model->ComputeForcingField(solver, geometry, config);
}


void CHybrid_Mediator::SetupResolvedFlowSolver(const CGeometry* geometry,
                                               CSolver **solver_container,
                                               const unsigned long iPoint) {

  if (fluct_stress_model) {

    const su2double* primvar =
        solver_container[FLOW_SOL]->average_node[iPoint]->GetPrimitive();
    //solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    const su2double* turbvar =
        solver_container[TURB_SOL]->node[iPoint]->GetSolution();
    const su2double mean_eddy_visc =
        solver_container[FLOW_SOL]->node[iPoint]->GetEddyViscosity();

    fluct_stress_model->SetPrimitive(primvar);
    fluct_stress_model->SetTurbVar(turbvar);
    fluct_stress_model->CalculateEddyViscosity(geometry, config, iPoint,
                                               mean_eddy_visc,
                                               aniso_eddy_viscosity);

    /*--- XXX: This is an ad-hoc correction
     * Rescale the eddy viscosity to rapidly remove fluctuations where
     * resolved turbulence has been transported into regions with
     * insufficient resolution.
     *
     * This rescaling has been temporarily removed to improve the
     * of the model for fully-developed channel flow problems.  It is
     * unclear how necessary it is more more complex problems. ---*/

    // const su2double avg_resolution_adequacy =
    //   solver_container[FLOW_SOL]->average_node[iPoint]->GetResolutionAdequacy();
    // assert(avg_resolution_adequacy >= 0);
    // assert(avg_resolution_adequacy == avg_resolution_adequacy);
    // //const su2double factor = pow(min(avg_resolution_adequacy, 10.0), 4.0/3);
    // const su2double factor = max(min(avg_resolution_adequacy, 10.0), 1.0);
    // for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    //   for (unsigned short jDim = 0; jDim < nDim; jDim++) {
    //      aniso_eddy_viscosity[iDim][jDim] *= factor;
    //   }
    // }

    solver_container[FLOW_SOL]->node[iPoint]->SetAnisoEddyViscosity(aniso_eddy_viscosity);

  }
}

void CHybrid_Mediator::SetupResolvedFlowNumerics(const CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             const unsigned long iPoint,
                                             const unsigned long jPoint) {

  CAvgGrad_Hybrid* numerics = dynamic_cast<CAvgGrad_Hybrid*>(visc_numerics);

  /*--- We assume that the averaging start time from the cfg file is
   * dimensional ---*/
  const su2double time = config->GetCurrent_UnstTime();

  if (time > config->GetAveragingStartTime()) {
    su2double* primvar_i =
        solver_container[FLOW_SOL]->average_node[iPoint]->GetPrimitive();
    su2double* primvar_j =
        solver_container[FLOW_SOL]->average_node[jPoint]->GetPrimitive();

    su2double** primvar_grad_i =
        solver_container[FLOW_SOL]->average_node[iPoint]->GetGradient_Primitive();
    su2double** primvar_grad_j =
        solver_container[FLOW_SOL]->average_node[jPoint]->GetGradient_Primitive();

    assert(primvar_i != NULL);
    assert(primvar_j != NULL);

    numerics->SetPrimitive_Average(primvar_i, primvar_j);
    numerics->SetPrimVarGradient_Average(primvar_grad_i, primvar_grad_j);

  } else {

    /*--- Since averages are only calculated once per timestep, prevent
     * the averages from lagging behind the fluctuating quantities.
     * Otherwise, nonzero "fluctuations" will appear. ---*/

    su2double* primvar_i =
        solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
    su2double* primvar_j =
        solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();

    su2double** primvar_grad_i =
        solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
    su2double** primvar_grad_j =
        solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();


    numerics->SetPrimitive_Average(primvar_i, primvar_j);
    numerics->SetPrimVarGradient_Average(primvar_grad_i, primvar_grad_j);

  }

  su2double** aniso_viscosity_i =
      solver_container[FLOW_SOL]->node[iPoint]->GetAnisoEddyViscosity();
  su2double** aniso_viscosity_j =
      solver_container[FLOW_SOL]->node[jPoint]->GetAnisoEddyViscosity();

  const su2double alpha_i =
      solver_container[FLOW_SOL]->average_node[iPoint]->GetKineticEnergyRatio();
  const su2double alpha_j =
      solver_container[FLOW_SOL]->average_node[jPoint]->GetKineticEnergyRatio();

  numerics->SetAniso_Eddy_Viscosity(aniso_viscosity_i, aniso_viscosity_j);
  numerics->SetKineticEnergyRatio(alpha_i, alpha_j);
}


void CHybrid_Mediator::ComputeInvLengthTensor(CVariable* flow_vars,
                                              CVariable* flow_avgs,
                                              CVariable* turb_vars,
                                              const su2double val_alpha,
                                              const int short hybrid_res_ind) {

  unsigned short iDim, jDim, kDim;
  su2double Sd[3][3], Sd_avg[3][3], Om[3][3], delta[3][3], Pij[3][3];
  su2double Gt[3][3], Gp[3][3], Gpd[3][3], tauSGET[3][3];
  su2double div_vel, div_avg_vel, div_fluct;

  // Not intended for use in 2D
  if (nDim == 2) {
    SU2_MPI::Error("The RDELTA resolution adequacy option is not implemented for 2D!", CURRENT_FUNCTION);
  }

  // 1) Preliminaries: Extract data, convenience calcs...

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

  const su2double rho =  flow_avgs->GetDensity();

  su2double alpha = max(val_alpha, 1e-8);

  // Get primative gradients
  su2double** val_gradprimvar =  flow_vars->GetGradient_Primitive();
  su2double** val_gradprimavg =  flow_avgs->GetGradient_Primitive();

  // Get eddy viscosities
  /*--- Use the turb_vars (rather than flow_vars) for the eddy viscosity.
   * When the eddy viscosity is computed, it is immediately stored in the
   * turbulence solver nodes. It is only copied over to the flow solver
   * at the start of the next iteration.  So retrieving the eddy viscosity
   * from the turbulence solver ensures that the most up-to-date value is
   * used, rather than a time-lagged value. ---*/
  su2double eddy_viscosity = turb_vars->GetmuT();
  su2double** aniso_viscosity = flow_vars->GetAnisoEddyViscosity();

  // Compute the fluctuating velocity gradient tensor
  for (iDim=0; iDim<nDim; iDim++) {
    for (jDim=0; jDim<nDim; jDim++) {
      Gt[iDim][jDim] = val_gradprimvar[iDim+1][jDim];
      Gp[iDim][jDim] = val_gradprimvar[iDim+1][jDim] - val_gradprimavg[iDim+1][jDim];
    }
  }

  // Compute divergences
  div_vel = div_avg_vel = div_fluct = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++) {
    div_vel     += val_gradprimvar[iDim+1][iDim];
    div_avg_vel += val_gradprimavg[iDim+1][iDim];
    div_fluct   += Gp[iDim][iDim];
  }

  // Evaluate deviatoric part of fluctuating velocity gradient
  for (iDim =0; iDim < nDim; iDim++) {
    for (jDim =0; jDim < nDim; jDim++) {
      Gpd[iDim][jDim] = Gp[iDim][jDim] - delta[iDim][jDim]*div_fluct/3;
    }
  }

  // The deviatoric part of the strain rate tensors (mean and resolved)
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Sd[iDim][jDim]     = 0.5*(val_gradprimvar[iDim+1][jDim] + val_gradprimvar[jDim+1][iDim]) -
                           delta[iDim][jDim]*div_vel/3.0;
      Sd_avg[iDim][jDim] = 0.5*(val_gradprimavg[iDim+1][jDim] + val_gradprimavg[jDim+1][iDim]) -
                           delta[iDim][jDim]*div_avg_vel/3.0;
    }
  }

  // The rotation rate tensor (resolved)
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Om[iDim][jDim] = 0.5*(val_gradprimvar[iDim+1][jDim] - val_gradprimvar[jDim+1][iDim]);
    }
  }


  // Evaluate tauSGET
  for (iDim =0; iDim < nDim; iDim++) {
    for (jDim =0; jDim < nDim; jDim++) {
      tauSGET[iDim][jDim] = 0;
      for (kDim =0; kDim < nDim; kDim++) {
        tauSGET[iDim][jDim] += aniso_viscosity[iDim][kDim]*Gpd[jDim][kDim] +
                               aniso_viscosity[jDim][kDim]*Gpd[iDim][kDim];
      }
    }
  }

  // Make tauSGET trace = 0
  su2double traceTauSGET = 0;
  for (iDim =0; iDim < nDim; iDim++) {
    traceTauSGET += tauSGET[iDim][iDim];
  }
  for (iDim = 0; iDim < nDim; iDim++) {
    tauSGET[iDim][iDim] -= traceTauSGET/3;
  }


  const su2double ktot = max(turb_vars->GetSolution(0),1e-8);

  // v2 here is *subgrid*, so must multiply by alpha
  su2double v2_sgs;
  if (config->GetKind_Turb_Model() == KE) {
    v2_sgs = alpha*max(turb_vars->GetSolution(2),1e-8);
  } else if (config->GetKind_Turb_Model() == SST) {
    v2_sgs = alpha*2.0*max(turb_vars->GetSolution(0),1e-8)/3.0;
  } else {
    SU2_MPI::Error("The RDELTA resolution adequacy option is only implemented for KE and SST turbulence models!", CURRENT_FUNCTION);
  }

  // 2) tauSGRS contribution.  NB: Neglecting divergence contribution
  // here.  TODO: Add divergence contribution.

  /*--- Testing on the WMH indicates that scaling the whole stress by
   * alpha*(2-alpha) improves the model performance.  That change would
   * make the turbulent kinetic energy inconsistent, so it is avoided here.
   * But that indicates there's some other issue. ---*/

  su2double alpha_fac = alpha*(2.0 - alpha);
  alpha_fac = max(alpha_fac, 1e-8);
  alpha_fac = min(alpha_fac, 1.0);

  // Strain only part
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Pij[iDim][jDim] = 0.0;
      for (kDim = 0; kDim < nDim; kDim++) {
        // NB: Assumes that eddy_viscosity is *full* eddy
        // viscosity---i.e., not multiplied by alpha---so that
        // alpha is necessary here
        Pij[iDim][jDim] += 2.0*alpha_fac*eddy_viscosity*Sd_avg[iDim][kDim]*Sd[kDim][jDim];
      }
    }
  }

  switch (hybrid_res_ind) {
  case RDELTA_INDICATOR_FULLP_WITH_DIV:
    /*--- It's currently unclear whether the divergence of the velocity
     * should be included in the subgrid production.  If it is included,
     * negative eigenvalues could result, which would be interpreted as
     * "negative" lengthscales. ---*/
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        for (kDim = 0; kDim < nDim; kDim++) {
          // rotation contribution
          Pij[iDim][jDim] += 2.0*alpha_fac*eddy_viscosity*Sd_avg[iDim][kDim]*Om[kDim][jDim];
          // Isotropic contribution (from delta_ij div_vel)
          Pij[iDim][jDim] += 2.0*alpha_fac*eddy_viscosity*Sd_avg[iDim][jDim]*div_vel/3.0;
        }
        // rho*k contribtuion
	Pij[iDim][jDim] -= 2.0/3.0*alpha*rho*ktot *
            (Sd[iDim][jDim] + Om[iDim][jDim] + div_vel*delta[iDim][jDim]/3.0);
      }
    }
    break;
  case RDELTA_INDICATOR_FULLP_VELCON:
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        for (kDim = 0; kDim < nDim; kDim++) {
          // rotation contribution
          Pij[iDim][jDim] += 2.0*alpha_fac*eddy_viscosity*Sd_avg[iDim][kDim]*Om[kDim][jDim];
          // Contribution from div of velocity is (intentionally) omitted
        }
        // rho*k contribtuion
	Pij[iDim][jDim] -= 2.0*alpha*rho*ktot*(Sd[iDim][jDim]+Om[iDim][jDim])/3.0;
      }
    }
    break;
  default:
    SU2_MPI::Error("Unrecognized RDELTA option.", CURRENT_FUNCTION);
    break;
  }

  // 3) tauSGET contribution
  
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      for (kDim = 0; kDim < nDim; kDim++) {
        Pij[iDim][jDim] += tauSGET[iDim][kDim]*Gt[kDim][jDim];
      }
    }
  }


  // 4) Compute inverse length scale tensor from production tensor
  // NB: Typo in AIAA paper
  const su2double t0 = 1.5*sqrt(1.5);
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      invLengthTensor[iDim][jDim] =
        0.5*(Pij[iDim][jDim] + Pij[jDim][iDim]) / (t0*rho*v2_sgs*sqrt(v2_sgs));
    }
  }

#ifndef NDEBUG
  // check for nans
  bool found_nan = false;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (invLengthTensor[iDim][jDim] != invLengthTensor[iDim][jDim]) {
        found_nan = true;
      }
    }
  }
  if (found_nan) {
    std::cout << "alpha = " << alpha << std::endl;
    std::cout << "alpha_fac = " << alpha_fac << std::endl;
    for (iDim = 0; iDim < nDim; iDim++) {
      for (jDim = 0; jDim < nDim; jDim++) {
        std::cout << "tauSGET[" << iDim << "][" << jDim << "] = "
                  << tauSGET[iDim][jDim] << std::endl;
        std::cout << "Gt[" << iDim << "][" << jDim << "] = "
                  << Gt[iDim][jDim] << std::endl;
        std::cout << "muSGET[" << iDim << "][" << jDim << "] = "
                  << aniso_viscosity[iDim][jDim] << std::endl;
      }
    }
    SU2_MPI::Error("invLengthTensor has NaN!", CURRENT_FUNCTION);
  }

#endif

}


void CHybrid_Mediator::SolveEigen(const su2double* const* M,
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

void CHybrid_Mediator::SolveGeneralizedEigen(const su2double* const* A,
                                             const su2double* const* B,
					     vector<su2double> &eigvalues,
					     vector<vector<su2double> > &eigvectors) {

  eigvalues.resize(nDim);
  eigvectors.resize(nDim, std::vector<su2double>(nDim));

#ifdef HAVE_LAPACK

  unsigned short iDim, jDim;

#ifndef NDEBUG
  const su2double tolerance = 1E-10;
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim< nDim; jDim++) {
      const su2double tolerance = 1E-10;
      if (abs(A[iDim][jDim] - A[jDim][iDim]) > tolerance) {
        SU2_MPI::Error("Matrix A is not symmetric!", CURRENT_FUNCTION);
      }
      if (abs(B[iDim][jDim] - B[jDim][iDim]) > tolerance) {
        SU2_MPI::Error("Matrix A is not symmetric!", CURRENT_FUNCTION);
      }
    }
  }
#endif // NDEBUG

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


void CHybrid_Mediator::SetFluctuatingStress(CFluctuatingStress* fluct_stress) {
  fluct_stress_model = fluct_stress;
}

void CHybrid_Mediator::SetForcingModel(CHybridForcingAbstractBase* forcing) {
  forcing_model = forcing;
}

const su2double* CHybrid_Mediator::GetForcingVector(unsigned long iPoint) {
  return forcing_model->GetForcingVector(iPoint);
}


CHybrid_Dummy_Mediator::CHybrid_Dummy_Mediator(unsigned short nDim,
                                               CConfig* config)
    : nDim(nDim) {

  zero_vector = new su2double[nDim];

  zero_tensor = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    zero_vector[iDim] = 0;
    zero_tensor[iDim] = new su2double[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      zero_tensor[iDim][jDim] = 0;
    }
  }
}

CHybrid_Dummy_Mediator::~CHybrid_Dummy_Mediator() {


  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] zero_tensor[iDim];
  }
  delete [] zero_tensor;
  delete [] zero_vector;
}

void CHybrid_Dummy_Mediator::SetupRANSNumerics(CSolver **solver_container,
                                               CNumerics* rans_numerics,
                                               unsigned long iPoint) {

  rans_numerics->SetKineticEnergyRatio(1.0);
  rans_numerics->SetResolvedTurbStress(zero_tensor);
}

void CHybrid_Dummy_Mediator::ComputeResolutionAdequacy(const CGeometry* geometry,
                                                       CSolver **solver_container,
                                                       unsigned long iPoint) {
  // Set resolution adequacy in the CNSVariables class
  solver_container[FLOW_SOL]->node[iPoint]->SetResolutionAdequacy(1.0);
}

void CHybrid_Dummy_Mediator::ComputeForcingField(CSolver** solver,
                                                 CGeometry *geometry,
                                                 CConfig *config) {
  // nothing to compute
}


void CHybrid_Dummy_Mediator::SetupResolvedFlowSolver(const CGeometry* geometry,
                                               CSolver **solver_container,
                                               unsigned long iPoint) {
  solver_container[FLOW_SOL]->node[iPoint]->SetAnisoEddyViscosity(zero_tensor);
}

void CHybrid_Dummy_Mediator::SetupResolvedFlowNumerics(const CGeometry* geometry,
                                             CSolver **solver_container,
                                             CNumerics* visc_numerics,
                                             unsigned long iPoint,
                                             unsigned long jPoint) {

  /*--- Set the "average" variables to be identical to the resolved
   * variables. ---*/

  su2double* primvar_i =
      solver_container[FLOW_SOL]->node[iPoint]->GetPrimitive();
  su2double* primvar_j =
      solver_container[FLOW_SOL]->node[jPoint]->GetPrimitive();

  su2double** primvar_grad_i =
      solver_container[FLOW_SOL]->node[iPoint]->GetGradient_Primitive();
  su2double** primvar_grad_j =
      solver_container[FLOW_SOL]->node[jPoint]->GetGradient_Primitive();

  /*--- This should be the zero tensor, if SetupResolvedFlowSolver was
   * called properly ---*/

  su2double** aniso_viscosity_i =
      solver_container[FLOW_SOL]->node[iPoint]->GetAnisoEddyViscosity();
  su2double** aniso_viscosity_j =
      solver_container[FLOW_SOL]->node[jPoint]->GetAnisoEddyViscosity();

  CAvgGrad_Hybrid* numerics = dynamic_cast<CAvgGrad_Hybrid*>(visc_numerics);
  numerics->SetPrimitive_Average(primvar_i, primvar_j);
  numerics->SetPrimVarGradient_Average(primvar_grad_i, primvar_grad_j);
  numerics->SetAniso_Eddy_Viscosity(aniso_viscosity_i, aniso_viscosity_j);
  numerics->SetKineticEnergyRatio(1.0, 1.0);
}

void CHybrid_Dummy_Mediator::SetFluctuatingStress(CFluctuatingStress* fluct_stress) {
  /*--- We won't use it, so just delete it. ---*/
  delete fluct_stress;
}

void CHybrid_Dummy_Mediator::SetForcingModel(CHybridForcingAbstractBase* forcing) {
  /*--- We won't use it, so just delete it. ---*/
  delete forcing;
}

const su2double* CHybrid_Dummy_Mediator::GetForcingVector(unsigned long iPoint) {
  // dummy mediator has no forcing
  return zero_vector;
}

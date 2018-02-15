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
  : rans_weight(1.0), CHybrid_Visc_Anisotropy(nDim) {
}

void CHybrid_Aniso_Q::CalculateViscAnisotropy() {
#ifndef NDEBUG
  if (rans_weight < 0 || rans_weight > 1) {
    cout << "ERROR: Isotropic weighting in hybrid RANS/LES was not in range [0,1]" << endl;
    cout << "       weight = " << rans_weight << endl;
    exit(EXIT_FAILURE);
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
    cout << "ERROR: Improper number of scalars passed to anisotropy model!" << endl;
    cout << "       Expected: 1" << endl;
    cout << "       Found:    " << val_scalars.size() << endl;
    exit(EXIT_FAILURE);
  }
#endif
  rans_weight = val_scalars[0];
}



CHybrid_Mediator::CHybrid_Mediator(int nDim, CConfig* config, string filename)
   : nDim(nDim), C_sf(config->Get_Hybrid_Model_Const()), config(config) {

  hybrid_forcing = new CHybridForcing(nDim);

}

CHybrid_Mediator::~CHybrid_Mediator() {

  delete hybrid_forcing;

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

  // FIXME: Update with new r_delta
  su2double r_delta = CalculateRk(NULL, NULL);
  solver_container[HYBRID_SOL]->node[iPoint]->SetResolutionAdequacy(r_delta);

}

void CHybrid_Mediator::SetupHybridParamNumerics(CGeometry* geometry,
                                                CSolver **solver_container,
                                                CNumerics *hybrid_param_numerics,
                                                unsigned short iPoint,
                                                unsigned short jPoint) {

  /*--- Find and store turbulent length and timescales ---*/
  su2double turb_T =
      solver_container[TURB_SOL]->node[iPoint]->GetTurbTimescale();

#ifndef NDEBUG
  if (turb_T <= 0) {
    cout << "Error: Turbulent timescale was " << turb_T << std::endl;
    exit(EXIT_FAILURE);
  }
#endif

  hybrid_param_numerics->SetTurbTimescale(turb_T);

  /*--- Pass resolution adequacy into the numerics object ---*/

  su2double r_k = solver_container[HYBRID_SOL]->node[iPoint]->GetResolutionAdequacy();
  hybrid_param_numerics->SetResolutionAdequacy(r_k);
}

void CHybrid_Mediator::SetupStressAnisotropy(CGeometry* geometry,
                                             CSolver **solver_container,
                                             CHybrid_Visc_Anisotropy* hybrid_anisotropy,
                                             unsigned short iPoint) {
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

su2double CHybrid_Mediator::GetProjResolution(su2double** resolution_tensor,
                                              vector<su2double> direction) {

#ifndef NDEBUG
  su2double magnitude_squared = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    magnitude_squared += direction[iDim]*direction[iDim];
  if (abs(magnitude_squared - 1.0) > 1e-7) {
    cout << "ERROR: The unit vector for the projected resolution calc had a ";
    cout << "magnitude greater than 1!" << endl;
    cout << "    Magnitude: " << sqrt(magnitude_squared) << endl;
    exit(EXIT_FAILURE);
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

void CHybrid_Mediator::SolveEigen(su2double** M,
                                  vector<su2double> &eigvalues,
                                  vector<vector<su2double> > &eigvectors) {
unsigned short iDim, jDim;

#ifndef NDEBUG
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (M[iDim][jDim] != M[iDim][jDim]) {
        cout << "ERROR: SolveEigen received a matrix with NaN!" << endl;
        exit(EXIT_FAILURE);
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
    cout << "ERROR: SolveEigen received an empty matrix!" << endl;
    exit(EXIT_FAILURE);
  }

  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (iDim == jDim) continue;
      if (abs(M[iDim][jDim] - M[jDim][iDim]) > 1e-10 &&
          abs(M[iDim][jDim] - M[jDim][iDim])/abs(M[iDim][jDim]) > 1e-6) {
        cout << "ERROR: SolveEigen received a non-symmetric matrix!" << endl;
        cout << "    The difference between elements" << endl;
        cout << "      [" << iDim << ", " << jDim << "] and [" << jDim << ", " << iDim << "]" << endl;
        cout << "      was: " << M[iDim][jDim] - M[jDim][iDim] << endl;
        cout << "    Matrix:" << endl;
        cout << "      [[" << M[0][0] << ", " << M[0][1] << ", " << M[0][2] << "]" << endl;
        cout << "       [" << M[1][0] << ", " << M[1][1] << ", " << M[1][2] << "]" << endl;
        cout << "       [" << M[2][0] << ", " << M[2][1] << ", " << M[2][2] << "]]" << endl;
        exit(EXIT_FAILURE);
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
    cout << "ERROR: The solver failed to compute eigenvalues." << endl;
    exit(EXIT_FAILURE);
  }

  /*--- Rewrite arrays to eigenvalues output ---*/

  for (iDim = 0; iDim < nDim; iDim++)
    eigvalues[iDim] = eigval[iDim];

  /*--- Check the values ---*/

  for (iDim = 0; iDim < nDim; iDim++) {
    if (eigval[iDim] < 0.0 && eigval[iDim] > -1e-4) {
        eigvalues[iDim] = 0.0;
    } else if (eigval[iDim] < -1e-4) {
      cout << "ERROR: The solver returned a large negative eigenvalue!" << endl;
      cout << "    Eigenvalues: [";
      cout << eigval[0] << ", ";
      cout << eigval[1] << ", ";
      cout << eigval[2] << "]" << endl;
      exit(EXIT_FAILURE);
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
  cout << "Eigensolver without LAPACK not implemented; please use LAPACK." << endl;
  exit(EXIT_FAILURE);
#endif

}

void CHybrid_Mediator::AdjustEddyViscosity(CSolver **solver_container,
                                           unsigned short iPoint,
                                           su2double *eddy_viscosity) {
  su2double* alpha = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
  *eddy_viscosity *= alpha[0]*alpha[0];
}



CHybrid_Dummy_Mediator::CHybrid_Dummy_Mediator(int nDim, CConfig* config) : nDim(nDim) {

  /*--- Set the default value of the hybrid parameter ---*/
  dummy_alpha = new su2double[1];
  dummy_alpha[0] = 1.0;
}

CHybrid_Dummy_Mediator::~CHybrid_Dummy_Mediator() {
  delete [] dummy_alpha;
}

void CHybrid_Dummy_Mediator::SetupRANSNumerics(CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* rans_numerics,
                                         unsigned short iPoint,
                                         unsigned short jPoint) {
  // This next line is just here for testing purposes.
  su2double* alpha = solver_container[HYBRID_SOL]->node[iPoint]->GetSolution();
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

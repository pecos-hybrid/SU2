/*!
 * \file hybrid_RANS_LES_model.hpp
 * \brief 
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

#pragma once

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include "fluctuating_stress.hpp"
#include "numerics_structure.hpp"
#include "hybrid_RANS_LES_forcing.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

#include "stdio.h"
#include "math.h"

// Forward declarations to resolve circular dependencies
class CSolver;
class CHybridForcingAbstractBase;

using namespace std;

/*!
 * \class CAbstract_Hybrid_Mediator
 * \brief Base abstract class for a hybrid model mediator object.
 *
 * In order to decouple the RANS model, the energy transfer model,
 * and the resolved flow, a mediator object is necessary.  This allows the
 * RANS, energy transfer model, and resolved flow equations to follow the
 * single responsibility principle, while this class makes sure they
 * have the information they need.
 *
 * The main purpose of this class is to pass variables to where they
 * need to go.
 *
 * \author: C. Pederson
 */
class CAbstract_Hybrid_Mediator {
 public:
  /*!
   * \brief A virtual destructor.
   */
  virtual ~CAbstract_Hybrid_Mediator() {};

  /**
   * \brief Retrieve and pass along all necessary info for the RANS model.
   *
   * \param[in] solver_container - An array of solvers
   * \param[in] rans_numerics - The source numerics for the turb. solver
   * \param[in] iPoint - The number of the node being evaluated
   */
  virtual void SetupRANSNumerics(CSolver **solver_container,
                                 CNumerics* rans_numerics,
                                 unsigned long iPoint) = 0;

  /**
   * \brief Compute instantaneous resolution adequacy and store with variables
   *
   * \param[in] geometry - The grid information
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The number of the node being evaluated
   */
  virtual void ComputeResolutionAdequacy(const CGeometry* geometry,
                                         CSolver **solver_container,
                                         unsigned long iPoint) = 0;

  /**
   * \brief Evaluate forcing field
   */
  virtual void ComputeForcingField(CSolver** solver, CGeometry *geometry,
                                   CConfig *config) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for the resolved flow.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The number of the node being evaluated
   */
  virtual void SetupResolvedFlowSolver(const CGeometry* geometry,
                                       CSolver **solver_container,
                                       unsigned long iPoint) = 0;

  /**
   * \brief Retrieve and pass along all necessary info for resolved numerics.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] visc_numerics - The viscous numerics for the resolved flow solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  virtual void SetupResolvedFlowNumerics(const CGeometry* geometry,
                                         CSolver **solver_container,
                                         CNumerics* visc_numerics,
                                         unsigned long iPoint,
                                         unsigned long jPoint) = 0;

  /*!
   * \brief Set the fluctuating stress model to be used.
   *
   * This will take ownership of the fluctuating stress model.
   *
   * \param[in] fluct_stress - A model for the fluctuating stress
   */
  virtual void SetFluctuatingStress(CFluctuatingStress* fluct_stress) = 0;

  virtual void SetForcingModel(CHybridForcingAbstractBase* forcing_model) = 0;
  virtual const su2double* GetForcingVector(unsigned long iPoint) = 0;
};


/*!
 * \class CHybrid_Mediator
 * \brief Mediator object for the model-split hybrid RANS/LES model.
 * \author: C. Pederson
 */
class CHybrid_Mediator : public CAbstract_Hybrid_Mediator {
 protected:

  unsigned short nDim;
  su2double C_zeta; /*!> \brief Scaling constant for the transformation tensor zeta */
  su2double **Q,        /*!> \brief An approximate 2nd order structure function tensor */
            **Qapprox;  /*!> \brief An approximate 2nd order structure function tensor (used for temporary calculations) */
  su2double **invLengthTensor; /*!> \brief Inverse length scale tensor formed from production and v2 (or tke, depending on availability) */
  su2double **aniso_eddy_viscosity; /*!> \brief A 2D array used to hold the value of the anisotropic eddy viscosity during calculations. */
  std::vector<std::vector<su2double> > constants;
  CFluctuatingStress* fluct_stress_model;
  CConfig* config;
  CHybridForcingAbstractBase* forcing_model;


  /*--- Data structures for LAPACK ---*/
#ifdef HAVE_LAPACK
  int info, lwork;
  su2double wkopt;
  su2double* work;
  su2double eigval[3], eigvec[9], vr[9], wi[3];
  su2double mat[9], matb[9];
  int num_found;
  int isupp[3];
#endif

  /*!
   * \brief Calculates the resolution inadequacy parameter
   * \param[in] Q - The approximate 2nd order structure function
   * \param[in] v2 - The v2 value from Durbin's k-eps-v2-f model
   * \return The resolution inadequacy parameter
   */
  su2double CalculateRk(const su2double* const* Q, su2double v2);

  /*!
   * \brief Projects the resolution on a specific vector
   * \param[in] resolution_tensor - The tensor representing separation distances
   * \param[in] direction - The direction vector (assumed to be normalized)
   * \return The magnitude of the resolution tensor projected on the direction.
   */
  su2double GetProjResolution(const su2double* const *resolution_tensor,
                              const vector<su2double>& direction);

  /**
   * \brief Uses a resolution tensor and a gradient-gradient tensor to build an
   *        approximate two-point structure function tensor
   * \param[in] val_ResolutionTensor - A tensor representing cell-cell distances
   * \param[in] val_PrimVar_Grad - The gradient in the resolved velocity field.
   * \param[out] val_Q - An approximate resolution-scale two-point second-order
   *                 structure function.
   */
  template <class T>
  void CalculateApproxStructFunc(T val_ResolutionTensor,
                                 const su2double* const* val_PrimVar_Grad,
                                 su2double** val_Q);

  /*!
   * \brief Loads the model fit constants from *.dat files
   * \param[in] filename - The base name for the files (e.g. [filename]0.dat)
   * \return The 3 sets of constants pulled from the files.
   */
  vector<vector<su2double> > LoadConstants(const string& filename);

  /*!
   * \brief Solve for the eigenvalues of Q, given eigenvalues of M
   *
   * Solves for the eigenvalues of the expected value of the contracted velocity
   * differences at grid resolution (Q), given the eigenvalues of the
   * resolution tensor.
   *
   * \param[in] eig_values_M - Eigenvalues of the resolution tensor.
   * \return The eigenvalues of the expected value of the approx SF tensor.
   */
  vector<su2double> GetEigValues_Q(const vector<su2double>& eig_values_M);

  /*!
   * \brief Calculates the eigenvalues of a modified resolution tensor.
   * \param eig_values_M - Eigenvalues of the grid-based resolution tensor.
   * \return Eigenvalues of a transformation mapping the approximate velocity
   *         differences at grid resolution to the two-point second-order
   *         structure function.
   */
  vector<su2double> GetEigValues_Zeta(const vector<su2double>& eig_values_M);

 public:
  /*!
   * \brief Builds a transformation for the approximate structure function.
   * \param[in] values_M - The cell-to-cell distances in the "principal
   *            "directions"
   */
  vector<vector<su2double> > BuildZeta(const su2double* values_M,
                                       const su2double* const* vectors_M);


  /**
   * \brief Constructor for the hybrid mediator object.
   * \param[in] nDim - The number of dimensions of the problem
   * \param[in] CConfig - The configuration for the current zone
   */
  CHybrid_Mediator(unsigned short nDim, CConfig* config, const string& filename="");

  /**
   * \brief Destructor for the hybrid mediator object.
   */
  ~CHybrid_Mediator();

  /*!
   * \brief Calculates the production-based inverse length scale tensor
   * \param[in] flow_vars - Pointer to mean flow variables
   * \param[in] turb_vars - Pointer to turbulence model variables
   * \param[in] hybr_vars - Pointer to hybrid model variables
   */
  void ComputeInvLengthTensor(CVariable* flow_vars,
                              CVariable* flow_avgs,
                              CVariable* turb_vars,
                              su2double val_alpha,
                              int short hybrid_res_ind);

  su2double GetInvLengthScale(unsigned short ival, unsigned short jval) {
    return invLengthTensor[ival][jval];
  }


  /**
   * \brief Retrieve and pass along all necessary info for the RANS model.
   *
   * \param[in] solver_container - An array of solvers
   * \param[in] rans_numerics - The source numerics for the turb. solver
   * \param[in] iPoint - The number of the node being evaluated
   */
  void SetupRANSNumerics(CSolver **solver_container,
                         CNumerics* rans_numerics,
                         unsigned long iPoint);

  void ComputeResolutionAdequacy(const CGeometry* geometry,
                                 CSolver **solver_container,
                                 unsigned long iPoint);

  void ComputeForcingField(CSolver** solver, CGeometry *geometry,
                           CConfig *config);

  /**
   * \brief Retrieve and pass along all necessary info to calculate
   *        variables needed for the flow solver.
   *
   * \param[in] geometry - The grid information
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The number of the node being evaluated
   */
  void SetupResolvedFlowSolver(const CGeometry* geometry,
                               CSolver **solver_container,
                               unsigned long iPoint);

  /**
   * \brief Retrieve and pass along all necessary info for resolved numerics.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] visc_numerics - The viscous numerics for the resolved flow solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupResolvedFlowNumerics(const CGeometry* geometry,
                             CSolver **solver_container,
                             CNumerics* visc_numerics,
                             unsigned long iPoint,
                             unsigned long jPoint);

  /**
   * \brief Returns the constants for the numerical fit for the resolution tensor.
   * \return Constants for the numerical fit for the resolution tensor.
   */
  vector<vector<su2double> > GetConstants();


  void SolveEigen(const su2double* const* M, vector<su2double> &eigvalues,
                  vector<vector<su2double> > &eigvectors);

  void SolveGeneralizedEigen(const su2double* const* A, const su2double* const* B,
			     vector<su2double> &eigvalues,
			     vector<vector<su2double> > &eigvectors);

  /*!
   * \brief Set the fluctuating stress model to be used.
   *
   * This will take ownership of the fluctuating stress model.
   *
   * \param[in] fluct_stress - A model for the fluctuating stress
   */
  void SetFluctuatingStress(CFluctuatingStress* fluct_stress);

  void SetForcingModel(CHybridForcingAbstractBase* forcing);
  const su2double* GetForcingVector(unsigned long iPoint);
};

/*!
 * \class CHybrid_Dummy_Mediator
 * \brief Mediator object for RANS-only operation; isotropic stress and dummy hybrid parameter
 * \author: C. Pederson
 */
class CHybrid_Dummy_Mediator : public CAbstract_Hybrid_Mediator {
 protected:

  unsigned short nDim;
  su2double*  zero_vector; /*!< \brief A zero nDim x nDim tensor */
  su2double** zero_tensor; /*!< \brief A zero nDim x nDim tensor */


 public:

  /**
   * \brief Constructor for the hybrid mediator object.
   * \param[in] nDim - The number of dimensions of the problem
   * \param[in] CConfig - The configuration for the current zone
   */
  CHybrid_Dummy_Mediator(unsigned short nDim, CConfig* config);

  /**
   * \brief Destructor for the hybrid mediator object.
   */
  ~CHybrid_Dummy_Mediator();

  /**
   * \brief Retrieve and pass along all necessary info for the RANS model.
   *
   * \param[in] solver_container - An array of solvers
   * \param[in] rans_numerics - The source numerics for the turb. solver
   * \param[in] iPoint - The number of the node being evaluated
   */
  void SetupRANSNumerics(CSolver **solver_container,
                         CNumerics* rans_numerics,
                         unsigned long iPoint);

  void ComputeResolutionAdequacy(const CGeometry* geometry,
                                 CSolver **solver_container,
                                 unsigned long iPoint);

  void ComputeForcingField(CSolver** solver, CGeometry *geometry,
                           CConfig *config);

  /**
   * \brief Retrieve and pass along all necessary info to calculate
   *        variables needed for the flow solver.
   *
   * \param[in] geometry - The grid information
   * \param[in] solver_container - An array of solvers
   * \param[in] iPoint - The number of the node being evaluated
   */
  void SetupResolvedFlowSolver(const CGeometry* geometry,
                               CSolver **solver_container,
                               unsigned long iPoint);

  /**
   * \brief Retrieve and pass along all necessary info for resolved numerics.
   *
   * \param[in] geometry - A pointer to the geometry
   * \param[in] solver_container - An array of solvers
   * \param[in] visc_numerics - The viscous numerics for the resolved flow solver
   * \param[in] iPoint - The number of the node being evaluated
   * \param[in] jPoint - The number of the opposite node
   */
  void SetupResolvedFlowNumerics(const CGeometry* geometry,
                             CSolver **solver_container,
                             CNumerics* visc_numerics,
                             unsigned long iPoint,
                             unsigned long jPoint);

  /*!
   * \brief Set the fluctuating stress model to be used.
   *
   * This will take ownership of the fluctuating stress model.
   *
   * \param[in] fluct_stress - A model for the fluctuating stress
   */
  void SetFluctuatingStress(CFluctuatingStress* fluct_stress);

  void SetForcingModel(CHybridForcingAbstractBase* forcing);
  const su2double* GetForcingVector(unsigned long iPoint);
};

/*--- Template definitions:
 * These must be placed with the template declarations.  They can be kept here
 * or moved to an *.inl file that is included in this header. ---*/

template <class T>
void CHybrid_Mediator::CalculateApproxStructFunc(T val_ResolutionTensor,
                                                 const su2double* const* val_PrimVar_Grad,
                                                 su2double** val_Q) {
  unsigned int iDim, jDim, kDim, lDim, mDim;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      val_Q[iDim][jDim] = 0.0;

  for (iDim = 0; iDim < nDim; iDim++)
    for (jDim = 0; jDim < nDim; jDim++)
      for (kDim = 0; kDim < nDim; kDim++)
        for (lDim = 0; lDim < nDim; lDim++)
          for (mDim = 0; mDim < nDim; mDim++)
            val_Q[iDim][jDim] += val_ResolutionTensor[iDim][mDim]*
                             val_PrimVar_Grad[kDim+1][mDim]*
                             val_PrimVar_Grad[kDim+1][lDim]*
                             val_ResolutionTensor[lDim][jDim];
}


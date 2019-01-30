/*!
 * \file numerics_structure.hpp
 * \brief Headers of the main subroutines for the dumerical definition of the problem.
 *        The subroutines and functions are in the <i>numerics_structure.cpp</i>,
 *        <i>numerics_convective.cpp</i>, <i>numerics_viscous.cpp</i>, and
 *        <i>numerics_source.cpp</i> files.
 * \author F. Palacios, T. Economon
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

#pragma once

#include "numerics_structure.hpp"

/*!
 * \class CUpwSca_TurbKE
 * \brief Upwind convective flux for the zeta-f KE turbulence model equations.
 * \ingroup ConvDiscr
 * \author S. Haering.
 * \version 5.0.0 "Raven"
 */
class CUpwSca_TurbKE : public CUpwScalar {
private:

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn();

  /*!
   * \brief KE specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                          su2double **Jacobian_j, CConfig *config);

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CUpwSca_TurbKE(unsigned short val_nDim,
                 unsigned short val_nVar,
                 CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CUpwSca_TurbKE(void);

};

/*!
 * \class CAvgGrad_TurbKE
 * \brief Computes viscous term using average of gradient (zeta-f KE model).
 * \ingroup ViscDiscr
 * \author S. Haering
 */
class CAvgGrad_TurbKE : public CAvgGrad_Scalar {
private:
  su2double sigma_k,                     /*!< \brief Constants for the viscous terms, k-w (1), k-eps (2)*/
  sigma_e,
  sigma_z;

  su2double diff_kine,                     /*!< \brief Diffusivity for viscous terms of tke eq */
    diff_epsi,                           /*!< \brief Diffusivity for viscous terms of epsi eq */
    diff_zeta,                           /*!< \brief Diffusivity for viscous terms of zeta eq */
    diff_f;                           /*!< \brief Diffusivity for viscous terms of f eq */

  /*!
   * \brief Adds any extra variables to AD
   */
  void ExtraADPreaccIn(void);

  /*!
   * \brief KE specific steps in the ComputeResidual method
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void FinishResidualCalc(su2double *val_residual, su2double **Jacobian_i,
                          su2double **Jacobian_j, CConfig *config);


public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] constants - KE model parameters
   * \param[in] correct_grad - Use corrected grad variant
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_TurbKE(unsigned short val_nDim, unsigned short val_nVar,
                  su2double* constants, bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_TurbKE(void);

};

/*!
 * \class CSourcePieceWise_TurbKE
 * \brief Compute source terms of the zeta-f KE turbulence model equations.
 * \ingroup SourceDiscr
 * \author S. Haering.
 */
class CSourcePieceWise_TurbKE : public CNumerics {
private:
  su2double Lm, ///< The turbulent lengthscale
            Tm; ///< The turbulent timescale

  su2double sigma_k,
  sigma_e,
  sigma_z,
  C_L,
  C_T,
  C_1,
  C_2,
  C_mu,
  C_2p,
  C_eta,
  C_e1o,
  C_e2;

  su2double** ResolvedTurbStress;
  su2double SGSProduction;
  su2double Pk;

  bool incompressible;

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CSourcePieceWise_TurbKE(unsigned short val_nDim, unsigned short val_nVar,
                          su2double* constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CSourcePieceWise_TurbKE(void);

  /*!
   * \brief Residual for source term integration.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian wrt node i soln (for implicit).
   * \param[out] val_Jacobian_j - Jacobian wrt node j soln (for implicit).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual,
                       su2double **val_Jacobian_i, su2double **val_Jacobian_j,
                       CConfig *config);

  /*!
   * \brief Set the resolved turbulent stress
   * \param[in] val_turb_stress - The turbulent stress to be used
   */
  void SetResolvedTurbStress(su2double** val_turb_stress) {
    ResolvedTurbStress = val_turb_stress;
  }

  void SetProduction(su2double val_production) {
    Pk = val_production;
  }

  /*!
   * \brief Get the SGS production.
   * \return The SGS production term.
   */
  su2double GetSGSProduction(void) const { return SGSProduction; }

  su2double GetProduction(void) const { return Pk; }

};

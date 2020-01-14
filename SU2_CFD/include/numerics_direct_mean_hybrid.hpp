/*!
 * \file numerics_structure.hpp
 * \brief Header for the viscous numerics in the model-split hybrid
 *        RANS/LES model.
 * \author C. Pederson
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
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author C. Pederson
 */
class CAvgGrad_Hybrid : public CAvgGrad_Base {
private:
  /*--- Temporary variables used in calculation ---*/
  su2double **conductivity;  /*!< \brief Anisotropic thermal conductivity */
  su2double **deviatoric;    /*!< \brief Deviatoric part of the velocity gradient tensor. */
  su2double* Mean_PrimVar_Average; /*!< \brief The primitive variables computed from average flow variables */
  su2double** Mean_GradPrimVar_Average; /*!< \brief The gradient of the primitve variables from the average flow */ 
  su2double** Mean_GradPrimVar_Fluct; /*!< \brief The gradient of the fluctuating primitive variables */
  su2double** Mean_Aniso_Eddy_Viscosity;  /*!< \brief The anisotropic eddy viscosity from the energy transfer model. */
  /*--- Values to be set externally for residual calculation ---*/
  su2double* PrimVar_Average_i; /*!< \brief The primitive variables, computed from average flow variables */
  su2double* PrimVar_Average_j; /*!< \brief The primitive variables, computed from average flow variables */
  su2double **Aniso_Eddy_Viscosity_i, **Aniso_Eddy_Viscosity_j; /*!< \brief The anisotropic eddy viscosity from the energy transfer model. */
  su2double **PrimVar_Grad_Average_i, **PrimVar_Grad_Average_j; /*!< \brief The gradient of the primitve variables from the average flow */
  su2double alpha_i, alpha_j; /*!< \brief The kinetic energy ratio (of modeled to total turbulent kinetic energy) */

public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] val_correct_grad - Apply a correction to the gradient
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Hybrid(unsigned short val_nDim, unsigned short val_nVar, bool val_correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Hybrid(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);


  /*!
   * \brief Compute the Jacobian of the heat flux vector
   *
   * This Jacobian is projected onto the normal vector, so it is of
   * dimension nVar.
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_gradprimvar - Mean value of the gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           const su2double val_laminar_viscosity,
                           const su2double val_eddy_viscosity,
                           const su2double val_dist_ij,
                           const su2double *val_normal);

  /*!
   * \brief Compute the heat flux due to only molecular viscosity.
   * \param[in] val_gradprimvar - The resolved gradient of primitive variables
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   */
  void SetLaminarHeatFlux(su2double **val_gradprimvar,
                          const su2double val_laminar_viscosity);

  /*!
   * \brief Add the SGS heat flux to a previously computed heat flux vector
   *
   * SetLaminarHeatFlux must be called before this in order to initialize
   * the heat flux vector.
   *
   * \param[in] val_gradprimvar - The gradient of the average primitive variables
   * \param[in] val_alpha - The ratio of modeled to total turbulent kinetic energy
   * \param[in] val_eddy_viscosity - The average (i.e. RANS) eddy viscosity
   */
  void AddSGSHeatFlux(su2double **val_gradprimvar,
                      const su2double val_alpha,
                      const su2double val_eddy_viscosity);

  /*!
   * \brief Add the heat flux from the energy transfer model to the heat flux
   *
   * SetLaminarHeatFlux must be called before this in order to initialize
   * the heat flux vector.
   *
   * \param[in] val_gradprimvar - The gradient of the fluctuating primitive variables
   * \param[in] val_eddy_viscosity - An anisotropic, flucutating eddy viscosity
   */
  void AddSGETHeatFlux(su2double** val_gradprimvar,
                       su2double** val_eddy_viscosity);

  /*!
   * \brief Set diffusion of total energy due to unresolved molecular
   * diffusion.
   *
   * \param[in] val_gradturbvar - Mean value of the gradient of the
   *     turbulence variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   */
  void SetLaminar_TKE_Diffusion(const su2double* const *val_gradturbvar,
                                su2double val_laminar_viscosity);

  /*!
   * \brief Set diffusion of total energy due to turbulent transport.
   *
   * \param[in] val_gradturbvar - Mean value of the gradient of the
   *     turbulence variables.
   * \param[in] val_alpha - Ratio of turbulent kinetic energy
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   */
  void AddSGS_TKE_Diffusion(const su2double* const *val_gradturbvar,
                            su2double val_alpha,
                            su2double val_eddy_viscosity);

  /*!
   * \brief Compute the viscous stress tensor due to only the molecular viscosity.
   * \param[in] val_gradprimvar - The resolved gradient of primitive variables
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   */
  void SetLaminarStressTensor(su2double **val_gradprimvar,
                              const su2double val_laminar_viscosity);

  /*!
   * \brief Compute the mean subfilter stress.
   *
   * SetLaminarSressTensor must be called before this function in order
   * to initialize the stress tensor.
   *
   * \param[in] val_primvar - The average primitive variables
   * \param[in] val_gradprimvar - The gradient of the average primitive variables
   * \param[in] val_alpha - The ratio of modeled to total turbulent kinetic energy
   * \param[in] val_turb_ke - The total (not subfilter) turbulent kinetic energy
   * \param[in] val_eddy_viscosity - The average (i.e. RANS) eddy viscosity
   */
  void AddTauSGS(const su2double *val_primvar,
                 su2double **val_gradprimvar,
                 const su2double val_alpha,
                 const su2double val_turb_ke,
                 const su2double val_eddy_viscosity);

  /*!
   * \brief Add the stress tensor from the energy transfer model to the
   *        viscous stress
   * \param[in] val_gradprimvar - The gradient of the fluctuating primitive variables
   * \param[in] val_eddy_viscosity - An anisotropic, flucutating eddy viscosity
   */
  void AddTauSGET(su2double **val_gradprimvar,
                  su2double **val_eddy_viscosity);

  void AddTauSGETJacobian(su2double *val_Mean_PrimVar,
			  su2double **mu,
			  const su2double val_dist_ij,
			  const su2double *nvec);

  void AddSGETHeatFluxJacobian(su2double *val_Mean_PrimVar,
			       su2double **mu,
			       const su2double val_dist_ij,
			       const su2double *val_normal);



  /*!
   * \brief Set the anisotropic eddy viscosity from the energy transfer model.
   * \param[in] aniso_eddy_viscosity_i - The eddy viscosity at point i
   * \param[in] aniso_eddy_viscosity_j - The eddy viscosity at point j
   */
  void SetAniso_Eddy_Viscosity(su2double** aniso_eddy_viscosity_i,
                               su2double** aniso_eddy_viscosity_j);

  /*!
   * \brief Set the average primitive variables.
   *
   * This sets the primitive variables from the average flow, which is not
   * the same as the face values (computed using an average) of the resolved
   * flow.
   *
   * \param[in] val_primvar_average_i - The average primitive variables at point i 
   * \param[in] val_primvar_average_j - The average primitive variables at point j 
   */
  void SetPrimitive_Average(su2double* val_primvar_average_i,
                            su2double* val_primvar_average_j);

  /*!
   * \brief Set the gradient of the primitive variables for the average.
   *
   * This sets the gradient of the primitive variables from the average
   * flow, which is not the same as the face values (computed using an
   * average) of the resolved flow.
   *
   * \param[in] val_primvar_grad_i - Gradient of the primitive variable at point i.
   * \param[in] val_primvar_grad_j - Gradient of the primitive variable at point j.
   */
  void SetPrimVarGradient_Average(su2double **val_primvar_grad_i,
                                  su2double **val_primvar_grad_j);

  /*!
   * \brief Set the ratio of modeled to total turbulent kinetic energy
   * \param[in] val_alpha_i - The kinetic energy ratio at point i
   * \param[in] val_alpha_j - The kinetic energy ratio at point j
   */
  void SetKineticEnergyRatio(const su2double val_alpha_i,
                             const su2double val_alpha_j);
};

#include "numerics_direct_mean_hybrid.inl"

/*!
 * \file numerics_direct_mean.hpp
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
 * \class CAvgGrad_Base
 * \brief Base class for viscous numerics
 * \ingroup ViscDiscr
 * \author C. Pederson, A. Bueno, and F. Palacios
 */
class CAvgGrad_Base : public CNumerics {
protected:
  const su2double Prandtl_Lam;        /*!< \brief Laminar Prandtl's number. */
  const su2double Prandtl_Turb;    /*!< \brief Turbulent Prandtl's number. */
  su2double *Mean_PrimVar,          /*!< \brief Mean primitive variables. */
  *PrimVar_i, *PrimVar_j,        /*!< \brief Primitives variables at point i and 1. */
  *Edge_Vector,                  /*!< \brief Vector form point i to point j. */
  **Mean_GradPrimVar, *Proj_Mean_GradPrimVar_Edge,  /*!< \brief Mean value of the gradient. */
  Mean_Laminar_Viscosity,                /*!< \brief Mean value of the viscosity. */
  Mean_Eddy_Viscosity,                   /*!< \brief Mean value of the eddy viscosity. */
  Mean_turb_ke,        /*!< \brief Mean value of the turbulent kinetic energy. */
  dist_ij;            /*!< \brief Length of the edge and face. */
  bool implicit; /*!< \brief Implicit calculus. */
  const bool correct_gradient;
  su2double
  **tau, /*!< \brief Viscous + turbulent stress tensor. */
  *viscous_heat_flux,  /*!< \brief Flux of total energy due to molecular and turbulent diffusion */
  **tau_jacobian_i, /*!< \brief Jacobian of the viscous + turbulent stress tensor, projected onto the normal vector. */
  *heat_flux_jac_i; /*!< \brief Jacobian of the molecular + turbulent heat flux vector, projected onto the normal vector. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Base(const unsigned short val_nDim,
                const unsigned short val_nVar,
                const bool correct_grad, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Base(void);
  
  /*!
   * \brief Calculate the viscous and turbulent stress tensor
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_turb_ke - Turbulent kinetic energy
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_qcr
   */
  void GetTau(const su2double *val_primvar, su2double **val_gradprimvar,
              const su2double val_turb_ke,
              const su2double val_laminar_viscosity,
              const su2double val_eddy_viscosity, const bool val_qcr);

  /**
   * \brief Calculate the jacobian of the viscous and turbulent stress tensor
   *
   * This Jacobian is projected onto the normal vector, so it is of dimension
   * [nDim][nVar]
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void GetTauJacobian(const su2double *val_Mean_PrimVar,
                      const su2double val_laminar_viscosity,
                      const su2double val_eddy_viscosity,
                      const su2double val_dist_ij,
                      const su2double *val_normal);

  /*!
   * \brief Compute the projection of the viscous fluxes into a direction.
   * \param[in] val_primvar - Primitive variables.
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void GetViscousProjFlux(const su2double *val_primvar, su2double **val_gradprimvar,
                          const su2double *val_normal);


  /*!
   * \brief TSL-Approximation of Viscous NS Jacobians.
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_dS - Area of the face between two nodes.
   * \param[in] val_Proj_Visc_Flux - Pointer to the projected viscous flux.
   * \param[out] val_Proj_Jac_Tensor_i - Pointer to the projected viscous Jacobian at point i.
   * \param[out] val_Proj_Jac_Tensor_j - Pointer to the projected viscous Jacobian at point j.
   */
  void GetViscousProjJacs(const su2double *val_Mean_PrimVar,
                          const su2double val_dS,
                          const su2double *val_Proj_Visc_Flux,
                          su2double **val_Proj_Jac_Tensor_i,
                          su2double **val_Proj_Jac_Tensor_j);
};

/*!
 * \class CAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author A. Bueno, and F. Palacios
 */
class CAvgGrad_Flow : public CAvgGrad_Base {
public:

  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Apply a correction to the gradients
   * \param[in] config - Definition of the particular problem.
   */
  CAvgGrad_Flow(const unsigned short val_nDim, const unsigned short val_nVar,
                const bool correct_grad, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CAvgGrad_Flow(void);

  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

  /*!
   * \brief Compute the heat flux due to molecular and turbulent diffusivity
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   */
  void GetViscousHeatFlux(su2double **val_gradprimvar,
                          const su2double val_laminar_viscosity,
                          const su2double val_eddy_viscosity);

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
  void GetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           const su2double val_laminar_viscosity,
                           const su2double val_eddy_viscosity,
                           const su2double val_dist_ij,
                           const su2double *val_normal);
};

/*!
 * \class CGeneralAvgGrad_Flow
 * \brief Class for computing viscous term using the average of gradients.
 * \ingroup ViscDiscr
 * \author M.Pini, S. Vitale
 * \version 3.2.1 "eagle"
 */

class CGeneralAvgGrad_Flow : public CAvgGrad_Base {
private:
  su2double *Mean_SecVar,                  /*!< \brief Mean primitive variables. */
  Mean_Thermal_Conductivity,             /*!< \brief Mean value of the thermal conductivity. */
  Mean_Cp;                               /*!< \brief Mean value of the Cp. */
  
public:
  
  /*!
   * \brief Constructor of the class.
   * \param[in] val_nDim - Number of dimension of the problem.
   * \param[in] val_nVar - Number of variables of the problem.
   * \param[in] correct_grad - Apply a correction to the gradients
   * \param[in] config - Definition of the particular problem.
   */
  CGeneralAvgGrad_Flow(const unsigned short val_nDim,
                       const unsigned short val_nVar,
                       const bool correct_grad, CConfig *config);
  
  /*!
   * \brief Destructor of the class.
   */
  ~CGeneralAvgGrad_Flow(void);
  
  /*!
   * \brief Compute the viscous flow residual using an average of gradients.
   * \param[out] val_residual - Pointer to the total residual.
   * \param[out] val_Jacobian_i - Jacobian of the numerical method at node i (implicit computation).
   * \param[out] val_Jacobian_j - Jacobian of the numerical method at node j (implicit computation).
   * \param[in] config - Definition of the particular problem.
   */
  void ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config);

  /*!
   * \brief Compute the heat flux due to molecular and turbulent diffusivity
   * \param[in] val_gradprimvar - Gradient of the primitive variables.
   * \param[in] val_laminar_viscosity - Laminar viscosity.
   * \param[in] val_eddy_viscosity - Eddy viscosity.
   * \param[in] val_thermal_conductivity - Thermal Conductivity.
   * \param[in] val_heat_capacity_cp - Heat Capacity at constant pressure.
   */
  void GetViscousHeatFlux(su2double **val_gradprimvar,
                          const su2double val_laminar_viscosity,
                          const su2double val_eddy_viscosity,
                          const su2double val_thermal_conductivity,
                          const su2double val_heat_capacity_cp);

  /*!
   * \brief Compute the Jacobian of the heat flux vector
   *
   * This Jacobian is projected onto the normal vector, so it is of
   * dimension nVar.
   *
   * \param[in] val_Mean_PrimVar - Mean value of the primitive variables.
   * \param[in] val_gradprimvar - Mean value of the gradient of the primitive variables.
   * \param[in] val_Mean_SecVar - Mean value of the secondary variables.
   * \param[in] val_laminar_viscosity - Value of the laminar viscosity.
   * \param[in] val_eddy_viscosity - Value of the eddy viscosity.
   * \param[in] val_thermal_conductivity - Value of the thermal conductivity.
   * \param[in] val_heat_capacity_cp - Value of the specific heat at constant pressure.
   * \param[in] val_dist_ij - Distance between the points.
   * \param[in] val_normal - Normal vector, the norm of the vector is the area of the face.
   */
  void GetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                           su2double **val_gradprimvar,
                           const su2double *val_Mean_SecVar,
                           const su2double val_laminar_viscosity,
                           const su2double val_eddy_viscosity,
                           const su2double val_thermal_conductivity,
                           const su2double val_heat_capacity_cp,
                           const su2double val_dist_ij,
                           const su2double *val_normal);
};

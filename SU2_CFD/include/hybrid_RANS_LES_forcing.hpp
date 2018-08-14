/*!
 * \file hybrid_RANS_LES_forcing.hpp
 * \brief Interface for the hybrid RANS/LES forcing terms
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

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"
#include "../include/solver_structure.hpp"

/*!
 * \class CHybridForcing
 * \brief Class for defining the forcing needed in a hybrid RANS/LES model
 * \author C. Pederson
 * \version 5.0.0 "Raven"
 *
 * This class could inherit from CNumerics or CSolver, but there's a
 * **lot** of extra code in numerics and solvers that forcing doesn't need.
 * The choice is to have a simpler class with a clearer interface at the
 * cost of some code duplication.
 *
 * Ideally, both would inherit from the same abstract interfaces,
 * but that would require fundamental restructuring of the code.
 */

class CHybridForcing {
 public:
  CHybridForcing(const unsigned short nDim, const unsigned long nPoint,
                 const unsigned long nPointDomain);
  CHybridForcing(CGeometry* geometry, CConfig* config);
  ~CHybridForcing();

  /**
   * \brief Transform physical coordinates into coordinates used for forcing.
   *
   * This method is only left public for testing purposes. It is usually
   * not necessary to call it manually.  See CalculateForcing instead.
   *
   * \param[in] x - Physical coordinates
   * \param[in] iDim - The component to be transformed
   * \param[in] val_time - The current time
   * \return A component of the shifted coordinates.
   */
  su2double TransformCoords(const su2double x, const su2double mean_velocity,
                            const su2double time, const su2double timescale);
  void SetTGField(const su2double* x, const su2double* L, su2double* b);
  void SetStreamFunc(const su2double* x, const su2double* L, su2double* h);
  void SetForcing_Gradient_LS(CGeometry *geometry, CConfig *config);
  void Set_MPI_Forcing_Gradient(CGeometry *geometry, CConfig *config);
  su2double GetTargetProduction(const su2double k_sgs,
                                const su2double dissipation,
                                const su2double resolution_adequacy,
                                const su2double alpha,
                                const su2double laminar_viscosity);
  su2double ComputeScalingFactor(const su2double L, const su2double P_F,
                                 const su2double dt, const su2double* b);
  void ComputeForcingField(CSolver** solver, CGeometry *geometry,
                           CConfig *config);

 private:
  const unsigned short nDim,   /*!< \brief Number of dimensions of the problem. */
  nVar,                        /*!< \brief Number of variables of the problem. */
  nVarGrad;                    /*!< \brief Number of variables for deallocating the LS Cvector. */
  const unsigned long nPoint,  /*!< \brief Number of points of the computational grid. */
  nPointDomain;                /*!< \brief Number of points of the computational grid. */

  su2double** node;
  su2double*** Gradient;       /*!< \brief The indexing is Gradient[iPoint][iVar][iDim] */

  su2double **Smatrix,  /*!< \brief Auxiliary structure for computing gradients by least-squares */
  **Cvector;            /*!< \brief Auxiliary structure for computing gradients by least-squares */

  int*** LeviCivita;    /*!< \brief The alternating or permutation tensor. */
};

#include "hybrid_RANS_LES_forcing.inl"

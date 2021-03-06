/*!
 * \file fluctuating_stress.hpp
 * \brief Interface for the fluctuating stress for hybrid RANS/LES
 * \author C. Pederson
 * \version 6.2.0 "Falcon"
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

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"

/*!
 * \class CFluctuatingStress
 * \brief An abstract class, used to model the energy transfer from
 *        resolved to modeled scales.
 * \author C. Pederson
 */
class CFluctuatingStress {
 public:
  /*!
   * \brief Basic constructor
   * \param[in] val_nDim - Number of dimensions of the problem.
   */
  CFluctuatingStress(unsigned short val_nDim);

  /*!
   * \brief Virtual desctructor.
   */
  virtual ~CFluctuatingStress();

  /**
   * \brief Set the flow primitive variables
   * \param[in] val_primitive - An array containing primitive variables
   *     from the flow
   */
  void SetPrimitive(const su2double* val_primitive);

  /**
   * \brief Set the mean of turbulence variables
   * \param[in] turb_vars - An array containing the mean turbulence variables
   */
  void SetTurbVar(const su2double* val_turb_var);

  /**
   * \brief Calculate the eddy viscosity.
   *
   * This function is pure virtual and must be overloaded in derived classes.
   * The eddy viscosity computed is anisotropic and is a rank-2 tensor.
   *
   * \param[in] geometry - The geometry
   * \param[in] config - The config settings
   * \param[in] iPoint - The current point to be used
   * \param[in] mean_eddy_visc - The mean (i.e. RANS) eddy viscosity.
   * \param[out] eddy_viscosity - The computed anisotropic eddy viscosity
   */
  virtual void CalculateEddyViscosity(const CGeometry* geometry,
                                      CConfig* config,
                                      unsigned long iPoint,
                                      su2double mean_eddy_visc,
                                      su2double** eddy_viscosity) const = 0;
 protected:
  const unsigned short nDim; /*!< \brief The number of physical dimensions */
  const su2double* FlowPrimVar; /*!< \brief The flow primitive variables */
  const su2double* TurbVar; /*!< \brief The mean turbulence solution variables */
};

/*!
 * \class CM43Model
 * \brief Haering's M43 energy transfer model
 * \author S. Haering and C. Pederson
 */
class CM43Model : public CFluctuatingStress {
 public:
  /*!
   * \brief Basic constructor.
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] config - Definition of the particular problem.
   */
  CM43Model(unsigned short val_nDim,
            CConfig* config);
  /*!
   * \brief Destructor
   */
  ~CM43Model();

  /**
   * \brief Calculate the eddy viscosity.
   *
   * The eddy viscosity computed is anisotropic and is a rank-2 tensor.
   * The mean eddy viscosity is only used if the cell has too high of an
   * aspect ratio to use the M43 model. The M43 model will replace the
   * typical M43 eddy viscosity with the isotropic RANS eddy viscosity.
   *
   * \param[in] geometry - The geometry
   * \param[in] config - The config settings
   * \param[in] iPoint - The current point to be used
   * \param[in] mean_eddy_visc - The mean (i.e. RANS) eddy viscosity.
   * \param[out] eddy_viscosity - The computed anisotropic eddy viscosity
   */
  void CalculateEddyViscosity(const CGeometry* geometry,
                              CConfig* config,
                              unsigned long iPoint,
                              su2double mean_eddy_visc,
                              su2double** eddy_viscosity) const;
 private:
  su2double** delta; /*!< \brief The Kroneckor delta. */
};

/*!
 * \class CNoStressModel
 * \brief A model that always returns 0 for the eddy viscosity.
 *
 * This model serves as a type of "null object"
 *
 * \author C. Pederson
 */
class CNoStressModel : public CFluctuatingStress {
 public:
  /*!
   * \brief Basic constructor.
   * \param[in] val_nDim - Number of dimensions of the problem.
   */
  CNoStressModel(unsigned short val_nDim);

  /**
   * \brief Calculate the eddy viscosity (which is always 0)
   *
   * \param[in] geometry - The geometry
   * \param[in] config - The config settings
   * \param[in] iPoint - The current point to be used
   * \param[in] mean_eddy_visc - The mean (i.e. RANS) eddy viscosity.
   * \param[out] eddy_viscosity - The computed anisotropic eddy viscosity
   */
  void CalculateEddyViscosity(const CGeometry* geometry,
                              CConfig* config,
                              unsigned long iPoint,
                              su2double mean_eddy_visc,
                              su2double** eddy_viscosity) const;
};

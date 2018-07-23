/*!
 * \file fluctuating_stress.hpp
 * \brief Interface for the fluctuating stress for hybrid RANS/LES
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

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/mpi_structure.hpp"
#include "../include/solver_structure.hpp"
#include "../include/numerics_structure.hpp"

#pragma once

/*!
 * \class CFluctuatingStress
 * \brief An abstract class, used to model the energy transfer from
 *        resolved to modeled scales.
 * \author C. Pederson
 */
class CFluctuatingStress {
 public:
  /**
   * \brief Basic constructor
   * \param solver_container
   * \param geometry
   */
  CFluctuatingStress(CSolver **solver_container, CGeometry* geometry);
  virtual ~CFluctuatingStress();
  /**
   * \brief Set the flow primitive variables
   * \param[in] val_primitive - An array containing primitive variables
   *     from the flow
   */
  void SetPrimitive(const su2double* val_primitive);
  /**
   * \brief Set the mean of turbulence variables
   * \param turb_vars - An array containing the mean turbulence variables
   */
  void SetMeanTurbVars(const su2double* turb_vars);
  /**
   * \brief Calcualte the eddy viscosity.
   *
   * This function is pure virtual and must be overloaded in derived classes.
   * The eddy viscosity computed is anisotropic and is a rank-2 tensor.
   *
   * \param[in] geometry - The geometry
   * \param[in] config - The config settings
   * \param[in] iPoint - The current point to be used
   * \param[out] eddy_viscosity - The computed anisotropic eddy viscosity
   */
  virtual void CalculateEddyViscosity(const CGeometry* geometry,
                                      CConfig* config,
                                      unsigned short iPoint,
                                      su2double** eddy_viscosity) const = 0;
 protected:
  const unsigned short nDim;
  const su2double* FlowPrimVar;
  const su2double* MeanTurbVar;
};

class CM43Model : public CFluctuatingStress {
 public:
  CM43Model(CSolver **solver_container, CGeometry* geometry, CConfig* config);
  void CalculateEddyViscosity(const CGeometry* geometry,
                              CConfig* config,
                              unsigned short iPoint,
                              su2double** eddy_viscosity) const;
 private:
  unsigned short density_index;
};

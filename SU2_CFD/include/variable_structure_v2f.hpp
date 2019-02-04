/*!
 * \file variable_structure.hpp
 * \brief Headers of the main subroutines for storing all the variables for
 *        each kind of governing equation (direct, adjoint and linearized).
 *        The subroutines and functions are in the <i>variable_structure.cpp</i> file.
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

#include "variable_structure.hpp"

/*!
 * \class CTurbKEVariable
 * \brief Main class for defining the variables of the turbulence model.
 * \ingroup Turbulence_Model
 * \author S. Haering
 * \version 4.3.x "Cardinal"
 */
class CTurbKEVariable : public CTurbVariable {

protected:
  su2double sigma_e, sigma_k, sigma_z, C_e1o, C_e2, C1, C_2p, C_T, C_L, C_eta;
  su2double Tm,		/*!< \brief T_m k-eps. */
    Lm,		        /*!< \brief L_m k-eps */
    Re_T,
    Production;   /*!< \brief Production of TKE */

public:
  /*!
   * \brief Constructor of the class.
   */
  CTurbKEVariable(void);

  /*!
   * \overload
   * \param[in] val_rho_kine - Turbulent variable value (initialization value).
   * \param[in] val_rho_omega - Turbulent variable value (initialization value).
   * \param[in] val_muT - Turbulent variable value (initialization value).
   * \param[in] val_nDim - Number of dimensions of the problem.
   * \param[in] val_nvar - Number of variables of the problem.
   * \param[in] constants -
   * \param[in] config - Definition of the particular problem.
   */
  CTurbKEVariable(su2double val_rho_kine, su2double val_rho_epsi,
                  su2double val_zeta, su2double val_f,
                  su2double val_muT, su2double val_Tm, su2double val_Lm,
                  unsigned short val_nDim, unsigned short val_nvar,
                  su2double *constants, CConfig *config);

  /*!
   * \brief Destructor of the class.
   */
  ~CTurbKEVariable(void);

  /**
   * \brief Get the large-eddy timescale of the turbulence
   * \return The large-eddy timescale of the turbulence.
   */
  su2double GetTurbTimescale(void) const;

  /**
   * \brief Get the large-eddy lengthscale of the turbulence
   * \return The large-eddy lengthscale of the turbulence
   */
  su2double GetTurbLengthscale(void) const;

  /*!
   * \brief Get the component anisotropy ratio (max-to-min)
   * \return The Reynolds stress component anisotropy ratio (max-to-min)
   */
  su2double GetAnisoRatio(void);

  /**
   * \brief Sets the large-eddy lengthscale and the large-eddy timescale
   * \param[in] val_turb_T - Large eddy timescale of the turbulence
   * \param[in] val_turb_L - Large eddy lengthscale of the turbulence
   */
  void SetTurbScales(su2double val_turb_T, su2double val_turb_L);

  /*!
   * \brief Set the production of turbulent kinetic energy.
   * \param[in] val_production - Production of turbulent kinetic energy.
   */
  void SetProduction(su2double val_production);

  /*!
   * \brief Get the production of turbulent kinetic energy.
   * \return Production of turbulent kinetic energy.
   */
  su2double GetProduction(void) const;
};

#include "variable_structure_v2f.inl"

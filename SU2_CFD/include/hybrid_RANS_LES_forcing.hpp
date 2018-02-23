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

#include "../../Common/include/config_structure.hpp"

/*!
 * \class CHybridForcing
 * \brief Class for defining the forcing needed in a hybrid RANS/LES model
 * \author C. Pederson
 * \version 5.0.0 "Raven"
 *
 * This class could inherit from CNumerics, but there's a **lot** of extra
 * code in numerics that forcing doesn't need  The choice is to have a
 * simpler class with a clearer interface at the cost of some code duplication.
 *
 * This class is a bit of a heavyweight class, requiring a lot of getters
 * and setters for relatively little calculation.  If this class frequently
 * changes, then the single-responsibility principle  (each class should only
 * have one reason to change) justifies the existence of this class. But if
 * this class becomes more stable, it should be dismantled into simpler
 * inline functions.
 */

class CHybridForcing {
 public:
  CHybridForcing(unsigned short nDim);
  ~CHybridForcing();

  /*!
   * \brief Set coordinates of the point.
   * \param[in] val_coord - Coordinates of the point.
   */
  void SetCoord(su2double* val_coord);

  void SetUnstTime(su2double val_time);

  /*!
   * \brief Set the turbulent timescale
   * \param[in] val_turb_T - Turbulent timescale at point i
   */
  void SetTurbTimescale(su2double val_turb_T);

  /*!
   * \brief Set the turbulent timescale
   * \param[in] val_turb_T - Turbulent lengthscale at point i
   */
  void SetTurbLengthscale(su2double val_turb_L);

  /*!
   * \brief Set the value of the hybrid RANS/LES blending variable.
   * \param[in] val_hybrid_param - Value of the hybrid parameter(s)
   */
  void SetHybridParameter(su2double* val_hybrid_param);

  void SetPrimitive(su2double* val_prim_vars);

  void SetPrimVarGradient(su2double **val_primvar_grad);

  /*!
   * \brief Set the value of the turbulent variable.
   * \param[in] val_turbvar - Value of the turbulent variable
   */
  void SetTurbVar(su2double *val_turbvar);

  /*!
   * \brief Set the hybrid source terms necessary.
   * \param val_S_terms - In this case, S_alpha and S_cf
   */
  void SetHybridSourceTerms(su2double* val_S_terms);

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
  su2double TransformCoords(su2double* x, unsigned short iDim,
                            su2double val_time);

  /**
   * \brief Use the coordinates of a point to build the forcing stress at
   *        that point.
   *
   * This method is only left public for testing purposes. It is usually
   * not necessary to call it manually.  See CalculateForcing instead.
   *
   * \param[in] x - Physical coordinates
   * \param[out] val_tau_F - The unscaled forcing stress
   */
  void BuildForcingStress(su2double* x, su2double** val_tau_F);
  /**
   * \brief Calculate all of the necessary forcing terms.
   *
   * Instead of recomputing the forcing terms for the resolved flow,
   * RANS, and hybrid parameter equations, just compute them once by
   * calling this function.  This function serves as a wrapper to some of
   * the lower level functions, such as coordinate shifting and building
   * the stress tensor.
   *
   * You must set all the appropriate variables (such as Coord and
   * PrimVarGradient) before calling this function.
   */
  void CalculateForcing();

  /**
   * \brief Getter for the forcing stress.
   *
   * CalculateForcing() must be executed before this contains any valid data.
   *
   * \return The forcing stress.
   */
  su2double** GetStress();

  /**
   * \brief Getter for the forcing stress.
   *
   * CalculateForcing() must be executed before this contains any valid data.
   *
   * \return The production of turbulent kinetic energy due to forcing
   */
  su2double GetProduction();

  /**
   * \brief Getter for the ratio of the production due to forcing and a rescaled
   * production term.
   *
   * \return The ratio of P_F_unscaled / P_lim
   */
  su2double GetForcingRatio();

 private:
  unsigned short nVar, nDim;

  /*--- External variables ---*/
  su2double HybridParam;
  su2double TurbT;
  su2double TurbL;
  su2double* Original_Coord, *Shifted_Coord;
  su2double* PrimVars;
  su2double** PrimVar_Grad;
  su2double* TurbVar;
  su2double S_alpha;
  su2double S_cf;
  su2double time;

  /*--- Problem variables --*/
  su2double** tau_F, **tau_F_unscaled;
  su2double P_F, P_F_unscaled;
  su2double P_lim;

  su2double CalculateCoefficient();
};

#include "hybrid_RANS_LES_forcing.inl"

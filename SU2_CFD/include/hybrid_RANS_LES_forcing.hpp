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
//#include "solver_structure.hpp" // moved to hybrid_RANS_LES_forcing.cpp to avoid circular includes here

// Forward declarations to resolve circular dependencies
class CSolver;


/*!
 * \class CHybridForcingAbstractBase
 * \brief Abstract base class for hybrid forcing.
 * \author T. A. Oliver
 * \version 5.0.0 "Raven"
 *
 * This class provides a uniform interface to various forcing options,
 * which are implemented as derived classes.
 */
class CHybridForcingAbstractBase {
 public:
  CHybridForcingAbstractBase(unsigned short nDim,
                             unsigned long nPoint,
                             unsigned long nPointDomain);
  virtual ~CHybridForcingAbstractBase();

  /**
   * \brief Transform physical coordinates into coordinates used for forcing.
   *
   * This method is only left public for testing purposes. It is usually
   * not necessary to call it manually.
   *
   * \param[in] x - Physical coordinates
   * \param[in] iDim - The component to be transformed
   * \param[in] val_time - The current time
   * \return A component of the shifted coordinates.
   */
  su2double TransformCoords(su2double x, su2double mean_velocity,
                            su2double time) const;


  /**
   * \brief Evaluate forcing field (pure virtual)
   */
  virtual void ComputeForcingField(CSolver** solver, CGeometry *geometry,
                                   CConfig *config) = 0;

  virtual const su2double* GetForcingVector(unsigned long iPoint) const = 0;

 protected:
  const unsigned short nDim,   /*!< \brief Number of dimensions of the problem. */
  nVar,                        /*!< \brief Number of variables of the problem. */
  nVarGrad;                    /*!< \brief Number of variables for deallocating the LS Cvector. */
  const unsigned long nPoint,  /*!< \brief Number of points of the computational grid. */
  nPointDomain;                /*!< \brief Number of points of the computational grid. */
};


/*!
 * \class CHybridForcingTG0
 * \brief Class for defining the forcing needed in a hybrid RANS/LES model
 * \author T. Oliver, C. Pederson
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
class CHybridForcingTG0 : public CHybridForcingAbstractBase{
 public:
  CHybridForcingTG0(unsigned short nDim, unsigned long nPoint,
                    unsigned long nPointDomain);
  CHybridForcingTG0(CGeometry* geometry, const CConfig* config);
  ~CHybridForcingTG0();

  /**
   * \brief Evaluate baseline TG field at point.
   *
   * This method is only left public for testing purposes. It is usually
   * not necessary to call it manually.
   *
   * \param[in]  x - Forcing coordinates (i.e., result of TransformCoords)
   * \param[in]  Lsgs - SGS length scales
   * \param[in]  Lmesh - Mesh length scales
   * \param[in]  D - Domain lengths in periodic directions
   * \param[in]  dwall - Distance to nearest wall
   * \param[out] b - TG velocity at point.
   */
  void SetTGField(const su2double* x, su2double Lsgs,
                  const su2double* Lmesh, const su2double* D,
                  su2double dwall, su2double* h) const;

  void SetAxiTGField(const su2double* x, const su2double Lsgs,
		     const su2double* Lmesh, const su2double* D,
		     const su2double dwall, su2double* h) const;


  su2double GetTargetProduction(su2double v2,
                                su2double T,
                                su2double alpha) const;

  su2double ComputeScalingFactor(su2double Ftar,
                                 su2double resolution_adequacy,
                                 su2double alpha,
                                 su2double alpha_kol,
                                 su2double PFtest) const;

  void ComputeForcingField(CSolver** solver, CGeometry *geometry,
                           CConfig *config);

  const su2double* GetForcingVector(unsigned long iPoint) const;

 protected:

  const su2double forcing_scale;  /*!< \brief The forcing vortices will be of period N*L, where N is the forcing scale and L is the turbulent length scale. */
  const su2double C_F;  /*!< \brief The overall strength of the forcing. */

  su2double** node;
};

#include "hybrid_RANS_LES_forcing.inl"

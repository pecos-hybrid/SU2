/*!
 * \file bulk_viscosity_model.hpp
 * \brief Headers of the bulk viscosity model
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../../Common/include/mpi_structure.hpp"

#include <stdio.h>
#include <string.h>
#include <iostream>
#include <string>
#include <cmath>

using namespace std;

#include "../../Common/include/config_structure.hpp"
#include "transport_model.hpp"

/*!
 * \brief Abstract class for defining the Bulk Viscosity Model as
 *   a child class for each particular Model
 * \author: C. Pederson
 */
class CBulkViscosityModel {
protected:
 su2double kappa_ratio; /*!< \brief Ratio of bulk viscosity to laminar viscosity*/

public:

	  /*!
		 * \brief Constructor of the class.
		 */
		CBulkViscosityModel(void);

		/*!
		 * \brief Destructor of the class.
		 */
		virtual ~CBulkViscosityModel(void);

		/*!
		 * \brief Get the ratio of bulk viscosity to shear viscosity;
		 */
		su2double GetBulkViscosityRatio() const { return kappa_ratio; };

		/*!
		 * \brief Get ratio of the second viscosity coefficient to shear
     *     viscosity.
		 */
		su2double GetSecondViscosityCoeffRatio() const { return kappa - TWO3; };

    /*!
     * \brief Set Viscosity.
     */
    virtual void SetBulkViscosity(su2double T) = 0;
};


/*!
 * \brief Child class for defining a constant bulk viscosity model
 *
 * If the bulk viscosity ratio is set to 0, then this model is
 * equivalent to Stokes' assumption.
 *
 * \author: C. Pederson
 */
class CConstantBulkViscosity : public CBulkViscosity {
public:

	   /*!
		 * \brief Constructor of the class.
		 */
		CConstantBulkViscosity(su2double val_kappa_ratio);

		/*!
		 * \brief Destructor of the class.
		 */
		~CConstantBulkViscosity(void);

    /*!
     * \brief Set Viscosity.
     */
    void SetBulkViscosity() {};
};

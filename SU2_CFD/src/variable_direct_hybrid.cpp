/*!
 * \file variable_direct_hybrid.cpp
 * \brief Definition of the hybrid parameter(s) for hybrid RANS/LES
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
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/variable_structure.hpp"

CHybridVariable::CHybridVariable(void) {
}

CHybridVariable::~CHybridVariable(void) {};

CHybridVariable::CHybridVariable(unsigned short val_nDim,
                                     unsigned short val_nvar,
                                     CConfig *config)
  : CVariable(val_nDim, val_nvar, config) {

    unsigned short iVar;

    /*--- Allocate space for the limiter ---*/

    Limiter = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++)
      Limiter[iVar] = 0.0;

    Solution_Max = new su2double [nVar];
    Solution_Min = new su2double [nVar];
    for (iVar = 0; iVar < nVar; iVar++) {
      Solution_Max[iVar] = 0.0;
      Solution_Min[iVar] = 0.0;
    }

    /*--- Initialize hybrid variables in balanced RANS mode ---*/

    Resolution_Adequacy = 1;
    RANS_Weight = 1;
    Forcing_Ratio = 1.0;
    S_terms[0] = 0.0;
    S_terms[1] = 1.0;

}

CHybridConvVariable::CHybridConvVariable(su2double val_hybrid_param,
                                             unsigned short val_nDim,
                                             unsigned short val_nvar,
                                             CConfig *config)
  : CHybridVariable(val_nDim, val_nvar, config) {

    bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                      (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

    /*--- Initialization of hybrid parameter ---*/
    Solution[0] = val_hybrid_param;    Solution_Old[0] = val_hybrid_param;

    /*--- Allocate and initialize solution for the dual time strategy ---*/
    if (dual_time) {
      Solution_time_n[0]  = val_hybrid_param;
      Solution_time_n1[0] = val_hybrid_param;
    }

    /*--- Initialize resolution adequacy in balanced balanced RANS mode ---*/

    Resolution_Adequacy = 1;
    RANS_Weight = 1;

}

CHybridConvVariable::~CHybridConvVariable(void) {};



/*!
 * \file variable_direct_turbulent.cpp
 * \brief Definition of the solution fields.
 * \author F. Palacios, A. Bueno
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

#include "../include/variable_structure_v2f.hpp"

CTurbKEVariable::CTurbKEVariable(void) : CTurbVariable() { }

CTurbKEVariable::CTurbKEVariable(su2double val_kine, su2double val_epsi,
                                 su2double val_zeta, su2double val_f,
                                 su2double val_muT, su2double val_Tm,
                                 su2double val_Lm, unsigned short val_nDim,
                                 unsigned short val_nvar,
                                 su2double *constants, CConfig *config)
  :
  CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of variables ---*/
  Solution[0] = val_kine; Solution_Old[0] = val_kine;
  Solution[1] = val_epsi;	Solution_Old[1] = val_epsi;
  Solution[2] = val_zeta; Solution_Old[2] = val_zeta;
  Solution[3] = val_f;	Solution_Old[3] = val_f;
  Tm  = val_Tm;
  Lm  = val_Lm;

  /*--- Initialization of eddy viscosity ---*/  
  muT = val_muT;

  /*--- Allocate and initialize solution for the dual time strategy ---*/
  if (dual_time) {
    Solution_time_n[0]  = val_kine; Solution_time_n[1]  = val_epsi;
    Solution_time_n1[0] = val_kine; Solution_time_n1[1] = val_epsi;
    Solution_time_n[2]  = val_zeta; Solution_time_n[3]  = val_f;
    Solution_time_n1[2] = val_zeta; Solution_time_n1[3] = val_f;
  }

}

CTurbKEVariable::~CTurbKEVariable(void) {
}

su2double CTurbKEVariable::GetAnisoRatio(void) {
  // XXX: This floor is arbitrary.
  const su2double TKE_MIN = EPS;
  return TWO3*Solution[0]/max(TKE_MIN, Solution[2]);
}

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

#include "../include/variable_structure.hpp"
#include "../include/variable_structure_v2f.hpp"

/*--- Default values of static members
 * These can be changed by calling CTurbKEVariable::SetConstants() ---*/
su2double CTurbKEVariable::C_mu    = 0.22;
su2double CTurbKEVariable::C_T     = 6.0;
su2double CTurbKEVariable::C_eta   = 70.0;

CTurbKEVariable::CTurbKEVariable(void) : CTurbVariable() { }

CTurbKEVariable::CTurbKEVariable(su2double val_kine, su2double val_epsi,
                                 su2double val_zeta, su2double val_f,
                                 su2double val_muT, su2double val_Tm,
                                 su2double val_Lm, unsigned short val_nDim,
                                 unsigned short val_nvar, CConfig *config)
    : CTurbVariable(val_nDim, val_nvar, config) {

  bool dual_time = ((config->GetUnsteady_Simulation() == DT_STEPPING_1ST) ||
                    (config->GetUnsteady_Simulation() == DT_STEPPING_2ND));

  /*--- Initialization of variables ---*/
  Solution[0] = val_kine; Solution_Old[0] = val_kine;
  Solution[1] = val_epsi;	Solution_Old[1] = val_epsi;
  Solution[2] = val_zeta; Solution_Old[2] = val_zeta;
  Solution[3] = val_f;	Solution_Old[3] = val_f;
  timescale  = val_Tm;
  lengthscale  = val_Lm;
  typical_length = lengthscale;
  typical_time = timescale;

  /*--- Initialization of eddy viscosity ---*/  
  muT = val_muT;

  /*--- Initialize kolmogorov scales to arbitrary values ---*/

  alpha_kol = 0.01;
  kol_time = 0.01 * timescale;
  kol_length = 0.01 * lengthscale;

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

 void CTurbKEVariable::SetTurbScales(const su2double nu,
                                     const su2double S,
                                     const su2double VelMag,
                                     const su2double L_inf) {
  /*--- Scalars ---*/
  const su2double kine = Solution[0];
  const su2double epsi = Solution[1];
  const su2double v2   = Solution[2];

  /*--- Relevant scales ---*/
  const su2double scale = EPS;

  /*--- Clipping to avoid nonphysical quantities
   * We keep "tke_positive" in order to allow tke=0 but clip negative
   * values.  If tke is some small nonphysical quantity (but not zero),
   * then it is possible for L1 > L3 and T1 > T3 when the viscous
   * scales are smaller than this nonphysical limit. ---*/

  const su2double tke_lim = max(kine, scale*VelMag*VelMag);
  const su2double tke_positive = max(kine, 0.0);
  const su2double tdr_lim = max(epsi, scale*VelMag*VelMag*VelMag/L_inf);
  const su2double zeta_lim = max(v2/tke_lim, scale);
  const su2double S_FLOOR = scale;

  //--- Model time scale ---//

  typical_time = tke_positive/tdr_lim;
  // sqrt(3) instead of sqrt(6) because of sqrt(2) factor in S
  stag_time    = 0.6/max(sqrt(3.0)*C_mu*S*zeta_lim, S_FLOOR);
  kol_time     = C_T*sqrt(nu/tdr_lim);
  timescale = max(min(typical_time, stag_time), kol_time);

  //--- Model length scale ---//
  typical_length = pow(tke_positive,1.5)/tdr_lim;
  // sqrt(3) instead of sqrt(6) because of sqrt(2) factor in S
  stag_length    = sqrt(tke_positive)/max(sqrt(3.0)*C_mu*S*zeta_lim, S_FLOOR);
  kol_length     = C_eta*pow(pow(nu,3.0)/tdr_lim,0.25);
  //... multiply by C_L in source numerics
  lengthscale = max(min(typical_length, stag_length), kol_length);
}

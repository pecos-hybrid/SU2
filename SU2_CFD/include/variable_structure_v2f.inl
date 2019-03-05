/*!
 * \file variable_structure.inl
 * \brief In-Line subroutines of the <i>variable_structure.hpp</i> file.
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

inline su2double CTurbKEVariable::GetTurbTimescale() const{
  return timescale;
}

inline su2double CTurbKEVariable::GetTurbLengthscale() const {
 return lengthscale;
}

inline void CTurbKEVariable::SetProduction(su2double val_production) {
  Production = val_production;
}

inline su2double CTurbKEVariable::GetProduction(void) const {
  return Production;
}

inline void CTurbKEVariable::SetKolKineticEnergyRatio(const su2double nu) {
  const su2double TKE_MIN = 1.0E-8;
  const su2double ktot = max(Solution[0], TKE_MIN);
  const su2double tdr = Solution[1];
  const su2double Cnu = 1.0;
  alpha_kol = min(Cnu*std::sqrt(nu*tdr)/ktot, 1.0);
}

inline su2double CTurbKEVariable::GetAnisoRatio(void) {
  // XXX: This floor is arbitrary.
  const su2double TKE_MIN = EPS;
  return TWO3*Solution[0]/max(TKE_MIN, Solution[2]);
}

inline su2double CTurbKEVariable::GetTypicalLengthscale(void) const {
  return typical_length;
}

inline su2double CTurbKEVariable::GetTypicalTimescale(void) const {
  return typical_time;
}

inline su2double CTurbKEVariable::GetKolLengthscale(void) const {
  return kol_length;
}

inline su2double CTurbKEVariable::GetKolTimescale(void) const {
  return kol_time;
}

inline su2double CTurbKEVariable::GetKolKineticEnergyRatio(void) const {
  return alpha_kol;
}

inline void CTurbKEVariable::SetConstants(const su2double* constants)  {
  C_mu    = constants[0];
  C_T     = constants[8];
  C_eta   = constants[10];
}


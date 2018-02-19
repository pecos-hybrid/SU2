/*!
 * \file hybrid_RANS_LES_forcing.inl
 * \brief Inline functions for the hybrid RANS/LES forcing terms
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

inline void CHybridForcing::SetCoord(su2double* val_coord) {
  Original_Coord = val_coord;
}

inline void CHybridForcing::SetTurbLengthscale(su2double val_L_m) {
  TurbL = val_L_m;
}

inline void CHybridForcing::SetTurbTimescale(su2double val_T_m) {
  TurbT = val_T_m;
}

inline void CHybridForcing::SetHybridParameter(su2double* val_hybrid_param) {
  HybridParam = val_hybrid_param[0];
}

inline void CHybridForcing::SetPrimitive(su2double* val_prim_vars) {
    PrimVars = val_prim_vars;
}

inline void CHybridForcing::SetPrimVarGradient(su2double **val_primvar_grad) {
  PrimVar_Grad = val_primvar_grad;
}

inline void CHybridForcing::SetTurbVar(su2double *val_turbvar) {
  TurbVar = val_turbvar;
}

inline void CHybridForcing::SetHybridSourceTerms(su2double* val_S_terms) {
  S_alpha = val_S_terms[0];
  S_cf = val_S_terms[1];
}

inline su2double** CHybridForcing::GetStress() {
  return tau_F;
}

inline su2double CHybridForcing::GetProduction() { return P_F; };

inline su2double CHybridForcing::GetForcingRatio() {
  return P_F_unscaled / P_lim;
};

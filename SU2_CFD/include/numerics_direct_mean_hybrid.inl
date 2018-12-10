/*!
 * \file numerics_direct_mean_hybrid.inl
 * \brief In-Line subroutines of the <i>numerics_direct_mean_hybird.hpp</i> file.
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

inline void CAvgGrad_Hybrid::SetAniso_Eddy_Viscosity(su2double** aniso_eddy_viscosity_i,
                                                     su2double** aniso_eddy_viscosity_j) {
  Aniso_Eddy_Viscosity_i = aniso_eddy_viscosity_i;
  Aniso_Eddy_Viscosity_j = aniso_eddy_viscosity_j;
}


inline void CAvgGrad_Hybrid::SetPrimitive_Average(su2double* val_primvar_average_i,
                                                  su2double* val_primvar_average_j) {
  PrimVar_Average_i = val_primvar_average_i;
  PrimVar_Average_j = val_primvar_average_j;
}

inline void CAvgGrad_Hybrid::SetPrimVarGradient_Average(su2double **val_primvar_grad_i,
                                                        su2double **val_primvar_grad_j) {
  PrimVar_Grad_Average_i = val_primvar_grad_i;
  PrimVar_Grad_Average_j = val_primvar_grad_j;
}

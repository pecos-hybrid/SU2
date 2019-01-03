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

inline su2double CHybridForcingAbstractBase::TransformCoords(
                                               const su2double x,
                                               const su2double mean_velocity,
                                               const su2double time,
                                               const su2double timescale) {
  const su2double N_T = 4.0;
  return x + mean_velocity*fmod(time, N_T*timescale);
}

inline void CHybridForcingAbstractBase::SetTGField(const su2double* x,
                                                   const su2double* L,
                                                   su2double* b) {
  const su2double pi = atan(1.0)*4.0;
  const int A = 3, B = -1, C = -2;
  const su2double a[3] = {2*pi/L[0], 2*pi/L[1], 2*pi/L[2]};

  b[0] = A * cos(a[0]*x[0]) * sin(a[1]*x[1]) * sin(a[2]*x[2]);
  b[1] = B * sin(a[0]*x[0]) * cos(a[1]*x[1]) * sin(a[2]*x[2]);
  b[2] = C * sin(a[0]*x[0]) * sin(a[1]*x[1]) * cos(a[2]*x[2]);
}

inline void CHybridForcing::SetStreamFunc(const su2double* x,
                                          const su2double* L,
                                          su2double* h) {
  const su2double pi = atan(1.0)*4.0;
  const int D = -1, E = 1, F = -2;
  const su2double a[3] = {2*pi/L[0], 2*pi/L[1], 2*pi/L[2]};

  h[0] = D * sin(a[0]*x[0]) * cos(a[1]*x[1]) * cos(a[2]*x[2]);
  h[1] = E * cos(a[0]*x[0]) * sin(a[1]*x[1]) * cos(a[2]*x[2]);
  h[2] = F * cos(a[0]*x[0]) * cos(a[1]*x[1]) * sin(a[2]*x[2]);
}


inline su2double CHybridForcing::GetTargetProduction(const su2double k_sgs,
                                                     const su2double dissipation,
                                                     const su2double resolution_adequacy,
                                                     const su2double alpha,
                                                     const su2double laminar_viscosity) {
  const su2double C_kol = 10.0;
  const su2double alpha_kol =
      C_kol*alpha*sqrt(laminar_viscosity * dissipation)/k_sgs;
  const su2double S_r = tanh(log(resolution_adequacy));
  const su2double D_r = (1.0 + alpha_kol - alpha) * min(S_r, 0.0);
  const su2double F_r = alpha * max(S_r, 0.0);
  return -dissipation*(S_r - D_r - F_r);
}

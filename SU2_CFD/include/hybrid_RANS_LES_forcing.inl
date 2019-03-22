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
                                               const su2double time) const {
  return x + mean_velocity*time;
}

inline su2double CHybridForcingTG0::GetTargetProduction(const su2double v2,
                                                        const su2double Tsgs,
                                                        const su2double alpha) const {
  return C_F * std::sqrt(alpha*v2)/Tsgs;
}

inline const su2double* CHybridForcingTG0::GetForcingVector(unsigned long iPoint) const {
  return node[iPoint];
}

inline su2double CHybridForcingTG0::ComputeScalingFactor(
                     const su2double Ftar,
                     const su2double resolution_adequacy,
                     const su2double alpha,
                     const su2double alpha_kol,
                     const su2double PFtest) const {

  su2double eta = 0.0;

  // TODO: Compare this with Sigfried's improved version once channel
  // validation is successful.
  if ( (PFtest >= 0.0) && (resolution_adequacy < 1.0) ) {
    const su2double Sr = std::tanh(1.0 - 1.0/sqrt(resolution_adequacy));
    if (alpha <= alpha_kol) {
      eta = -Ftar * Sr * (alpha - alpha_kol);
    } else {
      eta = -Ftar * Sr;
    }
  }

  return eta;
}

inline void CHybridForcingTG0::SetTGField(
                const su2double* x, const su2double Lsgs,
                const su2double* Lmesh, const su2double* D,
                const su2double dwall, su2double* h) const {

  const su2double A = 1./3., B = -1.0, C = 2./3.;
  su2double a[3];

  for (unsigned int ii=0; ii<3; ii++) {
    const su2double ell = std::min(Lsgs, dwall);
    const su2double elllim = std::max(ell, 2.0*Lmesh[ii]);

    if (D[ii] > 0.0) {
      const su2double denom = round(D[ii]/std::min(elllim, D[ii]));
      a[ii] = M_PI/(D[ii]/denom);
    } else {
      a[ii] = M_PI/elllim;
    }
  }

  // h[0] = A * cos(a[0]*x[0]) * sin(a[1]*x[1]) * sin(a[2]*x[2]);
  // h[1] = B * sin(a[0]*x[0]) * cos(a[1]*x[1]) * sin(a[2]*x[2]);
  // h[2] = C * sin(a[0]*x[0]) * sin(a[1]*x[1]) * cos(a[2]*x[2]);

  //// CHANNEL HACK
  //h[0] = A * cos(a[0]*x[0]) * sin(a[1]*(x[1]-1.0)) * sin(a[2]*x[2]-M_PI/2);
  //h[1] = B * sin(a[0]*x[0]) * cos(a[1]*(x[1]-1.0)) * sin(a[2]*x[2]-M_PI/2);
  //h[2] = C * sin(a[0]*x[0]) * sin(a[1]*(x[1]-1.0)) * cos(a[2]*x[2]-M_PI/2);
  // CHANNEL HACK
  h[0] = A * cos(a[0]*x[0]) * sin(a[1]*(x[1]-1.0)) * sin(a[2]*x[2]);
  h[1] = B * sin(a[0]*x[0]) * cos(a[1]*(x[1]-1.0)) * sin(a[2]*x[2]);
  h[2] = C * sin(a[0]*x[0]) * sin(a[1]*(x[1]-1.0)) * cos(a[2]*x[2]);
}


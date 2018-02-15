/*!
 * \file hybrid_RANS_LES_forcing.cpp
 * \brief Implementation of the hybrid RANS/LES forcing terms
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

#include "../include/hybrid_RANS_LES_forcing.hpp"

CHybridForcing::CHybridForcing(unsigned short nDim) : nDim(nDim) {

  Shifted_Coord = new su2double[nDim];

  tau_F = new su2double*[nDim];
  tau_F_unscaled = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    tau_F[iDim] = new su2double[nDim];
    tau_F_unscaled[iDim] = new su2double[nDim];
  }

}

CHybridForcing::~CHybridForcing() {

  delete [] Shifted_Coord;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] tau_F[iDim];
    delete [] tau_F_unscaled[iDim];
  }
  delete [] tau_F;
  delete [] tau_F_unscaled;

}


su2double CHybridForcing::TransformCoords(su2double* x,
                                          unsigned short iDim,
                                          su2double val_time) {
  const su2double pi = std::atan(1.0)*4.0;
  return 2.0*pi/(std::pow(HybridParam, 1.5)*TurbL) *
      (x[iDim] + PrimVars[iDim+1] * std::fmod(val_time, (HybridParam*TurbT)));
}

void CHybridForcing::BuildForcingStress(su2double* x, su2double** val_tau_F) {

  unsigned short nDim = 3;
  su2double g[3];

  g[0] =  1.0*cos(x[0])*sin(x[1])*sin(x[2]);
  g[1] =  1.0*sin(x[0])*cos(x[1])*sin(x[2]);
  g[2] = -2.0*sin(x[0])*sin(x[1])*cos(x[2]);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      val_tau_F[iDim][jDim] = g[iDim]*g[jDim] +
          g[iDim]*PrimVars[jDim+1] + g[jDim]*PrimVars[iDim+1];
    }
  }
}

void CHybridForcing::CalculateForcing() {

  const su2double C_F = 1.0e-2;

  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Shifted_Coord[iDim] = TransformCoords(Original_Coord, iDim, time);

  BuildForcingStress(Shifted_Coord, tau_F_unscaled);

  P_F_unscaled = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      P_F_unscaled += tau_F_unscaled[iDim][jDim] *
          0.5*(PrimVar_Grad[iDim][jDim] + PrimVar_Grad[jDim][iDim]);
    }
  }

  int sign_P_F_unscaled = (P_F_unscaled < 0) ? -1 : 1;
  P_lim = sign_P_F_unscaled * max(abs(P_F_unscaled), C_F*TurbVar[1]);
  su2double coefficient = CalculateCoefficient();

  /*--- Calculate forcing stress tensor ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      tau_F[iDim][jDim] = coefficient * tau_F_unscaled[iDim][jDim];
    }
  }

  /*--- Calculate production due to forcing ---*/

  P_F = coefficient * P_F_unscaled;

}

su2double CHybridForcing::CalculateCoefficient() {

  // TODO: Pull C_alpha from the config file
  // XXX: This assumes that k is the first turbulent variable
  const su2double C_alpha = 0.1;
  return C_alpha * TurbVar[0] / (HybridParam*TurbT) * (S_alpha-S_cf) / P_lim;

}

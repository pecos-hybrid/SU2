/*!
 * \file hybrid_RANS_LES_forcing.inl
 * \brief Inline functions for the hybrid RANS/LES forcing terms
 * \author C. Pederson
 * \version 6.2.0 "Falcon"
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
  return C_F * sqrt(alpha*v2)/Tsgs;
}

inline const su2double* CHybridForcingTG0::GetForcingVector(unsigned long iPoint) const {
  return node[iPoint];
}

inline void CHybridForcingTG0::SetTGField(
                const su2double* x, const su2double Lsgs,
                const su2double* D,
                const su2double dwall, su2double* h) const {

  //const su2double A = 1./3., B = -1.0, C = 2./3.;
  const su2double A = 1./3., B = 2./3., C = -1.0;
  su2double a[3];

  if (dwall > 0 && Lsgs > 0) {
    for (unsigned int ii=0; ii<3; ii++) {
      const su2double ell = min(Lsgs, dwall);

      if (D[ii] > 0.0) {
        const su2double denom = round(D[ii]/min(ell, D[ii]));
        a[ii] = M_PI/D[ii]*denom;
      } else {
        a[ii] = M_PI/ell;
      }
    }

    h[0] = A * cos(a[0]*x[0]) * sin(a[1]*x[1]) * sin(a[2]*x[2]);
    h[1] = B * sin(a[0]*x[0]) * cos(a[1]*x[1]) * sin(a[2]*x[2]);
    h[2] = C * sin(a[0]*x[0]) * sin(a[1]*x[1]) * cos(a[2]*x[2]);
  } else {
    h[0] = 0.0;
    h[1] = 0.0;
    h[2] = 0.0;
  }

#ifndef NDEBUG
  const bool found_nan = ((h[0]!=h[0]) || (h[1]!=h[1]) || (h[2]!=h[2]) );
  if (found_nan) {
    std::cout << "Found a NaN in forcing!" << std::endl;
    std::cout << "xyz   = " << x[0] << " " << x[1] << " " << x[2] << std::endl;
    std::cout << "a     = " << a[0] << " " << a[1] << " " << a[2] << std::endl;
    std::cout << "dwall = " << dwall << std::endl;
    std::cout << "Lsgs  = " << Lsgs << std::endl;
    std::cout << "D     = " << D[0] << " " << D[1] << " " << D[2] << std::endl;
    std::cout << "h     = " << h[0] << " " << h[1] << " " << h[2] << std::endl;
    SU2_MPI::Error("Found a NaN in forcing.", CURRENT_FUNCTION);
  }
#endif

}

inline void CHybridForcingTG0::SetAxiTGField(
                const su2double* x, const su2double Lsgs,
                const su2double* D,
                const su2double dwall, su2double* h) const {

  // Convert incoming coords and lengths to cylindrical coords...
  // In these vectors, 0 corresponds to x, 1 corresponds to r, and 2
  // corresponds to theta
  // NB: Assume that user-specified D comes in in x,r,theta
  su2double r[3];
  r[0] = x[0]; r[1] = sqrt(x[1]*x[1] + x[2]*x[2]); r[2] = atan(x[2]/x[1]);

  su2double Rsgs[3]; // assume Lsgs in each direction
  Rsgs[0] = Lsgs;
  Rsgs[1] = Lsgs; //cos(r[2])*Lsgs + sin(r[2])*Lsgs;
  Rsgs[2] = Lsgs/r[1]; //(-sin(r[2])*Lsgs + cos(r[2])*Lsgs)/r[1];


  // Set forcing velocity field in x,r,theta coords

  //const su2double A = 1./3., B = 2./3., C = -1.0;
  const su2double A = 1.0, B = -1.0;
  su2double a[3];

  for (unsigned int ii=0; ii<3; ii++) {
    //const su2double ell = std::min(Rsgs[ii], dwall);
    const su2double ell = Rsgs[ii];

    if (D[ii] > 0.0) {
      const su2double denom = round(D[ii]/std::min(ell, D[ii]));
      a[ii] = M_PI/(D[ii]/denom);
    } else {
      a[ii] = M_PI/ell;
    }
  }

  su2double theta_deg = r[2]*180.0/M_PI;

  // su2double htmp[3];
  // if ( (theta_deg > 2.0) && (theta_deg < 13.0) ){
  //   htmp[0] = A * cos(a[0]*r[0]) * sin(a[1]*r[1]) * sin(a[2]*r[2]);
  //   htmp[1] = B * sin(a[0]*r[0]) * cos(a[1]*r[1]) * sin(a[2]*r[2]);
  //   htmp[2] = (B/a[2]) * sin(a[0]*r[0]) * sin(a[1]*r[1]) * cos(a[2]*r[2]);
  // } else {
  //   htmp[0] = 0.0;
  //   htmp[1] = 0.0;
  //   htmp[2] = 0.0;
  // }

  su2double htmp[3];
  htmp[0] = A * cos(a[0]*r[0]) * sin(a[1]*r[1]) * sin(a[2]*r[2]);
  htmp[1] = B * sin(a[0]*r[0]) * cos(a[1]*r[1]) * sin(a[2]*r[2]);
  htmp[2] = (B/a[2]) * sin(a[0]*r[0]) * sin(a[1]*r[1]) * cos(a[2]*r[2]);

  h[0] = htmp[0];
  h[1] = htmp[1]*cos(r[2]) - htmp[2]*sin(r[2]);
  h[2] = htmp[1]*sin(r[2]) + htmp[2]*cos(r[2]);

  bool found_nan = ((h[0]!=h[0]) || (h[1]!=h[1]) || (h[2]!=h[2]) );
  if (found_nan) {
    std::cout << "WTF!?! Found in forcing!" << std::endl;
    std::cout << "xyz = " << x[0] << " " << x[1] << " " << x[2] << std::endl;
    std::cout << "xrt = " << r[0] << " " << r[1] << " " << r[2] << std::endl;
    std::cout << "a   = " << a[0] << " " << a[1] << " " << a[2] << std::endl;
    std::cout << "htmp= " << htmp[0] << " " << htmp[1] << " " << htmp[2] << std::endl;
    std::cout << "h   = " << h[0] << " " << htmp[1] << " " << htmp[2] << std::endl;
    SU2_MPI::Error("Found a NaN in forcing.", CURRENT_FUNCTION);
  }

}


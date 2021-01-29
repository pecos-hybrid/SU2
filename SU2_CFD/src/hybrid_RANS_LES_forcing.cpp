/*!
 * \file hybrid_RANS_LES_forcing.cpp
 * \brief Implementation of the hybrid RANS/LES forcing terms
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

#include "../include/hybrid_RANS_LES_forcing.hpp"
#include "../include/solver_structure.hpp"
#include "../include/variable_structure_v2f.hpp"

CHybridForcingAbstractBase::CHybridForcingAbstractBase(
                              const unsigned short nDim,
                              const unsigned long nPoint,
                              const unsigned long nPointDomain)
  : nDim(nDim), nVar(nDim), nVarGrad(nDim),
    nPoint(nPoint), nPointDomain(nPointDomain) {
}

CHybridForcingAbstractBase::~CHybridForcingAbstractBase() {
}

CHybridForcingTG0::CHybridForcingTG0(const unsigned short nDim,
                               const unsigned long nPoint,
                               const unsigned long nPointDomain)
    : CHybridForcingAbstractBase(nDim, nPoint, nPointDomain),
      forcing_scale(4), C_F(8) {

  node = new su2double*[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new su2double[nVar];
  }

}


CHybridForcingTG0::CHybridForcingTG0(CGeometry* geometry, const CConfig* config)
    : CHybridForcingAbstractBase(geometry->GetnDim(),
                                 geometry->GetnPoint(),
                                 geometry->GetnPointDomain()),
      forcing_scale(config->GetHybrid_Forcing_Vortex_Length()),
      C_F(config->GetHybrid_Forcing_Strength()) {

  node = new su2double*[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new su2double[nVar];
  }

}

CHybridForcingTG0::~CHybridForcingTG0() {

  if (node != NULL) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      delete [] node[iPoint];
    }
    delete [] node;
  }
}

void CHybridForcingTG0::ComputeForcingField(CSolver** solver, CGeometry *geometry,
                                            CConfig *config) {

  const su2double TKE_MIN = 1.0E-8;
  const su2double V2_MIN = 2.0/3*TKE_MIN;

  const unsigned short kind_time_marching = config->GetUnsteady_Simulation();

  assert(kind_time_marching == TIME_STEPPING   ||
         kind_time_marching == DT_STEPPING_1ST ||
         kind_time_marching == DT_STEPPING_2ND );

  const su2double time = config->GetCurrent_UnstTimeND();
  assert(time >= 0);
  /*--- Timestep is used to check if forcing is physical ---*/
  const su2double delta_t = solver[FLOW_SOL]->node[0]->GetDelta_Time();
  assert(delta_t > 0);

  /*--- Allocate some scratch arrays to avoid continual reallocation ---*/
  assert(nDim == 3);
  su2double h[3]; // Initial TG vortex field.
  su2double Lsgs; // SGS length scale
  su2double Tsgs; // SGS time scale
  su2double x[3]; // Position
  su2double dwall; // distance to the wall
  su2double uf[3];

  // Domain lengths for periodic directions
  su2double *D = config->GetHybrid_Forcing_Periodic_Length();

  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    /*--- This velocity retrieval works only for compressible ---*/
    assert(config->GetKind_Regime() == COMPRESSIBLE);
    const su2double* prim_vars = solver[FLOW_SOL]->node[iPoint]->GetPrimitive();
    const su2double* prim_mean =
      solver[FLOW_SOL]->average_node[iPoint]->GetPrimitive();

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      /*--- Copy the values (and not the pointer), since we're changing them. ---*/
      x[iDim] = geometry->node[iPoint]->GetCoord(iDim);
      x[iDim] = TransformCoords(x[iDim], prim_mean[iDim+1], time);
    }

    /*--- Setup the TG vortex ---*/

    // total k, epsilon, and v2 from model
    const su2double ktot = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
    const su2double aniso_ratio = solver[TURB_SOL]->node[iPoint]->GetAnisoRatio();
    const su2double v2 = max(TWO3*ktot*aniso_ratio, V2_MIN);

    // ratio of modeled to total TKE
    const su2double beta = solver[FLOW_SOL]->average_node[iPoint]->GetKineticEnergyRatio();

    // average of r_M
    const su2double resolution_adequacy =
      solver[FLOW_SOL]->average_node[iPoint]->GetResolutionAdequacy();

    const su2double L_typical = solver[TURB_SOL]->node[iPoint]->GetTypicalLengthscale();
    const su2double L_kol = solver[TURB_SOL]->node[iPoint]->GetKolLengthscale();
    const su2double T_typical = solver[TURB_SOL]->node[iPoint]->GetTypicalTimescale();
    const su2double T_kol = solver[TURB_SOL]->node[iPoint]->GetKolTimescale();
    assert(L_typical >= 0);
    assert(L_kol > 0);
    assert(T_typical >= 0);
    assert(T_kol > 0);

    Lsgs = max(forcing_scale * pow(beta, 1.5) * L_typical, L_kol);
    Lsgs = max(Lsgs, L_kol);
    Tsgs = beta * T_typical;
    Tsgs = max(Tsgs, T_kol);

    // Get dwall
    dwall = geometry->node[iPoint]->GetWall_Distance();

    // Period of the forcing structure in each direction
    su2double periods[3];
    if (config->isHybrid_Forced_Axi()) {
      // Angular periodic version 
      this->SetAxiTGField(x, Lsgs, D, dwall, h, periods);
    } else {
      // Compute TG velocity at this point
      this->SetTGField(x, Lsgs, D, dwall, h, periods);
    }

    const su2double Ftar = this->GetTargetProduction(v2, Tsgs, beta);

    // Compute PFtest
    su2double PFtest = 0.0;
    for (unsigned short iDim=0; iDim<nDim; iDim++) {
      uf[iDim] = prim_vars[iDim+1] - prim_mean[iDim+1];
      PFtest += uf[iDim]*h[iDim];
    }
    PFtest *= Ftar;

    const su2double beta_kol = solver[TURB_SOL]->node[iPoint]->GetKolKineticEnergyRatio();

    su2double eta = this->ComputeScalingFactor(resolution_adequacy,
                                               beta, beta_kol, PFtest);

    // XXX: Turn off forcing around the shock and the separation region.
    const bool zonal_forcing = config->GetHybrid_Forcing_Cutoff();
    su2double forcing_clipping = 0.0;
    if (zonal_forcing) {
      // We modified "x" above, so retrieve original (global) value of x
      const su2double x_global = geometry->node[iPoint]->GetCoord(0);
      forcing_clipping = 0.5 - 0.5*tanh(10*(x_global - 0.35));
      eta *= forcing_clipping;
    }

    /*--- Store eta*h so we can compute the derivatives ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      node[iPoint][iDim] = Ftar*eta*h[iDim];
    }

    /*--- Save the forcing for output ---*/

    solver[FLOW_SOL]->node[iPoint]->SetForcingVector(node[iPoint]);
    if (config->GetWrt_Hybrid_Forcing()) {
      solver[FLOW_SOL]->node[iPoint]->SetForcingFactor(eta);
      solver[FLOW_SOL]->node[iPoint]->SetForcingClipping(forcing_clipping);
      solver[FLOW_SOL]->node[iPoint]->SetForcingLength(periods);
      solver[FLOW_SOL]->node[iPoint]->SetForcingProduction(PFtest);
      solver[FLOW_SOL]->node[iPoint]->SetForcingStructure(h);
    }

  } // end loop over points
}

su2double CHybridForcingTG0::ComputeScalingFactor(
                     const su2double resolution_adequacy,
                     const su2double beta,
                     const su2double beta_kol,
                     const su2double PFtest) const {

  su2double eta = 0.0;

  if (PFtest >= 0.0) {
    const su2double F_r = -tanh(1.0 - pow(min(resolution_adequacy, 1.0), -0.5));
    const su2double beta_hat = (1.0 - beta)/max(1.0 - beta_kol, 0.01) - 1;
    const su2double Dlim = F_r * (tanh(10*beta_hat)+1);
    eta = F_r - Dlim;
  }

  // Check for NaNs in debug mode.
  assert(eta == eta);

  return eta;
}

void CHybridForcingTG0::SetTGField(const su2double* x,
                                   const su2double Lsgs,
                                   const su2double* D,
                                   const su2double dwall,
                                   su2double* h,
                                   su2double* lmod) const {

  constexpr su2double A = 1./3., B = 2./3., C = -1.0;
  su2double a[3];

  /*--- Using the wall-distance as a max is overly conservative. The
     attached-eddy-hypothesis only says that the largest eddy diameter
     is proportional to the wall-normal distance; the constant of
     proportionality is undefined.  Using the mixing length paradigm,
     the constant of proportionality is the von Karman constant.---*/
  constexpr su2double kappa = 0.41;
  /*--- Relationship between k^{3/2}/epsilon and integral
   lengthscale for high Reynolds numbers. See Pope pg 244 ---*/
  constexpr su2double L_11_over_L = 0.43;
  if (dwall > 0 && Lsgs > 0) {
    for (unsigned int ii=0; ii<3; ii++) {
      const su2double ell = min(L_11_over_L*Lsgs, kappa*dwall);

      if (D[ii] > 0.0) {
        const su2double denom = round(D[ii]/min(ell, D[ii]));
        a[ii] = 2*M_PI/D[ii]*denom;
      } else {
        a[ii] = 2*M_PI/ell;
      }
      lmod[ii] = a[ii]/(2*M_PI);
    }

    h[0] = A * cos(a[0]*x[0]) * sin(a[1]*x[1]) * sin(a[2]*x[2]);
    h[1] = B * sin(a[0]*x[0]) * cos(a[1]*x[1]) * sin(a[2]*x[2]);
    h[2] = C * sin(a[0]*x[0]) * sin(a[1]*x[1]) * cos(a[2]*x[2]);
  } else {
    h[0] = 0.0;
    h[1] = 0.0;
    h[2] = 0.0;
    lmod[0] = 0;
    lmod[1] = 0;
    lmod[2] = 0;
  }

#ifndef NDEBUG
  const bool found_nan = ((h[0]!=h[0]) || (h[1]!=h[1]) || (h[2]!=h[2]) );
  if (found_nan) {
    std::cout << "Bad value found in forcing!" << std::endl;
    std::cout << "xyz   = " << x[0] << " " << x[1] << " " << x[2] << std::endl;
    std::cout << "dwall = " << dwall << std::endl;
    std::cout << "Lsgs  = " << Lsgs << std::endl;
    std::cout << "a     = " << a[0] << " " << a[1] << " " << a[2] << std::endl;
    std::cout << "D     = " << D[0] << " " << D[1] << " " << D[2] << std::endl;
    std::cout << "h     = " << h[0] << " " << h[1] << " " << h[2] << std::endl;
    SU2_MPI::Error("Found a NaN in forcing.", CURRENT_FUNCTION);
  }
#endif

}

void CHybridForcingTG0::SetAxiTGField(const su2double* x,
                                      const su2double Lsgs,
                                      const su2double* D,
                                      const su2double dwall,
                                      su2double* h,
                                      su2double* lmod) const {

  // Convert incoming coords and lengths to cylindrical coords...
  // In these vectors, 0 corresponds to x, 1 corresponds to r, and 2
  // corresponds to theta
  // NB: Assume that user-specified D comes in in x,r,theta
  su2double r[3];
  r[0] = x[0]; r[1] = sqrt(x[1]*x[1] + x[2]*x[2]); r[2] = atan2(x[2], x[1]);

  su2double Rsgs[3]; // assume Lsgs in each direction
  Rsgs[0] = Lsgs;
  Rsgs[1] = Lsgs; //cos(r[2])*Lsgs + sin(r[2])*Lsgs;
  // Angle (in radians) of a sector with arc-length Lsgs is Lsgs/r
  Rsgs[2] = Lsgs/r[1]; //(-sin(r[2])*Lsgs + cos(r[2])*Lsgs)/r[1];

  // Express wall-distance in equivalent radians
  su2double wall_dist[3] = {dwall, dwall, dwall/r[1]};

  // Set forcing velocity field in x,r,theta coords

  su2double a[3];
  /*--- Using the wall-distance as a max is overly permissive. The
     attached-eddy-hypothesis only says that the large-eddy diameter
     is proportional to the wall-normal distance; the constant of
     proportionality is undefined.  Using the mixing length paradigm,
     the constant of proportionality is the von Karman constant.---*/
  constexpr su2double kappa = 0.41;
  /*--- Relationship between k^{3/2}/epsilon and integral
   lengthscale for high Reynolds numbers. See Pope pg 244 ---*/
  constexpr su2double L_11_over_L = 0.43;
  if (dwall > 0 && Lsgs > 0) {
    for (unsigned int ii=0; ii<3; ii++) {
      const su2double ell = std::min(L_11_over_L*Rsgs[ii], kappa*wall_dist[ii]);

      if (D[ii] > 0.0) {
        const su2double denom = round(D[ii]/std::min(ell, D[ii]));
        a[ii] = 2*M_PI/(D[ii]/denom);
      } else {
        a[ii] = 2*M_PI/ell;
      }
      lmod[ii] = a[ii]/(2*M_PI);
    }

    const su2double A = (abs(a[1]) > 1E-16) ? 1.0 : 0;
    const su2double B = (abs(a[1]) > 1E-16) ? -A*a[0]/a[1] : 0;
    // Forcing functions in steamwise, wall-normal, and radial
    su2double htmp[3];
    htmp[0] = A * cos(a[0]*r[0]) * sin(a[1]*r[1]) * sin(a[2]*r[2]);
    htmp[1] = B * sin(a[0]*r[0]) * cos(a[1]*r[1]) * sin(a[2]*r[2]);
    /*--- Extra cosine in the next line is not accidental. To satisfy
          incompressibility, we need:
            htmp[1] + \partial htmp[2] / \partial r[2] = 0 
          See Pope Eq. 4.30, on page 110 ---*/
    htmp[2] = (B/a[2]) * sin(a[0]*r[0]) * cos(a[1]*r[1]) * cos(a[2]*r[2]);

    h[0] = htmp[0];
    h[1] = htmp[1]*cos(r[2]) - htmp[2]*sin(r[2]);
    h[2] = htmp[1]*sin(r[2]) + htmp[2]*cos(r[2]);
    bool found_nan = ((h[0]!=h[0]) || (h[1]!=h[1]) || (h[2]!=h[2]) );
    if (found_nan) {
      std::cout << "Bad value found in forcing!" << "\n";
      std::cout << "xyz = " << x[0] << " " << x[1] << " " << x[2] << "\n";
      std::cout << "xrt = " << r[0] << " " << r[1] << " " << r[2] << "\n";
      std::cout << "Rsgs = " << Rsgs[0] << " " << Rsgs[1] << " " << Rsgs[2] << "\n";
      std::cout << "wall_dist = " << wall_dist[0] << " " << wall_dist[1] << " " << wall_dist[2] << "\n";
      std::cout << "a   = " << a[0] << " " << a[1] << " " << a[2] << "\n";
      std::cout << "htmp= " << htmp[0] << " " << htmp[1] << " " << htmp[2] << "\n";
      std::cout << "h   = " << h[0] << " " << h[1] << " " << h[2] << std::endl;
      SU2_MPI::Error("Found a NaN in forcing.", CURRENT_FUNCTION);
    }
  } else {
    h[0] = 0;
    h[1] = 0;
    h[2] = 0;
    lmod[0] = 0;
    lmod[1] = 0;
    lmod[2] = 0;
  }


}

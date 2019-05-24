/*!
 * \file numerics_direct_turbulent.cpp
 * \brief This file contains all the convective term discretization.
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

#include "../../Common/include/config_structure.hpp"
#include "../include/numerics_structure.hpp"
#include "../include/numerics_structure_v2f.hpp"
#include <limits>

CUpwSca_TurbKE::CUpwSca_TurbKE(unsigned short val_nDim,
                               unsigned short val_nVar,
                               CConfig *config)
  :
  CUpwScalar(val_nDim, val_nVar, config) {
}

CUpwSca_TurbKE::~CUpwSca_TurbKE(void) {
}

void CUpwSca_TurbKE::ExtraADPreaccIn() {
  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+2);
    AD::SetPreaccIn(V_j, nDim+2);
  }
  else {
    AD::SetPreaccIn(V_i, nDim+3);
    AD::SetPreaccIn(V_j, nDim+3);
  }
}

void CUpwSca_TurbKE::FinishResidualCalc(su2double *val_residual,
                                        su2double **val_Jacobian_i,
                                        su2double **val_Jacobian_j,
                                        CConfig *config) {

  val_residual[0] = a0*Density_i*TurbVar_i[0]+a1*Density_j*TurbVar_j[0];
  val_residual[1] = a0*Density_i*TurbVar_i[1]+a1*Density_j*TurbVar_j[1];
  val_residual[2] = a0*Density_i*TurbVar_i[2]+a1*Density_j*TurbVar_j[2];
  val_residual[3] = 0.0; // no convection in f scalar


  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_i[0][1] = 0.0;
    val_Jacobian_i[0][2] = 0.0;
    val_Jacobian_i[0][3] = 0.0;

    val_Jacobian_i[1][0] = 0.0;
    val_Jacobian_i[1][1] = a0;
    val_Jacobian_i[1][2] = 0.0;
    val_Jacobian_i[1][3] = 0.0;

    val_Jacobian_i[2][0] = 0.0;
    val_Jacobian_i[2][1] = 0.0;
    val_Jacobian_i[2][2] = a0;
    val_Jacobian_i[2][3] = 0.0;

    val_Jacobian_i[3][0] = 0.0;
    val_Jacobian_i[3][1] = 0.0;
    val_Jacobian_i[3][2] = 0.0;
    val_Jacobian_i[3][3] = 0.0;


    val_Jacobian_j[0][0] = a1;
    val_Jacobian_j[0][1] = 0.0;
    val_Jacobian_j[0][2] = 0.0;
    val_Jacobian_j[0][3] = 0.0;

    val_Jacobian_j[1][0] = 0.0;
    val_Jacobian_j[1][1] = a1;
    val_Jacobian_j[1][2] = 0.0;
    val_Jacobian_j[1][3] = 0.0;

    val_Jacobian_j[2][0] = 0.0;
    val_Jacobian_j[2][1] = 0.0;
    val_Jacobian_j[2][2] = a1;
    val_Jacobian_j[2][3] = 0.0;

    val_Jacobian_j[3][0] = 0.0;
    val_Jacobian_j[3][1] = 0.0;
    val_Jacobian_j[3][2] = 0.0;
    val_Jacobian_j[3][3] = 0.0;
  }

}

CAvgGrad_TurbKE::CAvgGrad_TurbKE(unsigned short val_nDim,
                                 unsigned short val_nVar,
                                 su2double *constants,
                                 bool correct_grad,
                                 CConfig *config)
  :
  CAvgGrad_Scalar(val_nDim, val_nVar, correct_grad, config) {

  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];
}

CAvgGrad_TurbKE::~CAvgGrad_TurbKE(void) {
}

void CAvgGrad_TurbKE::ExtraADPreaccIn() {
  AD::SetPreaccIn(Volume);
}

void CAvgGrad_TurbKE::FinishResidualCalc(su2double *val_residual,
                                         su2double **Jacobian_i,
                                         su2double **Jacobian_j,
                                         CConfig *config) {

  su2double sigma_kine_i, sigma_kine_j, sigma_epsi_i, sigma_epsi_j;
  su2double sigma_zeta_i, sigma_zeta_j;
  su2double diff_i_kine, diff_i_epsi, diff_j_kine, diff_j_epsi;
  su2double diff_i_zeta, diff_j_zeta;


  /*--- Compute the blended constant for the viscous terms ---*/
  sigma_kine_i = sigma_k;
  sigma_kine_j = sigma_k;
  sigma_epsi_i = sigma_e;
  sigma_epsi_j = sigma_e;
  sigma_zeta_i = sigma_z;
  sigma_zeta_j = sigma_z;

  /*--- Compute mean effective viscosity ---*/
  diff_i_kine = Laminar_Viscosity_i + sigma_kine_i*Eddy_Viscosity_i;
  diff_j_kine = Laminar_Viscosity_j + sigma_kine_j*Eddy_Viscosity_j;
  diff_i_epsi = Laminar_Viscosity_i + sigma_epsi_i*Eddy_Viscosity_i;
  diff_j_epsi = Laminar_Viscosity_j + sigma_epsi_j*Eddy_Viscosity_j;
  diff_i_zeta = Laminar_Viscosity_i + sigma_zeta_i*Eddy_Viscosity_i;
  diff_j_zeta = Laminar_Viscosity_j + sigma_zeta_j*Eddy_Viscosity_j;

  // Could instead use weighted average!
  diff_kine = 0.5*(diff_i_kine + diff_j_kine);
  diff_epsi = 0.5*(diff_i_epsi + diff_j_epsi);
  diff_zeta = 0.5*(diff_i_zeta + diff_j_zeta);
  diff_f = 1.0;

  /*--- Compute vector going from iPoint to jPoint ---*/
  su2double s_mag=sqrt(dist_ij_2);

  su2double n_mag=0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    n_mag += Normal[iDim]*Normal[iDim];
  }

  n_mag = sqrt(n_mag);

  val_residual[0] = diff_kine*Proj_Mean_GradTurbVar[0];
  val_residual[1] = diff_epsi*Proj_Mean_GradTurbVar[1];
  val_residual[2] = diff_zeta*Proj_Mean_GradTurbVar[2];
  val_residual[3] = diff_f   *Proj_Mean_GradTurbVar[3];

  /*--- For Jacobians ->
    Use of TSL approx. to compute derivatives of the gradients ---*/
  if (implicit) {

    Jacobian_i[0][0] = -diff_kine*proj_vector_ij/Density_i;
    Jacobian_i[0][1] = 0.0;
    Jacobian_i[0][2] = 0.0;
    Jacobian_i[0][3] = 0.0;

    Jacobian_i[1][0] = 0.0;
    Jacobian_i[1][1] = -diff_epsi*proj_vector_ij/Density_i;
    Jacobian_i[1][2] = 0.0;
    Jacobian_i[1][3] = 0.0;

    Jacobian_i[2][0] = 0.0;
    Jacobian_i[2][1] = 0.0;
    Jacobian_i[2][2] = -diff_zeta*proj_vector_ij/Density_i;
    Jacobian_i[2][3] = 0.0;

    Jacobian_i[3][0] = 0.0;
    Jacobian_i[3][1] = 0.0;
    Jacobian_i[3][2] = 0.0;
    Jacobian_i[3][3] = -diff_f*proj_vector_ij;
    // if (correct_gradient) {
    //   Jacobian_i[3][3] -= 0.5 * diff_f * n_mag/Volume;
    //   Jacobian_i[3][3] += 0.5 * diff_f * proj_vector_ij*s_mag/Volume;
    // }

    Jacobian_j[0][0] = diff_kine*proj_vector_ij/Density_j;
    Jacobian_j[0][1] = 0.0;
    Jacobian_j[0][2] = 0.0;
    Jacobian_j[0][3] = 0.0;

    Jacobian_j[1][0] = 0.0;
    Jacobian_j[1][1] = diff_epsi*proj_vector_ij/Density_j;
    Jacobian_j[1][2] = 0.0;
    Jacobian_j[1][3] = 0.0;

    Jacobian_j[2][0] = 0.0;
    Jacobian_j[2][1] = 0.0;
    Jacobian_j[2][2] = diff_zeta*proj_vector_ij/Density_j;
    Jacobian_j[2][3] = 0.0;

    Jacobian_j[3][0] = 0.0;
    Jacobian_j[3][1] = 0.0;
    Jacobian_j[3][2] = 0.0;
    Jacobian_j[3][3] = diff_f*proj_vector_ij;
    // if (correct_gradient) {
    //   Jacobian_j[3][3] += 0.5 * diff_f * n_mag/Volume;
    //   Jacobian_j[3][3] -= 0.5 * diff_f * proj_vector_ij*s_mag/Volume;
    // }
  }

}


CSourcePieceWise_TurbKE::CSourcePieceWise_TurbKE(unsigned short val_nDim,
                                                 unsigned short val_nVar,
                                                 su2double *constants,
                                                 CConfig *config)
  :
  CNumerics(val_nDim, val_nVar, config) {

  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  /*--- Closure constants ---*/

  /* zeta f */
  C_mu    = constants[0];
  sigma_k = constants[1];
  sigma_e = constants[2];
  sigma_z = constants[3];
  C_e1o   = constants[4];
  C_e2    = constants[5];
  C_1     = constants[6];
  C_2p    = constants[7];
  C_T     = constants[8];
  C_L     = constants[9];
  C_eta   = constants[10];

}

CSourcePieceWise_TurbKE::~CSourcePieceWise_TurbKE(void) { }

void CSourcePieceWise_TurbKE::ComputeResidual(su2double *val_residual,
                                              su2double **val_Jacobian_i,
                                              su2double **val_Jacobian_j,
                                              CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+7);
  AD::SetPreaccIn(StrainMag_i);
  AD::SetPreaccIn(TurbVar_i, nVar);
  AD::SetPreaccIn(TurbVar_Grad_i, nVar, nDim);
  AD::SetPreaccIn(Volume); AD::SetPreaccIn(dist_i);
  AD::SetPreaccIn(PrimVar_Grad_i, nDim+1, nDim);

  // Pull variables out of V_i
  if (incompressible) {
    Density_i = V_i[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];
  }
  else {
    Density_i = V_i[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];
  }

  // for readability...
  const su2double tke = TurbVar_i[0];
  const su2double tdr = TurbVar_i[1];
  const su2double v20 = TurbVar_i[2];
  const su2double f   = TurbVar_i[3];

  // clip values to avoid non-physical quantities...

  su2double scale = 1.0e-14;
  su2double L_inf = config->GetLength_Reynolds();
  su2double* VelInf = config->GetVelocity_FreeStreamND();
  su2double VelMag = 0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    VelMag += VelInf[iDim]*VelInf[iDim];
  VelMag = sqrt(VelMag);

  const su2double tke_lim = max(tke, scale*VelMag*VelMag);
  const su2double tdr_lim = max(tdr, scale*VelMag*VelMag*VelMag/L_inf);

  // make sure v2 is well-behaved
  su2double zeta = max(v20/tke_lim, scale);
  // Extra max(..., 0) necessary in case v20 and tke are negative
  const su2double v2 = max(max(v20, zeta*tke), 0.0);

  // Grab other quantities for convenience/readability
  const su2double rho = Density_i;
  const su2double muT = Eddy_Viscosity_i;
  const su2double S   = StrainMag_i; //*sqrt(2.0) already included
  const su2double Vol = Volume;

  su2double diverg = 0.0;
  for (unsigned int iDim = 0; iDim < nDim; iDim++)
    diverg += PrimVar_Grad_i[iDim+1][iDim];

  /*--- TurbT and TurbL are passed in from the Solver class ---*/

  const su2double T1 = tke_lim / tdr_lim;
  const su2double Lsq = (C_L*TurbL)*(C_L*TurbL); // C_L isn't used till here

  //--- v2-f ---//

  // 4 equations.  For each equation, we identify production and
  // dissipation terms.  This is somewhat artificial for f.  The
  // Jacobian is abused to help keep k and epsilon positive.


  // TKE equation...
  su2double Pk_rk, Pk_re, Pk_rv2;
  su2double Dk, Dk_rk, Dk_re, Dk_rv2;

  //... production
  // NB: we ignore the Jacobian of production here

  //SGSProduction     = muT*S*S - 2.0/3.0*rho*tke*diverg;
  SGSProduction = muT*S*S;
  if (config->GetBoolDivU_inTKEProduction()) {
    SGSProduction -= (2.0/3.0)*rho*tke*diverg;
  }

  /*--- If using a hybrid method, include resolved production ---*/

  if (config->GetKind_HybridRANSLES() == MODEL_SPLIT) {
    /*--- Limit alpha to protect from imbalance in k_model vs k_resolved. ---*/
    if (KineticEnergyRatio >= 0  && KineticEnergyRatio < 1) {
      //SGSProduction *= KineticEnergyRatio;
      const su2double alpha = KineticEnergyRatio;
      const su2double alpha_fac = alpha*(2.0 - alpha);
      SGSProduction = alpha_fac*muT*S*S;
      if (config->GetBoolDivU_inTKEProduction()) {
	SGSProduction -= 2.0/3.0*rho*alpha*tke*diverg;
      }

      if (config->GetUse_Resolved_Turb_Stress()) {
        su2double Pk_resolved = 0;
        for (unsigned short iDim = 0; iDim < nDim; iDim++) {
          for (unsigned short jDim = 0; jDim < nDim; jDim++) {
            Pk_resolved += PrimVar_Grad_i[iDim+1][jDim] * ResolvedTurbStress[iDim][jDim];
          }
        }
        Pk = SGSProduction + Pk_resolved;
      }
      /*--- If using production instead of resolved turb. stress,
       * the production (`Pk`) has been set externally.  So we only need to
       * make sure that SGSProduction has been properly calculated, since
       * it will be pulled back out and passed to the averaging routine. ---*/
    } else {
      /*--- Don't allow resolved turb. stress to be added if alpha is not
       * valid, in case the resolved flow data are bad. ---*/
      Pk = SGSProduction;
    }
  } else {
    Pk = SGSProduction;
  }

  Pk_rk  = 0.0;
  Pk_re  = 0.0;
  Pk_rv2 = 0.0;

  //... dissipation
  Dk     = rho*tke/T1;

  Dk_rk  = 1.0/T1;
  Dk_re  = 0.0;
  Dk_rv2 = 0.0;


  // Dissipation equation...
  su2double Pe, Pe_rk, Pe_re, Pe_rv2;
  su2double De, De_rk, De_re, De_rv2;

  // NB: C_e1 depends on tke and v2 in v2-f
  const su2double C_e1 = C_e1o*(1.0+0.045*sqrt(1.0/zeta));

  // ... production
  Pe = C_e1*Pk/TurbT;

  Pe_rk  = 0.0;
  Pe_re  = 0.0;
  Pe_rv2 = 0.0;

  // ... dissipation
  De = C_e2*rho*tdr/TurbT;

  De_rk  = 0.0;
  De_re  = C_e2/TurbT;
  De_rv2 = 0.0;


  // v2 equation...
  su2double Pv2, Pv2_rk, Pv2_re, Pv2_rv2, Pv2_f;
  su2double Dv2, Dv2_rk, Dv2_re, Dv2_rv2, Dv2_f;

  // ... production
  // Limit production of v2 based on max zeta = 2/3
  Pv2 = rho * min( tke*f, 2.0*Pk/3.0/rho + 5.0*v2/T1 );

  Pv2_rk  = 0.0;
  Pv2_re  = 0.0;
  Pv2_rv2 = 0.0;
  Pv2_f   = 0.0;

  // ... dissipation
  Dv2     =  6.0*rho*v2/T1;

  Dv2_rk  = 0.0;
  Dv2_re  = 0.0;
  Dv2_rv2 = 6.0/T1;
  Dv2_f   = 0.0;


  // f equation...
  su2double Pf;
  su2double Df, Df_f;

  //... production
  const su2double C1m6 = C_1 - 6.0;
  const su2double ttC1m1 = (2.0/3.0)*(C_1 - 1.0);
  const su2double C_2f = C_2p;

  su2double Rf = 1.0/TurbT;
  if (config->GetBoolUse_v2f_Rf_mod()) {
    Rf = min(1.0/TurbT, S/(sqrt(2.0)*3.0));
  }

  //Pf = (C_2f*Pk/(rho*tke_lim) - (C1m6*zeta - ttC1m1)/TurbT) / Lsq;
  Pf = (C_2f*Pk/(rho*tke_lim) - Rf*(C1m6*zeta - ttC1m1)) / Lsq;

  // not keeping any derivatives of Pf

  //... dissipation
  Df = f/Lsq;

  Df_f = 1.0/Lsq;

  // check for nans
#ifndef NDEBUG
  bool found_nan = ((Pk!=Pk)         || (Dk!=Dk)         ||
                    (Pe!=Pe)         || (De!=De)         ||
                    (Pv2!=Pv2)       || (Dv2!=Dv2)       ||
                    (Pf!=Pf)         || (Df!=Df)         ||
                    (Pk_rk!=Pk_rk)   || (Pk_re!=Pk_re)   || (Pk_rv2!=Pk_rv2)   ||
                    (Pe_rk!=Pe_rk)   || (Pe_re!=Pe_re)   || (Pe_rv2!=Pe_rv2)   ||
                    (Pv2_rk!=Pv2_rk) || (Pv2_re!=Pv2_re) || (Pv2_rv2!=Pv2_rv2) ||
                    (Dk_rk!=Dk_rk)   || (Dk_re!=Dk_re)   || (Dk_rv2!=Dk_rv2)   ||
                    (De_rk!=De_rk)   || (De_re!=De_re)   || (De_rv2!=De_rv2)   ||
                    (Dv2_rk!=Dv2_rk) || (Dv2_re!=Dv2_re) || (Dv2_rv2!=Dv2_rv2) );


  if (found_nan) {
    std::cout << "WTF!?! Found a nan at x = " << Coord_i[0] << ", " << Coord_i[1] << std::endl;
    std::cout << "turb state = " << tke << ", " << tdr << ", " << v2 << ", " << f << std::endl;
    std::cout << "T          = " << TurbT  << ", C_e1 = " << C_e1 << std::endl;
    std::cout << "TKE eqn    = " << Pk << " - " << Dk << std::endl;
    std::cout << "TDR eqn    = " << Pe << " - " << De << std::endl;
    std::cout << "v2  eqn    = " << Pv2 << " - " << Dv2 << std::endl;
  }
#endif


  // form source term and Jacobian...

  // TKE
  val_residual[0]      = (Pk      - Dk     ) * Vol;

  val_Jacobian_i[0][0] = (Pk_rk   - Dk_rk  ) * Vol;
  val_Jacobian_i[0][1] = (Pk_re   - Dk_re  ) * Vol;
  val_Jacobian_i[0][2] = (Pk_rv2  - Dk_rv2 ) * Vol;
  val_Jacobian_i[0][3] = 0.0;

  // Dissipation
  val_residual[1]      = (Pe      - De     ) * Vol;

  val_Jacobian_i[1][0] = (Pe_rk   - De_rk  ) * Vol;
  val_Jacobian_i[1][1] = (Pe_re   - De_re  ) * Vol;
  val_Jacobian_i[1][2] = (Pe_rv2  - De_rv2 ) * Vol;
  val_Jacobian_i[1][3] = 0.0;

  // v2
  val_residual[2]      = (Pv2     - Dv2    ) * Vol;

  val_Jacobian_i[2][0] = (Pv2_rk  - Dv2_rk ) * Vol;
  val_Jacobian_i[2][1] = (Pv2_re  - Dv2_re ) * Vol;
  val_Jacobian_i[2][2] = (Pv2_rv2 - Dv2_rv2) * Vol;
  val_Jacobian_i[2][3] = (Pv2_f   - Dv2_f  ) * Vol;

  // f
  val_residual[3]      = (Pf      - Df     ) * Vol;

  val_Jacobian_i[3][0] = 0.0;
  val_Jacobian_i[3][1] = 0.0;
  val_Jacobian_i[3][2] = 0.0;
  val_Jacobian_i[3][3] = (        - Df_f   ) * Vol;


  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();
}

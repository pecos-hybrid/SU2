/*!
 * \file numerics_direct_mean.cpp
 * \brief This file contains the numerical methods for compressible flow.
 * \author F. Palacios, T. Economon
 * \version 6.2.0 "Falcon"
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
#include "../include/numerics_direct_mean_hybrid.hpp"

CAvgGrad_Hybrid::CAvgGrad_Hybrid(unsigned short val_nDim,
                             unsigned short val_nVar,
                             bool val_correct_grad,
                             CConfig *config)
    : CAvgGrad_Base(val_nDim, val_nVar, val_nDim+3, val_correct_grad, config) {

  conductivity = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    conductivity[iDim] = new su2double[nDim];
  }

  deviatoric = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    deviatoric[iDim] = new su2double[nDim];
  }

  Mean_PrimVar_Average = new su2double [nPrimVar];

  Mean_GradPrimVar_Fluct = new su2double* [nPrimVar];
  Mean_GradPrimVar_Average = new su2double* [nPrimVar];
  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    Mean_GradPrimVar_Average[iVar] = new su2double [nDim];
    Mean_GradPrimVar_Fluct[iVar] = new su2double [nDim];
  }

  Mean_Aniso_Eddy_Viscosity = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    Mean_Aniso_Eddy_Viscosity[iDim] = new su2double[nDim];
  }

  PrimVar_Average_i = NULL;
  PrimVar_Average_j = NULL;
  Aniso_Eddy_Viscosity_i = NULL;
  Aniso_Eddy_Viscosity_j = NULL;
  PrimVar_Grad_Average_i = NULL;
  PrimVar_Grad_Average_j = NULL;
}

CAvgGrad_Hybrid::~CAvgGrad_Hybrid(void) {

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] conductivity[iDim];
  }
  delete [] conductivity;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] deviatoric[iDim];
  }
  delete [] deviatoric;

  delete [] Mean_PrimVar_Average;

  for (unsigned short iVar = 0; iVar < nPrimVar; iVar++) {
    delete [] Mean_GradPrimVar_Fluct[iVar];
    delete [] Mean_GradPrimVar_Average[iVar];
  }
  delete [] Mean_GradPrimVar_Fluct;
  delete [] Mean_GradPrimVar_Average;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delete [] Mean_Aniso_Eddy_Viscosity[iDim];
  }
  delete [] Mean_Aniso_Eddy_Viscosity;

}

void CAvgGrad_Hybrid::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  // TODO: Re-include the AD stuff

  /*--- Check validity of input ---*/

  assert(PrimVar_Average_i != NULL);
  assert(PrimVar_Average_j != NULL);
  assert(Aniso_Eddy_Viscosity_i != NULL);
  assert(Aniso_Eddy_Viscosity_j != NULL);
  assert(PrimVar_Grad_Average_i != NULL);
  assert(PrimVar_Grad_Average_j != NULL);
  assert(beta_i == beta_i);  // beta_i is not NaN
  assert(beta_j == beta_j);  // beta_j is not NaN

  /*--- Normalized normal vector ---*/

  unsigned short iVar, jVar, iDim, jDim;

  Area = 0.0;
  for (iDim = 0; iDim < nDim; iDim++)
    Area += Normal[iDim]*Normal[iDim];
  Area = sqrt(Area);

  for (iDim = 0; iDim < nDim; iDim++)
    UnitNormal[iDim] = Normal[iDim]/Area;

  for (iVar = 0; iVar < nPrimVar; iVar++) {
    PrimVar_i[iVar] = V_i[iVar];
    PrimVar_j[iVar] = V_j[iVar];
    Mean_PrimVar[iVar] = 0.5*(PrimVar_i[iVar]+PrimVar_j[iVar]);
    Mean_PrimVar_Average[iVar] = 0.5*(PrimVar_Average_i[iVar]+PrimVar_Average_j[iVar]);
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
  }

  /*--- Laminar and Eddy viscosity ---*/

  Laminar_Viscosity_i = V_i[nDim+5]; Laminar_Viscosity_j = V_j[nDim+5];
  Eddy_Viscosity_i = V_i[nDim+6]; Eddy_Viscosity_j = V_j[nDim+6];

  /*--- Mean Viscosities and turbulent kinetic energy---*/

  Mean_Laminar_Viscosity = 0.5*(Laminar_Viscosity_i + Laminar_Viscosity_j);
  Mean_Eddy_Viscosity = 0.5*(Eddy_Viscosity_i + Eddy_Viscosity_j);
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      Mean_Aniso_Eddy_Viscosity[iDim][jDim] = 0.5*(Aniso_Eddy_Viscosity_i[iDim][jDim] + Aniso_Eddy_Viscosity_j[iDim][jDim]);
    }
  }
  Mean_turb_ke = 0.5*(turb_ke_i + turb_ke_j);

  const su2double Mean_Bulk_Viscosity = Mean_Laminar_Viscosity * config->GetBulkViscosityRatio();

  /*--- Limit beta to protect from imbalance in k_model vs k_resolved ---*/

  const su2double Mean_Beta = min(max(0.5*(beta_i + beta_j), 0.0), 1.0);

  /*--- Mean gradient approximation ---*/

  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Mean_GradPrimVar_Average[iVar][iDim] = 0.5*(PrimVar_Grad_Average_i[iVar][iDim] + PrimVar_Grad_Average_j[iVar][iDim]);
    }
  }
  for (iVar = 0; iVar < 1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradTurbVar[iVar][iDim] = 0.5*(TurbVar_Grad_i[iVar][iDim] + TurbVar_Grad_j[iVar][iDim]);
    }
  }

  /*--- Projection of the mean gradient in the direction of the edge ---*/

  if (correct_gradient && dist_ij_2 != 0.0) {
    CorrectGradient(Mean_GradPrimVar, PrimVar_i, PrimVar_j, Edge_Vector,
                    dist_ij_2, nDim+1);
    CorrectGradient(Mean_GradPrimVar_Average, PrimVar_Average_i,
                    PrimVar_Average_j, Edge_Vector, dist_ij_2, nDim+1);
  }

  /*--- Gradient of the fluctuating variables ---*/

  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar_Fluct[iVar][iDim] = Mean_GradPrimVar[iVar][iDim] - Mean_GradPrimVar_Average[iVar][iDim];
    }
  }

#ifndef NDEBUG
  // Check that fluctuations are zero if averaging is not to be performed
  const su2double time = config->GetCurrent_UnstTime();
  if (time < config->GetAveragingStartTime()) {
    for (iVar = 0; iVar < nDim+1; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        assert(abs(Mean_GradPrimVar_Fluct[iVar][iDim]) < EPS);
      }
    }
  }
#endif

  /*--- Get projected flux tensor ---*/

  SetLaminarStressTensor(Mean_GradPrimVar, Mean_Laminar_Viscosity, Mean_Bulk_Viscosity);
  AddTauSGS(Mean_PrimVar_Average, Mean_GradPrimVar_Average, Mean_Beta,
            Mean_turb_ke, Mean_Eddy_Viscosity);
  AddTauSGET(Mean_GradPrimVar_Fluct,
             Mean_Aniso_Eddy_Viscosity);

  SetLaminarHeatFlux(Mean_GradPrimVar, Mean_Laminar_Viscosity);
  AddSGSHeatFlux(Mean_GradPrimVar_Average, Mean_Beta, Mean_Eddy_Viscosity);
  AddSGETHeatFlux(Mean_GradPrimVar_Fluct, Mean_Aniso_Eddy_Viscosity);

  if (config->GetUse_TKE_Diffusion()) {
    SetLaminar_TKE_Diffusion(Mean_GradTurbVar, Mean_Beta, Mean_Laminar_Viscosity);
    AddSGS_TKE_Diffusion(Mean_GradTurbVar, Mean_Beta, Mean_Eddy_Viscosity);
  }

  GetViscousProjFlux(Mean_PrimVar, Normal);

  /*--- Update viscous residual ---*/

  for (iVar = 0; iVar < nVar; iVar++)
    val_residual[iVar] = Proj_Flux_Tensor[iVar];

  /*--- Compute the implicit part ---*/

  if (implicit) {

    if (dist_ij_2 == 0.0) {
      for (iVar = 0; iVar < nVar; iVar++) {
        for (jVar = 0; jVar < nVar; jVar++) {
          val_Jacobian_i[iVar][jVar] = 0.0;
          val_Jacobian_j[iVar][jVar] = 0.0;
        }
      }
    } else {
      const su2double dist_ij = sqrt(dist_ij_2);
      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Bulk_Viscosity, 0,
                     dist_ij, UnitNormal);
      AddTauSGETJacobian(Mean_PrimVar, Mean_Aniso_Eddy_Viscosity,
			 dist_ij, UnitNormal);

      SetHeatFluxJacobian(Mean_PrimVar, Mean_Laminar_Viscosity,
                          0, dist_ij, UnitNormal);

      AddSGETHeatFluxJacobian(Mean_PrimVar, Mean_Aniso_Eddy_Viscosity,
			      dist_ij, UnitNormal);

      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor,
                         val_Jacobian_i, val_Jacobian_j);
    }

  }

}

void CAvgGrad_Hybrid::SetLaminarStressTensor(su2double **val_gradprimvar,
                                           const su2double val_laminar_viscosity,
                                           const su2double val_bulk_viscosity) {

  unsigned short iDim, jDim;

  su2double div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] = val_laminar_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
                        - TWO3*val_laminar_viscosity*div_vel*delta[iDim][jDim]
                        + val_bulk_viscosity * div_vel * delta[iDim][jDim];
    }
  }
}

void CAvgGrad_Hybrid::AddTauSGS(const su2double *val_primvar,
                               su2double **val_gradprimvar,
                               const su2double val_beta,
                               const su2double val_turb_ke,
                               const su2double val_eddy_viscosity) {

  assert(val_beta >= 0.0);
  assert(val_beta <= 1.0);
  unsigned short iDim, jDim;
  const su2double Density = val_primvar[nDim+2];

  const su2double val_alpha = pow(val_beta, 1.7);
  const su2double alpha_fac = val_alpha*(2.0 - val_alpha);
  const su2double mut_sgs = alpha_fac*val_eddy_viscosity;

  su2double div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
       tau[iDim][jDim] += mut_sgs*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
                         - TWO3*mut_sgs*div_vel*delta[iDim][jDim]
                         - TWO3*Density*val_beta*val_turb_ke*delta[iDim][jDim];
    }
  }
}

void CAvgGrad_Hybrid::AddTauSGET(su2double **val_gradprimvar,
                               su2double **val_eddy_viscosity) {
  su2double div_vel = 0.0;
  for (unsigned short iDim = 0 ; iDim < nDim; iDim++) {
    div_vel += val_gradprimvar[iDim+1][iDim];
  }
  for (unsigned short iDim =0; iDim < nDim; iDim++) {
    for (unsigned short jDim =0; jDim < nDim; jDim++) {
      deviatoric[iDim][jDim] = val_gradprimvar[iDim+1][jDim] - delta[iDim][jDim]*div_vel/3;
    }
  }
  for (unsigned short iDim =0; iDim < nDim; iDim++) {
    for (unsigned short jDim =0; jDim < nDim; jDim++) {
      for (unsigned short kDim =0; kDim < nDim; kDim++) {
        tau[iDim][jDim] += val_eddy_viscosity[iDim][kDim]*deviatoric[jDim][kDim] +
                           val_eddy_viscosity[jDim][kDim]*deviatoric[iDim][kDim];
        tau[iDim][iDim] -= TWO3*val_eddy_viscosity[jDim][kDim]*deviatoric[jDim][kDim];
      }
    }
  }
}

void CAvgGrad_Hybrid::AddTauSGETJacobian(su2double *val_Mean_PrimVar,
					 su2double **mu,
					 const su2double val_dist_ij,
					 const su2double *nvec) {

  /*--- QCR and wall functions are **not** accounted for here ---*/

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double xi = 1.0/(Density*val_dist_ij);

  su2double tauSGET_mom[3][3];

  // NB: Assumes normal and edge are aligned

  // Jacobian w.r.t. momentum
  for (unsigned short i = 0; i < nDim; i++) { // tau \cdot nvec component
    for (unsigned short j = 0; j < nDim; j++) { // momentum component
      tauSGET_mom[i][j] = 0;
      for (unsigned short k = 0; k < nDim; k++) { // summmation index
	for (unsigned short m = 0; m < nDim; m++) { // summation index
	  tauSGET_mom[i][j] += -xi*mu[m][k]*nvec[m]*nvec[k]*delta[i][j];
	}
	tauSGET_mom[i][j] += -xi*mu[i][k]*nvec[j]*nvec[k];
	tauSGET_mom[i][j] -= -xi*mu[i][k]*nvec[j]*nvec[k]/3.0;
	tauSGET_mom[i][j] -= -xi*mu[k][i]*nvec[j]*nvec[k]/3.0;
        tauSGET_mom[i][j] -= -xi*TWO3*mu[j][k]*nvec[i]*nvec[k];
        tauSGET_mom[i][j] += -xi*2*mu[k][k]*nvec[i]*nvec[j]/9.0;
      }
      tau_jacobian_i[i][j+1] += tauSGET_mom[i][j];
    }

    for (unsigned short j = 0; j < nDim; j++) {
      tau_jacobian_i[i][0] -= tauSGET_mom[i][j]*val_Mean_PrimVar[j+1];
    }
  
    // Jacobian w.r.t. energy
    // nothing to add
  }
}


void CAvgGrad_Hybrid::SetLaminarHeatFlux(su2double **val_gradprimvar,
                                         const su2double val_laminar_viscosity) {

  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  const su2double heat_flux_factor = Cp * (val_laminar_viscosity/Prandtl_Lam);

  /*--- Gradient of primitive variables -> [Temp vel_x vel_y vel_z Pressure] ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] = heat_flux_factor*val_gradprimvar[0][iDim];
  }
}

void CAvgGrad_Hybrid::SetLaminar_TKE_Diffusion(const su2double* const* val_gradturbvar,
                                               const su2double val_beta,
                                               const su2double val_laminar_viscosity) {
  const su2double laminar_diffusivity = val_beta * val_laminar_viscosity;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    TKE_diffusion[iDim] = laminar_diffusivity * val_gradturbvar[0][iDim];
  }
}

void CAvgGrad_Hybrid::AddSGSHeatFlux(su2double **val_gradprimvar,
                                   const su2double val_beta,
                                   const su2double val_eddy_viscosity) {

  assert(val_beta >= 0.0);
  assert(val_beta <= 1.0);

  const su2double val_alpha = pow(val_beta, 1.7);
  const su2double alpha_fac = val_alpha*(2.0 - val_alpha);

  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  //const su2double heat_flux_factor = Cp * (val_alpha*val_eddy_viscosity/Prandtl_Turb);
  const su2double heat_flux_factor = Cp * (alpha_fac*val_eddy_viscosity/Prandtl_Turb);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] += heat_flux_factor*val_gradprimvar[0][iDim];
  }

}

void CAvgGrad_Hybrid::AddSGS_TKE_Diffusion(const su2double * const* val_gradturbvar,
                                           const su2double val_beta,
                                           const su2double val_eddy_viscosity) {
 /*--- XXX: Figure out what to do with this term in hybrid.
  *
  * The limits are clear.  In RANS, it should be (mu + mu_t/sigma_k) dk/dx
  * In DNS, it should go to zero.  But the intermediate scaling is
  * not clear.  How do the modeled terms (molecular diffusion and
  * turbulent transport) scale with beta?
  *
  * It is also debateable how important these choices are.  Wilcox argues
  * that these terms are only important in hypersonic flows.
  *
  * Several possibilities:
  * 1. Scale with beta.  Simplest and easiest.
  * 2. The gradient should be w.r.t. (beta k), or k_{modeled}. This
  *    may be better, but does mean that sharp gradients in beta will
  *    produce high diffusion, even if beta ~= 0.  Is this the
  *    desired behavior? Does it present numerical difficulties?
  * 3. Scale mu_T with beta, as done for the stress. This is not a
  *    complete solution in and of itself, since the (mu * dk/dx) term
  *    should also vanish in the limit of beta -> 0. But it could
  *    be combined with the other two approaches.
  *
  * Choice 1 is taken for now.
  * ---*/
 const su2double sigma_k = 1.0;
 const su2double diffusivity = val_beta * val_eddy_viscosity/sigma_k;
 for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    TKE_diffusion[iDim] += diffusivity * val_gradturbvar[0][iDim];
  }
}

void CAvgGrad_Hybrid::AddSGETHeatFlux(su2double** val_gradprimvar,
                                      su2double** val_eddy_viscosity) {

  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  // TODO: Add this to the config file.
  const su2double Prandtl_sgs = 0.90;

  /*--- Anisotropic thermal conductivity model ---*/

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      conductivity[iDim][jDim] = Cp * val_eddy_viscosity[iDim][jDim]/Prandtl_sgs;
    }
  }

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      heat_flux_vector[iDim] += conductivity[iDim][jDim]*val_gradprimvar[0][jDim];
    }
  }
}

void CAvgGrad_Hybrid::SetHeatFluxJacobian(const su2double *val_Mean_PrimVar,
                                        const su2double val_laminar_viscosity,
                                        const su2double val_eddy_viscosity,
                                        const su2double val_dist_ij,
                                        const su2double *val_normal) {
  su2double sqvel = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double Pressure = val_Mean_PrimVar[nDim+1];
  const su2double phi = Gamma_Minus_One/Density;

  /*--- R times partial derivatives of temp. ---*/

  const su2double R_dTdu0 = -Pressure/(Density*Density) + 0.5*sqvel*phi;
  const su2double R_dTdu1 = -phi*val_Mean_PrimVar[1];
  const su2double R_dTdu2 = -phi*val_Mean_PrimVar[2];

  const su2double heat_flux_factor = val_laminar_viscosity/Prandtl_Lam + val_eddy_viscosity/Prandtl_Turb;
  const su2double cpoR = Gamma/Gamma_Minus_One; // cp over R
  const su2double conductivity_over_Rd = cpoR*heat_flux_factor/val_dist_ij;

  heat_flux_jac_i[0] = conductivity_over_Rd * R_dTdu0;
  heat_flux_jac_i[1] = conductivity_over_Rd * R_dTdu1;
  heat_flux_jac_i[2] = conductivity_over_Rd * R_dTdu2;

  if (nDim == 2) {

    const su2double R_dTdu3 = phi;
    heat_flux_jac_i[3] = conductivity_over_Rd * R_dTdu3;

  } else {

    const su2double R_dTdu3 = -phi*val_Mean_PrimVar[3];
    const su2double R_dTdu4 = phi;
    heat_flux_jac_i[3] = conductivity_over_Rd * R_dTdu3;
    heat_flux_jac_i[4] = conductivity_over_Rd * R_dTdu4;

  }
}

void CAvgGrad_Hybrid::AddSGETHeatFluxJacobian(su2double *val_Mean_PrimVar,
					      su2double **mu,
					      const su2double val_dist_ij,
					      const su2double *val_normal) {

  const su2double Prandtl_sgs = 0.90;

  su2double sqvel = 0.0;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    sqvel += val_Mean_PrimVar[iDim+1]*val_Mean_PrimVar[iDim+1];
  }

  const su2double Density = val_Mean_PrimVar[nDim+2];
  const su2double Pressure = val_Mean_PrimVar[nDim+1];
  const su2double phi = Gamma_Minus_One/Density;

  /*--- R times partial derivatives of temp. ---*/

  const su2double R_dTdu0 = -Pressure/(Density*Density) + 0.5*sqvel*phi;
  const su2double R_dTdu1 = -phi*val_Mean_PrimVar[1];
  const su2double R_dTdu2 = -phi*val_Mean_PrimVar[2];
  const su2double R_dTdu3 = -phi*val_Mean_PrimVar[3];
  const su2double R_dTdu4 = phi;

  const su2double cpoR = Gamma/Gamma_Minus_One; // cp over R
  const su2double cpoRdPr = cpoR/val_dist_ij/Prandtl_sgs;

  for (unsigned short i = 0; i < nDim; i++) { 
    for (unsigned short k = 0; k < nDim; k++) {
      heat_flux_jac_i[0] += val_normal[i]*val_normal[k]*mu[i][k]*cpoRdPr*R_dTdu0;
      heat_flux_jac_i[1] += val_normal[i]*val_normal[k]*mu[i][k]*cpoRdPr*R_dTdu1;
      heat_flux_jac_i[2] += val_normal[i]*val_normal[k]*mu[i][k]*cpoRdPr*R_dTdu2;
      heat_flux_jac_i[3] += val_normal[i]*val_normal[k]*mu[i][k]*cpoRdPr*R_dTdu3;
      heat_flux_jac_i[4] += val_normal[i]*val_normal[k]*mu[i][k]*cpoRdPr*R_dTdu4;
    }
  }

}

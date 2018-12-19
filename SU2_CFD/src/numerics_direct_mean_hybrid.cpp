/*!
 * \file numerics_direct_mean.cpp
 * \brief This file contains the numerical methods for compressible flow.
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
  assert(alpha_i == alpha_i);  // alpha_i is not NaN
  assert(alpha_j == alpha_j);  // alpha_j is not NaN

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

  /*--- Limit alpha to protect from imbalance in k_model vs k_resolved ---*/

  const su2double Mean_Alpha = min(max(0.5*(alpha_i + alpha_j), 0.0), 1.0);

  /*--- Mean gradient approximation ---*/

  for (iVar = 0; iVar < nDim+1; iVar++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_GradPrimVar[iVar][iDim] = 0.5*(PrimVar_Grad_i[iVar][iDim] + PrimVar_Grad_j[iVar][iDim]);
      Mean_GradPrimVar_Average[iVar][iDim] = 0.5*(PrimVar_Grad_Average_i[iVar][iDim] + PrimVar_Grad_Average_j[iVar][iDim]);
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
        assert(std::abs(Mean_GradPrimVar_Fluct[iVar][iDim]) < EPS);
      }
    }
  }
#endif

  /*--- Get projected flux tensor ---*/

  SetLaminarStressTensor(Mean_GradPrimVar, Mean_Laminar_Viscosity);
  AddTauSGS(Mean_PrimVar_Average, Mean_GradPrimVar_Average, Mean_Alpha,
            Mean_turb_ke, Mean_Eddy_Viscosity);
  AddTauSGET(Mean_GradPrimVar_Fluct,
             Mean_Aniso_Eddy_Viscosity);

  SetLaminarHeatFlux(Mean_GradPrimVar, Mean_Laminar_Viscosity);
  AddSGSHeatFlux(Mean_GradPrimVar_Average, Mean_Alpha, Mean_Eddy_Viscosity);
  AddSGETHeatFlux(Mean_GradPrimVar_Fluct, Mean_Aniso_Eddy_Viscosity);


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
      SetTauJacobian(Mean_PrimVar, Mean_Laminar_Viscosity, Mean_Eddy_Viscosity,
                     dist_ij, UnitNormal);
      SetHeatFluxJacobian(Mean_PrimVar, Mean_Laminar_Viscosity,
                          Mean_Eddy_Viscosity, dist_ij, UnitNormal);
      GetViscousProjJacs(Mean_PrimVar, Area, Proj_Flux_Tensor,
                         val_Jacobian_i, val_Jacobian_j);
    }

  }

}

void CAvgGrad_Hybrid::SetLaminarStressTensor(su2double **val_gradprimvar,
                                           const su2double val_laminar_viscosity) {

  unsigned short iDim, jDim;

  su2double div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] = val_laminar_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
                        - TWO3*val_laminar_viscosity*div_vel*delta[iDim][jDim];
    }
  }
}

void CAvgGrad_Hybrid::AddTauSGS(const su2double *val_primvar,
                               su2double **val_gradprimvar,
                               const su2double val_alpha,
                               const su2double val_turb_ke,
                               const su2double val_eddy_viscosity) {

  assert(val_alpha >= 0.0);
  assert(val_alpha <= 1.0);
  unsigned short iDim, jDim;
  const su2double Density = val_primvar[nDim+2];

  su2double div_vel = 0.0;
  for (iDim = 0 ; iDim < nDim; iDim++)
    div_vel += val_gradprimvar[iDim+1][iDim];

  for (iDim = 0 ; iDim < nDim; iDim++) {
    for (jDim = 0 ; jDim < nDim; jDim++) {
      tau[iDim][jDim] += val_alpha*val_eddy_viscosity*( val_gradprimvar[jDim+1][iDim] + val_gradprimvar[iDim+1][jDim] )
                        - TWO3*val_alpha*val_eddy_viscosity*div_vel*delta[iDim][jDim]
                        - TWO3*Density*val_alpha*val_turb_ke*delta[iDim][jDim];
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
        for (unsigned short lDim = 0; lDim < nDim; lDim++) {
          tau[iDim][jDim] -= TWO3*val_eddy_viscosity[kDim][lDim]*deviatoric[kDim][lDim]*delta[iDim][jDim];
        }
      }
    }
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

void CAvgGrad_Hybrid::AddSGSHeatFlux(su2double **val_gradprimvar,
                                   const su2double val_alpha,
                                   const su2double val_eddy_viscosity) {

  assert(val_alpha >= 0.0);
  assert(val_alpha <= 1.0);

  const su2double Cp = (Gamma / Gamma_Minus_One) * Gas_Constant;
  const su2double heat_flux_factor = Cp * (val_alpha*val_eddy_viscosity/Prandtl_Turb);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    heat_flux_vector[iDim] += heat_flux_factor*val_gradprimvar[0][iDim];
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

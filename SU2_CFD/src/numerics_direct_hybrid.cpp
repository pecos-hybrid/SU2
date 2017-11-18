/*!
 * \file numerics_direct_hybrid.cpp
 * \brief This file contains all the discretizations for hybrid parameter(s).
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
 * Copyright (C) 2012-2017 SU2, the open-source CFD code.
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

#include "../include/numerics_structure.hpp"

CUpwSca_HybridConv::CUpwSca_HybridConv(unsigned short val_nDim,
                                       unsigned short val_nVar,
                                       CConfig *config)
    : CNumerics(val_nDim, val_nVar, config) {

  implicit        = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible  = (config->GetKind_Regime() == INCOMPRESSIBLE);
  grid_movement   = config->GetGrid_Movement();

  Velocity_i = new su2double [nDim];
  Velocity_j = new su2double [nDim];

}

CUpwSca_HybridConv::~CUpwSca_HybridConv(void) {

  delete [] Velocity_i;
  delete [] Velocity_j;

}

void CUpwSca_HybridConv::ComputeResidual(su2double *val_residual, su2double **val_Jacobian_i, su2double **val_Jacobian_j, CConfig *config) {

  q_ij = 0.0;

  AD::StartPreacc();
  AD::SetPreaccIn(V_i, nDim+1); AD::SetPreaccIn(V_j, nDim+1);
  AD::SetPreaccIn(HybridParameter_i[0]); AD::SetPreaccIn(HybridParameter_j[0]);
  AD::SetPreaccIn(Normal, nDim);
  if (grid_movement) {
    AD::SetPreaccIn(GridVel_i, nDim); AD::SetPreaccIn(GridVel_j, nDim);
  }

  if (grid_movement) {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1] - GridVel_i[iDim];
      Velocity_j[iDim] = V_j[iDim+1] - GridVel_j[iDim];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  } else {
    for (iDim = 0; iDim < nDim; iDim++) {
      Velocity_i[iDim] = V_i[iDim+1];
      Velocity_j[iDim] = V_j[iDim+1];
      q_ij += 0.5*(Velocity_i[iDim]+Velocity_j[iDim])*Normal[iDim];
    }
  }

  a0 = 0.5*(q_ij+fabs(q_ij));
  a1 = 0.5*(q_ij-fabs(q_ij));
  val_residual[0] = a0*HybridParameter_i[0]+a1*HybridParameter_j[0];

  if (implicit) {
    val_Jacobian_i[0][0] = a0;
    val_Jacobian_j[0][0] = a1;
  }

  AD::SetPreaccOut(val_residual[0]);
  AD::EndPreacc();
}

CSourcePieceWise_HybridConv::CSourcePieceWise_HybridConv(unsigned short val_nDim, unsigned short val_nVar,
                                                 CConfig *config) : CNumerics(val_nDim, val_nVar, config) {

}

CSourcePieceWise_HybridConv::~CSourcePieceWise_HybridConv(void) { }

void CSourcePieceWise_HybridConv::ComputeResidual(su2double *val_residual,
                                                  su2double **val_Jacobian_i,
                                                  su2double **val_Jacobian_j,
                                                  CConfig *config) {

#ifndef NDEBUG
  if (Resolution_Adequacy < 0) {
    cout << "ERROR: Resolution adequacy < 0!" << endl;
    cout << "       Resolution adequacy: " << Resolution_Adequacy << endl;
    exit(EXIT_FAILURE);
  }
#endif

  /*--- Raising and lowering alpha based on resolution adequacy ---*/

  /*--- Turbulent timescale, with correction factor of 4 ---*/
  su2double T_alpha = 4.0*TurbT;

  // Gentle switch between raising and lowering alpha based on resolution
  su2double S_r;
  if (Resolution_Adequacy >= 1.0)
    S_r = tanh(Resolution_Adequacy - 1.0);
  else
    S_r = tanh(1.0 - 1.0/Resolution_Adequacy);

  /*--- D_r causes alpha > 1 to return rapidly to 1 in RANS regions.
   * But if RANS_Weight << 1, alpha can go above 1. ---*/
  su2double D_r = RANS_Weight * fmax(HybridParameter_i[0] + S_r - 1, 0);

  /*--- F_r causes alpha to return to 0 if it goes briefly negative.
   * F_r functions as a built-in clipping which maintains C1 continuity  ---*/
  su2double F_r = fmin(HybridParameter_i[0] + S_r, 0);

  val_residual[0] = (S_r - D_r - F_r)/T_alpha * Volume;

  /*--- Raising and lowering alpha based on grid coarsening/refinement ---*/

  su2double udMdx; // Dot product of resolved velocity and resolution tensor gradient
  su2double max_udMdx; // Max value of the convection of the resolution tensor
  for (unsigned int iDim = 0; iDim < nDim; iDim++) {
    for (unsigned int jDim = 0; jDim < nDim; jDim++) {
      udMdx = 0.0;
      for (unsigned int kDim = 0; kDim < nDim; kDim++) {
        // This is NOT upwinded.
        udMdx += V_i[kDim+1] * Resolution_Tensor_Gradient[kDim][iDim][jDim];
      }
      if (iDim == 0 and jDim == 0) max_udMdx = udMdx;
      else if (udMdx > max_udMdx) max_udMdx = udMdx;
    }
  }

  if (abs(max_udMdx) > EPS) {
    // Timescale for grid refinement/coarsening
    su2double T_c = TurbL/max_udMdx;
    if (S_r >= 0.0 && T_c >= 0.0) {
      val_residual[0] += 1.0/T_c * Volume;
    }
  }

  /*--- Jacobian of \alpha with respect to \alpha ---*/
  val_Jacobian_i[0][0] = 0.0;
  if (D_r != 0)
    val_Jacobian_i[0][0] += -RANS_Weight/T_alpha * Volume;
  if (F_r != 0)
    val_Jacobian_i[0][0] += -1.0/T_alpha * Volume;
}

CAvgGrad_HybridConv::CAvgGrad_HybridConv(unsigned short val_nDim,
                                         unsigned short val_nVar,
                                         bool correct_grad,
                                         CConfig *config)
    : CNumerics(val_nDim, val_nVar, config),
      sigma_alpha(10.),
      correct_gradient(correct_grad) {

  implicit = (config->GetKind_TimeIntScheme_Turb() == EULER_IMPLICIT);
  incompressible = (config->GetKind_Regime() == INCOMPRESSIBLE);

  Edge_Vector = new su2double [nDim];
  Proj_Mean_Grad_Normal = new su2double [nVar];
  Proj_Mean_Grad_Edge = new su2double [nVar];
  Proj_Mean_Grad = new su2double [nVar];
  Mean_Grad = new su2double* [nVar];
  for (iVar = 0; iVar < nVar; iVar++)
    Mean_Grad[iVar] = new su2double [nDim];

}

CAvgGrad_HybridConv::~CAvgGrad_HybridConv(void) {

  delete [] Edge_Vector;
  delete [] Proj_Mean_Grad_Normal;
  delete [] Proj_Mean_Grad_Edge;
  delete [] Proj_Mean_Grad;
  for (iVar = 0; iVar < nVar; iVar++)
    delete [] Mean_Grad[iVar];
  delete [] Mean_Grad;

}


void CAvgGrad_HybridConv::ComputeResidual(su2double *val_residual,
                                          su2double **Jacobian_i,
                                          su2double **Jacobian_j,
                                          CConfig *config) {

  AD::StartPreacc();
  AD::SetPreaccIn(Coord_i, nDim); AD::SetPreaccIn(Coord_j, nDim);
  AD::SetPreaccIn(Normal, nDim);
  AD::SetPreaccIn(HybridParam_Grad_i, nVar, nDim);
  AD::SetPreaccIn(HybridParam_Grad_j, nVar, nDim);

  if (incompressible) {
    AD::SetPreaccIn(V_i, nDim+5); AD::SetPreaccIn(V_j, nDim+5);

    Density_i = V_i[nDim+1];            Density_j = V_j[nDim+1];
    Laminar_Viscosity_i = V_i[nDim+3];  Laminar_Viscosity_j = V_j[nDim+3];
    Eddy_Viscosity_i = V_i[nDim+4];     Eddy_Viscosity_j = V_j[nDim+4];
  }
  else {
    AD::SetPreaccIn(V_i, nDim+7); AD::SetPreaccIn(V_j, nDim+7);

    Density_i = V_i[nDim+2];            Density_j = V_j[nDim+2];
    Laminar_Viscosity_i = V_i[nDim+5];  Laminar_Viscosity_j = V_j[nDim+5];
    Eddy_Viscosity_i = V_i[nDim+6];     Eddy_Viscosity_j = V_j[nDim+6];
  }

  /*--- Compute vector going from iPoint to jPoint ---*/

  dist_ij_2 = 0; proj_vector_ij = 0;
  for (iDim = 0; iDim < nDim; iDim++) {
    Edge_Vector[iDim] = Coord_j[iDim]-Coord_i[iDim];
    dist_ij_2 += Edge_Vector[iDim]*Edge_Vector[iDim];
    proj_vector_ij += Edge_Vector[iDim]*Normal[iDim];
  }
  if (dist_ij_2 == 0.0) proj_vector_ij = 0.0;
  else proj_vector_ij = proj_vector_ij/dist_ij_2;

  /*--- Mean gradient approximation ---*/
  for (iVar = 0; iVar < nVar; iVar++) {
    Proj_Mean_Grad_Normal[iVar] = 0.0;
    Proj_Mean_Grad_Edge[iVar] = 0.0;
    for (iDim = 0; iDim < nDim; iDim++) {
      Mean_Grad[iVar][iDim] = 0.5*(HybridParam_Grad_i[iVar][iDim] +
                                   HybridParam_Grad_j[iVar][iDim]);
      Proj_Mean_Grad_Normal[iVar] += Mean_Grad[iVar][iDim] *
                                     Normal[iDim];
      if (correct_gradient)
        Proj_Mean_Grad_Edge[iVar] += Mean_Grad[iVar][iDim]*Edge_Vector[iDim];
    }
    Proj_Mean_Grad[iVar] = Proj_Mean_Grad_Normal[iVar];
    if (correct_gradient) {
      Proj_Mean_Grad[iVar] -= Proj_Mean_Grad_Edge[iVar]*proj_vector_ij -
      (HybridParameter_j[iVar]-HybridParameter_i[iVar])*proj_vector_ij;
    }
  }

  /*--- Compute mean effective viscosity ---*/
  nu_i = Laminar_Viscosity_i + Eddy_Viscosity_i/sigma_alpha;
  nu_j = Laminar_Viscosity_j + Eddy_Viscosity_j/sigma_alpha;

  val_residual[0] = 0.5*(nu_i + nu_j)*Proj_Mean_Grad[0];

  /*--- For Jacobians -> Use of TSL approx. to compute derivatives of the gradients ---*/

  if (implicit) {
    Jacobian_i[0][0] = -0.5*(nu_i + nu_j)*proj_vector_ij;
    Jacobian_j[0][0] =  0.5*(nu_i + nu_j)*proj_vector_ij;
  }

  AD::SetPreaccOut(val_residual, nVar);
  AD::EndPreacc();


}


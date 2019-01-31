/*!
 * \file fluctuating_stress.cpp
 * \brief Implementation of the fluctuating stress for RANS/LES hybrid
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

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../include/fluctuating_stress.hpp"

CFluctuatingStress::CFluctuatingStress(unsigned short val_nDim)
    : nDim(val_nDim) {
}

CFluctuatingStress::~CFluctuatingStress() {
}

void CFluctuatingStress::SetTurbVar(const su2double* val_turb_var) {
  TurbVar = val_turb_var;
}

void CFluctuatingStress::SetPrimitive(const su2double* val_primitive) {
  FlowPrimVar = val_primitive;
}

CM43Model::CM43Model(unsigned short val_nDim,
                     CConfig* config)
    : CFluctuatingStress(val_nDim) {

  if (config->GetKind_Regime() != COMPRESSIBLE)
    SU2_MPI::Error("The M43 model has not been set up for incompressible flow!", CURRENT_FUNCTION);

  delta = new su2double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    delta[iDim] = new su2double [nDim];
  }

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      if (iDim == jDim) delta[iDim][jDim] = 1.0;
      else delta[iDim][jDim]=0.0;
    }
  }
}

CM43Model::~CM43Model() {
  if (delta != NULL) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      delete [] delta[iDim];
    }
    delete [] delta;
  }
}

void CM43Model::CalculateEddyViscosity(const CGeometry* geometry,
                                       CConfig* config,
                                       const unsigned long iPoint,
                                       const su2double mean_eddy_visc,
                                       su2double** eddy_viscosity) const {

  /*-- Retrieve necessary variables ---*/

  const su2double C_M = geometry->node[iPoint]->GetResolutionCoeff();
  const su2double* const* M43 = geometry->node[iPoint]->GetResolutionTensor43();

  const su2double density = FlowPrimVar[nDim+2];
  su2double dissipation = 0;
  switch (config->GetKind_Turb_Model()) {
    case KE:
     dissipation = max(0.0, TurbVar[1]); break;
    default:
      SU2_MPI::Error("The M43 model has not been set up for your RANS model.", CURRENT_FUNCTION);
  }

  /*--- Check the preconditions ---*/

  assert(density > 0);
  assert(M43 != NULL);
  assert(M43[0] != NULL);
  assert(eddy_viscosity != NULL);
  assert(C_M > 0);

  /*--- C_M0 is an overall coefficient used to calibrate the model to match
   * isotropic resolution ---*/
  const su2double C_M0 = 0.13;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      eddy_viscosity[iDim][jDim] =
          density * C_M0*C_M * pow(dissipation, 1.0/3) * M43[iDim][jDim];
    }
  }

  /*--- Blend with the RANS eddy viscosity for high AR cells ---*/

  /*--- Calculate the maximum aspect ratio ---*/
  // TODO: Move this aspect ratio calc to the dual grid class.
//  const su2double* resolution_values = geometry->node[iPoint]->GetResolutionValues();
//  su2double min_distance = resolution_values[0];
//  su2double max_distance = resolution_values[0];
//  for (unsigned short iDim = 1; iDim < nDim; iDim++) {
//    min_distance = min(min_distance, resolution_values[iDim]);
//    max_distance = max(max_distance, resolution_values[iDim]);
//  }
//  const su2double aspect_ratio = max_distance / min_distance;
//  assert(aspect_ratio >= 1.00);
//
//  const su2double AR_switch = 50;
//  if (aspect_ratio > AR_switch) {
//    const su2double blending = tanh((aspect_ratio - AR_switch)/10.0);
//    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
//      for (unsigned short jDim = 0; jDim < nDim; jDim++) {
//        eddy_viscosity[iDim][jDim] *= (1 - blending);
//        eddy_viscosity[iDim][jDim] += blending*delta[iDim][jDim]*mean_eddy_visc;
//      }
//    }
//  }
}

CNoStressModel::CNoStressModel(unsigned short val_nDim)
  : CFluctuatingStress(val_nDim) {
}

void CNoStressModel::CalculateEddyViscosity(const CGeometry* geometry,
                                            CConfig* config,
                                            const unsigned long iPoint,
                                            const su2double mean_eddy_visc,
                                            su2double** eddy_viscosity) const {

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      eddy_viscosity[iDim][jDim] = 0;
    }
  }
}

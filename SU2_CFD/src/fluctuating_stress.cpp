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
}

void CM43Model::CalculateEddyViscosity(const CGeometry* geometry,
                                       CConfig* config,
                                       unsigned short iPoint,
                                       su2double** eddy_viscosity) const {

  /*-- Retrieve necessary variables ---*/

  const su2double C_M = geometry->node[iPoint]->GetResolutionCoeff();
  su2double** M43 = geometry->node[iPoint]->GetResolutionTensor43();

  const su2double density = FlowPrimVar[nDim+2];
  su2double dissipation = 0;
  switch (config->GetKind_Turb_Model()) {
    case KE:
     dissipation = TurbVar[1]; break;
    default:
      SU2_MPI::Error("The M43 model has not been set up for your RANS model.", CURRENT_FUNCTION);
  }

  /*--- Check the preconditions ---*/

  assert(density > 0);
  assert(dissipation >= 0);
  assert(M43 != NULL);
  assert(M43[0] != NULL);
  assert(C_M > 0);

  /*--- C_M0 is an overall coefficient used to calibrate the model to match
   * isotropic resolution ---*/
  const su2double C_M0 = 0.12;

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      eddy_viscosity[iDim][jDim] =
          density * C_M0*C_M * pow(dissipation, 1.0/3) * M43[iDim][jDim];
    }
  }
}
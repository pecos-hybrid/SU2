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

CFluctuatingStress::CFluctuatingStress(CSolver **solver_container,
                                       CGeometry* geometry)
    : nDim(geometry->GetnDim()) {

  /*--- Allocate space for primitive variables ---*/

  unsigned short nPrimVar = solver_container[FLOW_SOL]->GetnPrimVar();
  FlowPrimVar = new su2double[nPrimVar];

  /*--- Allocate space for the turbulence variables ---*/

  unsigned short nTurbVar = solver_container[TURB_SOL]->GetnVar();
  MeanTurbVar = new su2double[nPrimVar];
}

CFluctuatingStress::~CFluctuatingStress() {
  if (FlowPrimVar != NULL) delete [] FlowPrimVar;
  if (MeanTurbVar != NULL) delete [] MeanTurbVar;
}

void CFluctuatingStress::SetMeanTurbVars(const su2double* turb_vars) {
  MeanTurbVar = turb_vars;
}

void CFluctuatingStress::SetPrimitive(const su2double* val_primitive) {
  FlowPrimVar = val_primitive;
}

CM43Model::CM43Model(CSolver **solver_container, CGeometry* geometry,
                     CConfig* config)
    : CFluctuatingStress(solver_container, geometry) {

  if (config->GetKind_Regime() == COMPRESSIBLE) {
    density_index = 0;
  } else {
    SU2_MPI::Error("The M43 model has not been set up for incompressible flow!", CURRENT_FUNCTION);
  }
}

void CM43Model::CalculateEddyViscosity(const CGeometry* geometry,
                                       CConfig* config,
                                       unsigned short iPoint,
                                       su2double** eddy_viscosity) const {

  /*-- Retrieve necessary variables ---*/

  const su2double C_M = geometry->node[iPoint]->GetResolutionCoeff();
  const su2double** M43 = geometry->node[iPoint]->GetResolutionTensor43();

  const su2double density = FlowPrimVar[density_index];
  su2double dissipation = 0;
  switch (config->GetKind_Turb_Model()) {
    case KE:
     dissipation = MeanTurbVar[1]; break;
    default:
      SU2_MPI::Error("The M43 model has not been set up for your RANS model.", CURRENT_FUNCTION);
  }

  /*--- Check the preconditions ---*/

  assert(density > 0);
  assert(dissipation >= 0);
  assert(M43 != NULL);
  assert(M43[0] != NULL);
  assert(C_M > 0);

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      eddy_viscosity[iDim][jDim] =
          density * C_M * pow(dissipation, 1.0/3) * M43[iDim][jDim];
    }
  }
}

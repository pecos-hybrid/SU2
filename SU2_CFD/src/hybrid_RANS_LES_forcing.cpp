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

CHybridForcingAbstractBase::CHybridForcingAbstractBase(
                              const unsigned short nDim,
                              const unsigned long nPoint,
                              const unsigned long nPointDomain)
  : nDim(nDim), nVar(nDim), nVarGrad(nDim),
    nPoint(nPoint), nPointDomain(nPointDomain) {

  LeviCivita = new int**[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    LeviCivita[iDim] = new int*[nDim];
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      LeviCivita[iDim][jDim] = new int[nDim];
      for (unsigned short kDim = 0; kDim < nDim; kDim++) {
        LeviCivita[iDim][jDim][kDim] = 0;
      }
    }
  }
  /*--- These could be assigned programmatically, but the straightforward
   * way is less error-prone. ---*/
  LeviCivita[0][1][2] = 1;
  LeviCivita[1][2][0] = 1;
  LeviCivita[2][0][1] = 1;
  LeviCivita[0][2][1] = -1;
  LeviCivita[2][1][0] = -1;
  LeviCivita[1][0][2] = -1;
}

CHybridForcingAbstractBase::~CHybridForcingAbstractBase() {

  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      delete [] LeviCivita[iDim][jDim];
    }
    delete [] LeviCivita[iDim];
  }
  delete [] LeviCivita;

}

CHybridForcingTGSF::CHybridForcingTGSF(const unsigned short nDim,
                               const unsigned long nPoint,
                               const unsigned long nPointDomain)
  : CHybridForcingAbstractBase(nDim, nPoint, nPointDomain) {

  node = NULL;
  Gradient = NULL;

  Smatrix = new su2double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];

  Cvector = new su2double* [nVarGrad];
  for (unsigned short iVar = 0; iVar < nVarGrad; iVar++)
    Cvector[iVar] = new su2double [nDim];
}


CHybridForcingTGSF::CHybridForcingTGSF(CGeometry* geometry, CConfig* config)
  : CHybridForcingAbstractBase(geometry->GetnDim(),
                               geometry->GetnPoint(),
                               geometry->GetnPointDomain()) {

  node = new su2double*[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new su2double[nVar];
  }

  Gradient = new su2double**[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
     Gradient[iPoint] = new su2double*[nVar];
     for (unsigned short iVar = 0; iVar < nVar; iVar++) {
       Gradient[iPoint][iVar] = new su2double[nDim];
     }
  }

  Smatrix = new su2double* [nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    Smatrix[iDim] = new su2double [nDim];

  Cvector = new su2double* [nVarGrad];
  for (unsigned short iVar = 0; iVar < nVarGrad; iVar++)
    Cvector[iVar] = new su2double [nDim];

}

CHybridForcingTGSF::~CHybridForcingTGSF() {

  if (node != NULL) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      delete [] node[iPoint];
    }
    delete [] node;
  }

  if (Gradient != NULL) {
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
       for (unsigned short iVar = 0; iVar < nVar; iVar++) {
         delete [] Gradient[iPoint][iVar];
       }
       delete [] Gradient[iPoint];
    }
    delete [] Gradient;
  }

  if (Smatrix != NULL) {
    for (unsigned short iDim = 0; iDim < nDim; iDim++)
      delete [] Smatrix[iDim];
    delete [] Smatrix;
  }

  if (Cvector != NULL) {
    for (unsigned short iVar = 0; iVar < nVarGrad; iVar++)
      delete [] Cvector[iVar];
    delete [] Cvector;
  }

}

su2double CHybridForcingTGSF::ComputeScalingFactor(const su2double L,
                                               const su2double P_F,
                                               const su2double dt,
                                               const su2double* b) {
  const su2double A = 3;

  su2double k_b = 0.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    k_b += b[iDim]*b[iDim];
  }
  k_b /= A*A;

  const su2double C_alpha = -1*min(P_F, 0.0)*k_b; // FIXME: If statements?

  const su2double pi = atan(1.0)*4.0;
  const su2double ax = 2*pi/L; // TODO: Check differences in ax, ay, az.
  const su2double ay = 2*pi/L;
  const su2double az = 2*pi/L;

  const su2double C_P = 1.0; // FIXME: Not defined in model doc.
  return C_P * sqrt(2*abs(C_alpha)*dt)/max(max(ax, ay), az);
}

void CHybridForcingTGSF::ComputeForcingField(CSolver** solver, CGeometry *geometry,
                                         CConfig *config) {

  const unsigned short kind_time_marching = config->GetUnsteady_Simulation();

  if (kind_time_marching == TIME_STEPPING ||
      kind_time_marching == DT_STEPPING_1ST ||
      kind_time_marching == DT_STEPPING_2ND) {

    // FIXME: Current UnstTime or Total UnstTime???
    const su2double time = config->GetTotal_UnstTimeND();
    const su2double dt = config->GetDelta_UnstTimeND();

    /*--- Allocate some scratch arrays to avoid continual reallocation ---*/

    su2double b[nDim]; // Initial TG vortex field.
    su2double h[nDim]; // Stream function from initial TG vortex field
    su2double L[nDim]; // Length scales
    su2double x[nDim]; // Position

    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

      const su2double timescale =
          solver[TURB_SOL]->node[iPoint]->GetTurbTimescale();
      /*--- This velocity retrieval works only for compressible ---*/
      assert(config->GetKind_Regime() == COMPRESSIBLE);
      const su2double* prim_vars =
          solver[FLOW_SOL]->node[iPoint]->GetPrimitive();
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        /*--- Copy the values (and not the pointer), since we're changing them. ---*/
        x[iDim] = geometry->node[iPoint]->GetCoord(iDim);
        x[iDim] = TransformCoords(x[iDim], prim_vars[iDim+1], time, timescale);
      }

      /*--- Setup the TG vortex and its stream function ---*/

      // TODO: The length scales are defined in three different places.
      // Eliminate the repitition.
      const su2double N_L = 4.0;
      const su2double L_iso =
          N_L*solver[TURB_SOL]->node[iPoint]->GetTurbLengthscale();
      L[0] = L_iso; L[1] = L_iso; L[2] = L_iso;
      SetTGField(x, L, b);
      SetStreamFunc(x, L, h);

      /*--- Compute the scaling factor for the TG vortex field ---*/

      const su2double k_sgs = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
      // FIXME: Where is k_total stored?
      const su2double k_total = k_sgs;
      const su2double dissipation = solver[TURB_SOL]->node[iPoint]->GetSolution(1);
      // FIXME: Where is average r_M stored?
      const su2double resolution_adequacy = 1.0;
      const su2double alpha = k_sgs / k_total;
      const su2double density = solver[FLOW_SOL]->node[iPoint]->GetSolution(0);
      const su2double laminar_viscosity =
          solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity() / density;
      const su2double P_F = GetTargetProduction(k_sgs, dissipation,
                                                resolution_adequacy, alpha,
                                                laminar_viscosity);
      const su2double eta = ComputeScalingFactor(L_iso, P_F, dt, b);

      /*--- Store eta*h so we can compute the derivatives ---*/

      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        node[iPoint][iVar] = eta*h[iVar];
      }
    }

    SetForcing_Gradient_LS(geometry, config);

    // TODO: Decide if this should be passed in or stored in forcing class
    su2double** g;
    for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
      for (unsigned short iDim = 0; iDim < nDim; iDim++) {
        g[iPoint][iDim] = 0.0;
        for (unsigned short jDim = 0; jDim < nDim; jDim++) {
          for (unsigned short kDim = 0; kDim < nDim; kDim++) {
            g[iPoint][iDim] +=
              LeviCivita[iDim][jDim][kDim]*Gradient[iPoint][kDim][jDim];
          }
        }
      }
    }
  } else {
    SU2_MPI::Error("Hybrid forcing has not been set up for steady flows.",
                   CURRENT_FUNCTION);
  }

}

void CHybridForcingTGSF::SetForcing_Gradient_LS(CGeometry *geometry, CConfig *config) {

  unsigned short iVar, iDim, jDim, iNeigh;
  unsigned long iPoint, jPoint;
  su2double *Coord_i, *Coord_j, r11, r12, r13, r22, r23, r23_a,
  r23_b, r33, weight, product, z11, z12, z13, z22, z23, z33, detR2;
  bool singular;

  assert(node != NULL);
  assert(Gradient != NULL);

  /*--- Loop over points of the grid ---*/

  for (iPoint = 0; iPoint < nPointDomain; iPoint++) {

    /*--- Set the value of the singular ---*/
    singular = false;

    /*--- Get coordinates ---*/

    Coord_i = geometry->node[iPoint]->GetCoord();

    /*--- Get variables from array ---*/

    su2double* StreamFunc_i = node[iPoint];

    /*--- Inizialization of variables ---*/

    for (iVar = 0; iVar < nVarGrad; iVar++)
      for (iDim = 0; iDim < nDim; iDim++)
        Cvector[iVar][iDim] = 0.0;

    r11 = 0.0; r12 = 0.0;   r13 = 0.0;    r22 = 0.0;
    r23 = 0.0; r23_a = 0.0; r23_b = 0.0;  r33 = 0.0;

//    AD::StartPreacc();
//    AD::SetPreaccIn(PrimVar_i, nPrimVarGrad);
//    AD::SetPreaccIn(Coord_i, nDim);

    for (iNeigh = 0; iNeigh < geometry->node[iPoint]->GetnPoint(); iNeigh++) {
      jPoint = geometry->node[iPoint]->GetPoint(iNeigh);
      Coord_j = geometry->node[jPoint]->GetCoord();

      su2double* StreamFunc_j = node[jPoint];

//      AD::SetPreaccIn(Coord_j, nDim);
//      AD::SetPreaccIn(PrimVar_j, nPrimVarGrad);

      weight = 0.0;
      for (iDim = 0; iDim < nDim; iDim++)
        weight += (Coord_j[iDim]-Coord_i[iDim])*(Coord_j[iDim]-Coord_i[iDim]);

      /*--- Sumations for entries of upper triangular matrix R ---*/

      if (weight != 0.0) {

        r11 += (Coord_j[0]-Coord_i[0])*(Coord_j[0]-Coord_i[0])/weight;
        r12 += (Coord_j[0]-Coord_i[0])*(Coord_j[1]-Coord_i[1])/weight;
        r22 += (Coord_j[1]-Coord_i[1])*(Coord_j[1]-Coord_i[1])/weight;

        if (nDim == 3) {
          r13 += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r23_a += (Coord_j[1]-Coord_i[1])*(Coord_j[2]-Coord_i[2])/weight;
          r23_b += (Coord_j[0]-Coord_i[0])*(Coord_j[2]-Coord_i[2])/weight;
          r33 += (Coord_j[2]-Coord_i[2])*(Coord_j[2]-Coord_i[2])/weight;
        }

        /*--- Entries of c:= transpose(A)*b ---*/

        for (iVar = 0; iVar < nVarGrad; iVar++)
          for (iDim = 0; iDim < nDim; iDim++)
            Cvector[iVar][iDim] += (Coord_j[iDim]-Coord_i[iDim])*(StreamFunc_j[iVar]-StreamFunc_i[iVar])/weight;

      }

    }

    /*--- Entries of upper triangular matrix R ---*/

    if (r11 >= 0.0) r11 = sqrt(r11); else r11 = 0.0;
    if (r11 != 0.0) r12 = r12/r11; else r12 = 0.0;
    if (r22-r12*r12 >= 0.0) r22 = sqrt(r22-r12*r12); else r22 = 0.0;

    if (nDim == 3) {
      if (r11 != 0.0) r13 = r13/r11; else r13 = 0.0;
      if ((r22 != 0.0) && (r11*r22 != 0.0)) r23 = r23_a/r22 - r23_b*r12/(r11*r22); else r23 = 0.0;
      if (r33-r23*r23-r13*r13 >= 0.0) r33 = sqrt(r33-r23*r23-r13*r13); else r33 = 0.0;
    }

    /*--- Compute determinant ---*/

    if (nDim == 2) detR2 = (r11*r22)*(r11*r22);
    else detR2 = (r11*r22*r33)*(r11*r22*r33);

    /*--- Detect singular matrices ---*/

    if (abs(detR2) <= EPS) { detR2 = 1.0; singular = true; }

    /*--- S matrix := inv(R)*traspose(inv(R)) ---*/

    if (singular) {
      for (iDim = 0; iDim < nDim; iDim++)
        for (jDim = 0; jDim < nDim; jDim++)
          Smatrix[iDim][jDim] = 0.0;
    }
    else {
      if (nDim == 2) {
        Smatrix[0][0] = (r12*r12+r22*r22)/detR2;
        Smatrix[0][1] = -r11*r12/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = r11*r11/detR2;
      }
      else {
        z11 = r22*r33; z12 = -r12*r33; z13 = r12*r23-r13*r22;
        z22 = r11*r33; z23 = -r11*r23; z33 = r11*r22;
        Smatrix[0][0] = (z11*z11+z12*z12+z13*z13)/detR2;
        Smatrix[0][1] = (z12*z22+z13*z23)/detR2;
        Smatrix[0][2] = (z13*z33)/detR2;
        Smatrix[1][0] = Smatrix[0][1];
        Smatrix[1][1] = (z22*z22+z23*z23)/detR2;
        Smatrix[1][2] = (z23*z33)/detR2;
        Smatrix[2][0] = Smatrix[0][2];
        Smatrix[2][1] = Smatrix[1][2];
        Smatrix[2][2] = (z33*z33)/detR2;
      }
    }

    /*--- Computation of the gradient: S*c ---*/
    for (iVar = 0; iVar < nVarGrad; iVar++) {
      for (iDim = 0; iDim < nDim; iDim++) {
        product = 0.0;
        for (jDim = 0; jDim < nDim; jDim++) {
          product += Smatrix[iDim][jDim]*Cvector[iVar][jDim];
        }
        Gradient[iPoint][iVar][iDim] = product;
      }
    }

//    AD::SetPreaccOut(node[iPoint]->GetGradient_Primitive(), nPrimVarGrad, nDim);
//    AD::EndPreacc();
  }

  Set_MPI_Forcing_Gradient(geometry, config);

}

void CHybridForcingTGSF::Set_MPI_Forcing_Gradient(CGeometry *geometry, CConfig *config) {
//  unsigned short iVar, iDim, iMarker, iPeriodic_Index, MarkerS, MarkerR;
//  unsigned long iVertex, iPoint, nVertexS, nVertexR, nBufferS_Vector, nBufferR_Vector;
//  su2double rotMatrix[3][3], *angles, theta, cosTheta, sinTheta, phi, cosPhi, sinPhi, psi, cosPsi, sinPsi,
//  *Buffer_Receive_Gradient = NULL, *Buffer_Send_Gradient = NULL;
//
//  su2double **Gradient = new su2double* [nPrimVarGrad];
//  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
//    Gradient[iVar] = new su2double[nDim];
//
//#ifdef HAVE_MPI
//  int send_to, receive_from;
//  SU2_MPI::Status status;
//#endif
//
//  for (iMarker = 0; iMarker < nMarker; iMarker++) {
//
//    if ((config->GetMarker_All_KindBC(iMarker) == SEND_RECEIVE) &&
//        (config->GetMarker_All_SendRecv(iMarker) > 0)) {
//
//      MarkerS = iMarker;  MarkerR = iMarker+1;
//
//#ifdef HAVE_MPI
//      send_to = config->GetMarker_All_SendRecv(MarkerS)-1;
//      receive_from = abs(config->GetMarker_All_SendRecv(MarkerR))-1;
//#endif
//
//      nVertexS = geometry->nVertex[MarkerS];  nVertexR = geometry->nVertex[MarkerR];
//      nBufferS_Vector = nVertexS*nPrimVarGrad*nDim;        nBufferR_Vector = nVertexR*nPrimVarGrad*nDim;
//
//      /*--- Allocate Receive and send buffers  ---*/
//      Buffer_Receive_Gradient = new su2double [nBufferR_Vector];
//      Buffer_Send_Gradient = new su2double[nBufferS_Vector];
//
//      /*--- Copy the solution old that should be sended ---*/
//      for (iVertex = 0; iVertex < nVertexS; iVertex++) {
//        iPoint = geometry->vertex[MarkerS][iVertex]->GetNode();
//        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexS+iVar*nVertexS+iVertex] = node[iPoint]->GetGradient_Primitive(iVar, iDim);
//      }
//
//#ifdef HAVE_MPI
//
//      /*--- Send/Receive information using Sendrecv ---*/
//      SU2_MPI::Sendrecv(Buffer_Send_Gradient, nBufferS_Vector, MPI_DOUBLE, send_to, 0,
//                        Buffer_Receive_Gradient, nBufferR_Vector, MPI_DOUBLE, receive_from, 0, MPI_COMM_WORLD, &status);
//
//#else
//
//      /*--- Receive information without MPI ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] = Buffer_Send_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//      }
//
//#endif
//
//      /*--- Deallocate send buffer ---*/
//      delete [] Buffer_Send_Gradient;
//
//      /*--- Do the coordinate transformation ---*/
//      for (iVertex = 0; iVertex < nVertexR; iVertex++) {
//
//        /*--- Find point and its type of transformation ---*/
//        iPoint = geometry->vertex[MarkerR][iVertex]->GetNode();
//        iPeriodic_Index = geometry->vertex[MarkerR][iVertex]->GetRotation_Type();
//
//        /*--- Retrieve the supplied periodic information. ---*/
//        angles = config->GetPeriodicRotation(iPeriodic_Index);
//
//        /*--- Store angles separately for clarity. ---*/
//        theta    = angles[0];   phi    = angles[1];     psi    = angles[2];
//        cosTheta = cos(theta);  cosPhi = cos(phi);      cosPsi = cos(psi);
//        sinTheta = sin(theta);  sinPhi = sin(phi);      sinPsi = sin(psi);
//
//        /*--- Compute the rotation matrix. Note that the implicit
//         ordering is rotation about the x-axis, y-axis,
//         then z-axis. Note that this is the transpose of the matrix
//         used during the preprocessing stage. ---*/
//        rotMatrix[0][0] = cosPhi*cosPsi;    rotMatrix[1][0] = sinTheta*sinPhi*cosPsi - cosTheta*sinPsi;     rotMatrix[2][0] = cosTheta*sinPhi*cosPsi + sinTheta*sinPsi;
//        rotMatrix[0][1] = cosPhi*sinPsi;    rotMatrix[1][1] = sinTheta*sinPhi*sinPsi + cosTheta*cosPsi;     rotMatrix[2][1] = cosTheta*sinPhi*sinPsi - sinTheta*cosPsi;
//        rotMatrix[0][2] = -sinPhi;          rotMatrix[1][2] = sinTheta*cosPhi;                              rotMatrix[2][2] = cosTheta*cosPhi;
//
//        /*--- Copy conserved variables before performing transformation. ---*/
//        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            Gradient[iVar][iDim] = Buffer_Receive_Gradient[iDim*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//
//        /*--- Need to rotate the gradients for all conserved variables. ---*/
//        for (iVar = 0; iVar < nPrimVarGrad; iVar++) {
//          if (nDim == 2) {
//            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//          }
//          else {
//            Gradient[iVar][0] = rotMatrix[0][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[0][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][1] = rotMatrix[1][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[1][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//            Gradient[iVar][2] = rotMatrix[2][0]*Buffer_Receive_Gradient[0*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][1]*Buffer_Receive_Gradient[1*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex] + rotMatrix[2][2]*Buffer_Receive_Gradient[2*nPrimVarGrad*nVertexR+iVar*nVertexR+iVertex];
//          }
//        }
//
//        /*--- Store the received information ---*/
//        for (iVar = 0; iVar < nPrimVarGrad; iVar++)
//          for (iDim = 0; iDim < nDim; iDim++)
//            node[iPoint]->SetGradient_Primitive(iVar, iDim, Gradient[iVar][iDim]);
//
//      }
//
//      /*--- Deallocate receive buffer ---*/
//      delete [] Buffer_Receive_Gradient;
//
//    }
//
//  }
//
//  for (iVar = 0; iVar < nPrimVarGrad; iVar++)
//    delete [] Gradient[iVar];
//  delete [] Gradient;

}


CHybridForcingTG0::CHybridForcingTG0(const unsigned short nDim,
                               const unsigned long nPoint,
                               const unsigned long nPointDomain)
  : CHybridForcingAbstractBase(nDim, nPoint, nPointDomain) {

  node = new su2double*[nPoint];
  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {
    node[iPoint] = new su2double[nVar];
  }

}


CHybridForcingTG0::CHybridForcingTG0(CGeometry* geometry, CConfig* config)
  : CHybridForcingAbstractBase(geometry->GetnDim(),
                               geometry->GetnPoint(),
                               geometry->GetnPointDomain()) {

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

  const unsigned short kind_time_marching = config->GetUnsteady_Simulation();

  assert(kind_time_marching == TIME_STEPPING   ||
         kind_time_marching == DT_STEPPING_1ST ||
         kind_time_marching == DT_STEPPING_2ND );

  // FIXME: Current UnstTime or Total UnstTime???
  const su2double time = config->GetTotal_UnstTimeND();
  const su2double dt = config->GetDelta_UnstTimeND();

  /*--- Allocate some scratch arrays to avoid continual reallocation ---*/
  su2double h[nDim]; // Initial TG vortex field.
  su2double Lsgs; // SGS length scale
  su2double x[nDim]; // Position
  su2double Lmesh[nDim]; // Mesh length scales in coord direction (computed from res tensor)
  su2double dwall; // distance to the wall
  su2double uf[nDim];

  // Domain lengths for periodic directions
  su2double *D = config->GetHybrid_Forcing_Periodic_Length();


  for (unsigned long iPoint = 0; iPoint < nPoint; iPoint++) {

    const su2double timescale =
      solver[TURB_SOL]->node[iPoint]->GetTurbTimescale();

    /*--- This velocity retrieval works only for compressible ---*/
    assert(config->GetKind_Regime() == COMPRESSIBLE);
    const su2double* prim_vars = solver[FLOW_SOL]->node[iPoint]->GetPrimitive();
    const su2double* prim_mean =
      solver[FLOW_SOL]->average_node[iPoint]->GetPrimitive();

    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      /*--- Copy the values (and not the pointer), since we're changing them. ---*/
      x[iDim] = geometry->node[iPoint]->GetCoord(iDim);
      x[iDim] = TransformCoords(x[iDim], prim_mean[iDim+1], time, timescale);
    }

    /*--- Setup the TG vortex ---*/

    // total k, epsilon, and v2 from model
    // NB: Assumes k-e-v2-f model for now
    // TODO: Generalize beyond v2-f
    const su2double ktot = solver[TURB_SOL]->node[iPoint]->GetSolution(0);
    const su2double tdr  = solver[TURB_SOL]->node[iPoint]->GetSolution(1);
    const su2double v2 = solver[TURB_SOL]->node[iPoint]->GetSolution(2);

    // resolved k, from averages
    const su2double kres =
      solver[FLOW_SOL]->average_node[iPoint]->GetResolvedKineticEnergy();


    // ratio of modeled to total TKE
    const su2double alpha = 1.0 - kres/ktot;

    // TODO: Need to limit alpha here?

    // FIXME: Where is average r_M stored?
    // I can't find it... looks like it isn't getting averaged?
    const su2double resolution_adequacy = 1.0;

    const su2double density = solver[FLOW_SOL]->node[iPoint]->GetSolution(0);
    const su2double nu =
      solver[FLOW_SOL]->node[iPoint]->GetLaminarViscosity() / density;

    const su2double Ltot = solver[TURB_SOL]->node[iPoint]->GetTurbLengthscale();
    const su2double Lsgs = alpha*Ltot;

    // FIXME: I think this is equivalent to repo version of CDP,but
    // not consistent with paper description, except for orthogonal
    // grids aligned with coordinate axes.  Check with Sigfried.
    su2double** ResolutionTensor = geometry->node[iPoint]->GetResolutionTensor();
    Lmesh[0] = ResolutionTensor[0][0];
    Lmesh[1] = ResolutionTensor[1][1];
    Lmesh[2] = ResolutionTensor[2][2];

    // Get dwall
    dwall = geometry->node[iPoint]->GetWall_Distance();

    // Compute TG velocity at this point
    this->SetTGField(x, Lsgs, Lmesh, D, dwall, h);

    const su2double Ttot = solver[TURB_SOL]->node[iPoint]->GetTurbTimescale();

    const su2double Ftar = this->GetTargetProduction(v2, Ttot, alpha);

    // Compute PFtest
    su2double PFtest = 0.0;
    for (unsigned short iDim=0; iDim<nDim; iDim++) {
      uf[iDim] = prim_vars[iDim+1] - prim_mean[iDim+1];
      PFtest += uf[iDim]*h[iDim];
    }
    PFtest *= Ftar*dt;

    const su2double Cnu = 1.0;
    const su2double alpha_kol = Cnu*std::sqrt(nu*tdr)/ktot;

    const su2double eta = this->ComputeScalingFactor(Ftar, resolution_adequacy,
                                                     alpha, alpha_kol, PFtest);

    /*--- Store eta*h so we can compute the derivatives ---*/
    for (unsigned short iDim = 0; iDim < nDim; iDim++) {
      node[iPoint][iDim] = eta*h[iDim];
    }

  } // end loop over points
}

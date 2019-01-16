/*!
 * \file forcing_test.cpp
 * \brief Test the forcing terms for the hybrid RANS/LES model
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

#define BOOST_TEST_MODULE Hybrid_Model
#include "boost/test/included/unit_test.hpp"
#include "MPI_global_fixture.hpp"

#include "../include/hybrid_RANS_LES_forcing.hpp"


void WriteQuadMeshFile () {
  /*--- Local variables ---*/
    int KindElem, KindBound, nDim;
    int iElem, iDim, jDim;
    int iNode, jNode, kNode;
    int iPoint;
    int num_Nodes;
    double xSpacing, ySpacing;

    std::ofstream Mesh_File;

    /*--- Set the VTK type for the interior elements and the boundary elements ---*/
    nDim = 2;
    KindElem  = 9; // Quadrilateral
    KindBound = 3; // Line

    /*--- Store the number of nodes in each direction ---*/
    iDim = 4;
    jDim = 4;

    /*--- The grid spacing in each direction ---*/
    xSpacing = 4.0;
    ySpacing = 2.0;

    /*--- Open .su2 grid file ---*/
    Mesh_File.open("test.su2", ios::out);
    Mesh_File << fixed << setprecision(15);

    /*--- Write the dimension of the problem and the number of interior elements ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Problem dimension" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NDIME= 2" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Inner element connectivity" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NELEM= " <<  (iDim-1)*(jDim-1) << std::endl;

    /*--- Write the element connectivity ---*/
    iElem = 0;
      for (jNode = 0; jNode < jDim-1; jNode++) {
        for (iNode = 0; iNode < iDim-1; iNode++) {
          Mesh_File << KindElem << "\t";
          Mesh_File << iNode     + (jNode*jDim)     << "\t";
          Mesh_File << (iNode+1) + (jNode*jDim)     << "\t";
          // NOTE: Reverse ordering here is essential
          Mesh_File << (iNode+1) + ((jNode+1)*jDim) << "\t";
          Mesh_File << iNode     + ((jNode+1)*jDim) << "\t";
          Mesh_File << iElem << std::endl;
          iElem++;
        }
      }


    /*--- Compute the number of nodes and write the node coordinates ---*/
    num_Nodes = iDim*jDim;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Node coordinates" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NPOIN= " << num_Nodes << std::endl;
    iPoint = 0;
      for (jNode = 0; jNode < jDim; jNode++) {
        for (iNode = 0; iNode < iDim; iNode++) {
          Mesh_File << iNode*xSpacing << "\t";
          Mesh_File << jNode*ySpacing << "\t";
          Mesh_File << iPoint << std::endl;
          iPoint++;
        }
      }



    /*--- Write the header information for the boundary markers ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Boundary elements" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NMARK= 4" << std::endl;

    /*--- Write the boundary information for each marker ---*/
    Mesh_File << "MARKER_TAG= lower" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1) << std::endl;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      Mesh_File << KindBound << "\t";
      Mesh_File << iNode       << "\t";
      Mesh_File << (iNode + 1) << std::endl;
    }
    Mesh_File << "MARKER_TAG= right" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1) << std::endl;
    for (jNode = 0; jNode < jDim-1; jNode++) {
      Mesh_File << KindBound << "\t";
      Mesh_File << (jNode+1)*iDim - 1 << "\t";
      Mesh_File << (jNode+2)*iDim - 1 << std::endl;
    }
    Mesh_File << "MARKER_TAG= upper" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1) << std::endl;
    for (iNode = jDim*(iDim-1); iNode < iDim*jDim-1; ++iNode) {
      Mesh_File << KindBound << "\t";
      Mesh_File << iNode       << "\t";
      Mesh_File << (iNode + 1) << std::endl;
    }
    Mesh_File << "MARKER_TAG= left" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1) << std::endl;
    for (jNode = 0; jNode < jDim-1; ++jNode) {
      Mesh_File << KindBound << "\t";
      Mesh_File << (jNode  )*iDim << "\t";
      Mesh_File << (jNode+1)*iDim << std::endl;
    }

    /*--- Close the mesh file and exit ---*/
    Mesh_File.close();
}

void WriteCfgFile(const unsigned short& nDim) {
  std::ofstream cfg_file;

  cfg_file.open("test.cfg", ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  if (nDim == 2)
    cfg_file << "MARKER_FAR= ( lower upper left right )"  << std::endl;
  else
    cfg_file << "MARKER_FAR= ( top bottom back front left right )"  << std::endl;
  cfg_file << "MESH_FILENAME= test.su2" << std::endl;
  cfg_file << "MESH_FORMAT= SU2" << std::endl;

  cfg_file.close();

}

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

BOOST_AUTO_TEST_CASE(Shifting_Returns_Original_Coordinates_When_Time_Is_Zero) {

  const unsigned short nDim = 3;

  CHybridForcingTGSF forcing(nDim, 0, 0);

  su2double time = 0.0;
  su2double T_m = 1.0;
  su2double x_original[nDim], x_shifted, x_exact[nDim];
  su2double velocity[nDim];

  x_original[0] = -1.0; x_original[1] = 1.0; x_original[2] = 2.0;
  for (unsigned short iDim = 0; iDim < nDim; iDim++) velocity[iDim] = 1.0;

  su2double tol = std::numeric_limits<su2double>::epsilon();
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    x_shifted =
        forcing.TransformCoords(x_original[iDim], velocity[iDim], time, T_m);
    BOOST_CHECK_CLOSE(x_shifted, x_original[iDim], tol);
  }
}

BOOST_AUTO_TEST_CASE(Coordinate_Shifting_With_Nonzero_Velocity) {

  const unsigned short nDim = 3;

  CHybridForcingTGSF forcing(nDim, 0, 0);

  const su2double time = 0.5;
  const su2double timescale = 2.0;
  const su2double x_original[nDim] = {-1.0, 1.0, 2.0};
  const su2double x_exact[nDim] =    {-0.5, 1.5, 2.5};
  const su2double velocity = 1;

  su2double tol = std::numeric_limits<su2double>::epsilon();
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    const su2double x_shifted =
        forcing.TransformCoords(x_original[iDim], velocity, time, timescale);
    BOOST_CHECK_CLOSE(x_shifted, x_exact[iDim], tol);
  }
}

BOOST_AUTO_TEST_CASE(Initial_Taylor_Green_Is_Periodic_In_L) {

  const unsigned short nDim = 3;
  const su2double x1[3] = {0.0, 0.0, 0.0};
  const su2double L[3] = {1.0, 2.0, 3.0};

  su2double b_at_x1[3], b_at_x2[3];

  CHybridForcingTGSF forcing(nDim, 0, 0);

  forcing.SetTGField(x1, L, b_at_x1);

  su2double tol = std::numeric_limits<su2double>::epsilon();
  su2double x2[3];
  for (unsigned short iDim = 0; iDim < nDim; iDim++) {
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      x2[jDim] = x1[jDim];
    }
    x2[iDim] += L[iDim];
    forcing.SetTGField(x2, L, b_at_x2);
    for (unsigned short jDim = 0; jDim < nDim; jDim++) {
      BOOST_CHECK_SMALL(std::abs(b_at_x2[jDim] - b_at_x1[jDim]), tol);
    }
  }
}


BOOST_AUTO_TEST_CASE(Smoke_Test_for_Forcing) {

}

BOOST_AUTO_TEST_CASE(Compute_Derivatives) {

  /*--- Arrange ---*/

  const unsigned short nDim = 3;
  WriteCfgFile(nDim);
  CConfig* config = new CConfig("test.cfg", SU2_CFD, 0, 1, 2, VERB_NONE);

  // The use of "geometry_aux" is necessary to duplicate a multigrid
  // configuration
  CGeometry *geometry_aux = new CPhysicalGeometry(config, 0, 1);
  CGeometry* geometry = new CPhysicalGeometry(geometry_aux, config);
  delete geometry_aux;

  // Initialize the geometry
  geometry->SetBoundaries(config);
  geometry->SetPoint_Connectivity();
  geometry->SetElement_Connectivity();
  geometry->SetBoundVolume();
  geometry->Check_IntElem_Orientation(config);
  geometry->Check_BoundElem_Orientation(config);
  geometry->SetEdges();
  geometry->SetVertex(config);
  geometry->SetCoord_CG();
  geometry->SetControlVolume(config, ALLOCATE);
  geometry->SetBoundControlVolume(config, ALLOCATE);

  CHybridForcingTGSF forcing(geometry, config);

  /*--- Act ---*/

  forcing.SetForcing_Gradient_LS(geometry, config);

  /*--- Assert ---*/

}


/*!
 * \file resolution_integration_test.cpp
 * \brief This test checks whether the resolution tensor is correctly set for a grid
 * of quadrilateral cells.
 * \author C. Pederson
 * \version 4.3.0 "Cardinal"
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

#define BOOST_TEST_DYN_LINK
#include <boost/test/unit_test.hpp>

#include "MPI_global_fixture.hpp"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits> // used to find machine epsilon
#include <cmath>  // std::abs

#include "../../Common/include/config_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"

namespace resolution_tensor_test {

/* ----------------------------------------------------------------------------
 *  Functions for Grid Setup
 * --------------------------------------------------------------------------*/


void WriteQuadMeshFile (const char* mesh_filename) {
  /*--- Local variables ---*/
    int KindElem, KindBound, nDim;
    int iElem, iDim, jDim;
    int iNode, jNode;
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
    Mesh_File.open(mesh_filename, ios::out);
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

void WriteTriangleMeshFile (const char* mesh_filename) {
  /*--- Local variables ---*/
    int KindElem, KindBound;
    int iElem, nNode, mNode;
    int iNode, jNode, nPoint, iPoint, jPoint, kPoint;
    ofstream Mesh_File;

    /*--- Set the VTK type for the interior elements and the boundary elements ---*/
    KindElem  = 5; // Triangle
    KindBound = 3; // Line

    /*--- Store the number of nodes and output mesh filename ---*/
    nNode     = 3;
    mNode     = 3;

    /*--- Open .su2 grid file ---*/
    Mesh_File.precision(15);
    Mesh_File.open(mesh_filename, ios::out);

    /*--- Write the dimension of the problem and the number of interior elements ---*/
    Mesh_File << "%" << endl;
    Mesh_File << "% Problem dimension" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NDIME= 2" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "% Inner element connectivity" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NELEM= " <<  2*(nNode-1)*(mNode-1) << endl;

    /*--- Write the element connectivity ---*/
    iElem = 0;
    for (jNode = 0; jNode < mNode-1; jNode++) {
      for (iNode = 0; iNode < nNode-1; iNode++) {
        iPoint = jNode*nNode + iNode;
        jPoint = jNode*nNode + iNode + 1;
        kPoint = (jNode + 1)*nNode + iNode;
        Mesh_File << KindElem << "\t" << iPoint << "\t" << jPoint << "\t" << kPoint << "\t" << iElem << endl;
        iElem ++;
        iPoint = jNode*nNode + (iNode + 1);
        jPoint = (jNode + 1)*nNode + (iNode + 1);
        kPoint = (jNode + 1)*nNode + iNode;
        Mesh_File << KindElem << "\t" << iPoint << "\t" << jPoint << "\t" << kPoint << "\t" << iElem << endl;
        iElem++;
      }
    }

    /*--- Compute the number of nodes and write the node coordinates ---*/
    nPoint = (nNode)*(mNode);
    Mesh_File << "%" << endl;
    Mesh_File << "% Node coordinates" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NPOIN= " << nNode*mNode << endl;
    iPoint = 0;
    for (jNode = 0; jNode < mNode; jNode++) {
      for (iNode = 0; iNode < nNode; iNode++) {
        Mesh_File << ((double)iNode)/((double)(nNode-1)) << "\t" << ((double)jNode)/((double)(mNode-1)) << "\t" << iPoint << endl;
        iPoint++;
      }
    }

    /*--- Write the header information for the boundary markers ---*/
    Mesh_File << "%" << endl;
    Mesh_File << "% Boundary elements" << endl;
    Mesh_File << "%" << endl;
    Mesh_File << "NMARK= 4" << endl;

    /*--- Write the boundary information for each marker ---*/
    Mesh_File << "MARKER_TAG= lower" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (nNode-1) << endl;
    for (iNode = 0; iNode < nNode-1; iNode++) {
      Mesh_File << KindBound << "\t" << iNode << "\t" << (iNode + 1) << endl;
    }
    Mesh_File << "MARKER_TAG= right" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (mNode-1) << endl;
    for (jNode = 0; jNode < mNode-1; jNode++) {
      Mesh_File << KindBound << "\t" << jNode*nNode + (nNode - 1) << "\t" << (jNode + 1)*nNode + (nNode - 1) << endl;
    }
    Mesh_File << "MARKER_TAG= upper" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (nNode-1) << endl;
    for (iNode = 0; iNode < nNode-1; iNode++) {
      Mesh_File << KindBound << "\t" << (nNode*mNode - 1) - iNode << "\t" << (nNode*mNode - 1) - (iNode + 1) << endl;
    }
    Mesh_File << "MARKER_TAG= left" << endl;
    Mesh_File << "MARKER_ELEMS= "<< (mNode-1) << endl;
    for (jNode = mNode-2; jNode > mNode-4; jNode--) {
      Mesh_File << KindBound << "\t" << (jNode + 1)*nNode << "\t" << jNode*nNode << endl;
    }

    /*--- Close the mesh file and exit ---*/
    Mesh_File.close();
}

void WriteHexMeshFile (su2double delta_x, su2double delta_y, su2double delta_z, const char* mesh_filename) {
    /*--- Local variables ---*/
    int KindElem, KindBound, nDim;
    int iElem, iDim, jDim, kDim;
    int iNode, jNode, kNode;
    int iPoint;
    int num_Nodes;

    std::ofstream Mesh_File;

    /*--- Set the VTK type for the interior elements and the boundary elements ---*/
    nDim = 3;
    KindElem  = 12; // Hexahedra
    KindBound = 9; // Quadrilateral

    /*--- Store the number of nodes in each direction ---*/
    iDim = 4;
    jDim = 4;
    kDim = 4;

    /*--- Open .su2 grid file ---*/
    Mesh_File.open(mesh_filename, ios::out);
    Mesh_File << fixed << setprecision(15);

    /*--- Write the dimension of the problem and the number of interior elements ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Problem dimension" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NDIME= 3" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Inner element connectivity" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NELEM= " <<  (iDim-1)*(jDim-1)*(kDim-1) << std::endl;

    /*--- Write the element connectivity ---*/
    iElem = 0;
    for (kNode = 0; kNode < kDim-1; kNode++) {
      for (jNode = 0; jNode < jDim-1; jNode++) {
        for (iNode = 0; iNode < iDim-1; iNode++) {
          Mesh_File << KindElem << "\t";
          // Proper ordering here is essential.
          // See VTK documentation for hexahedral cells for ordering.
          Mesh_File << iNode     + (jNode*jDim)     + (kNode*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + (jNode*jDim)     + (kNode*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + ((jNode+1)*jDim) + (kNode*jDim*iDim) << "\t";
          Mesh_File << iNode     + ((jNode+1)*jDim) + (kNode*jDim*iDim) << "\t";
          Mesh_File << iNode     + (jNode*jDim)     + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + (jNode*jDim)     + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << (iNode+1) + ((jNode+1)*jDim) + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << iNode     + ((jNode+1)*jDim) + ((kNode+1)*jDim*iDim) << "\t";
          Mesh_File << iElem << std::endl;
          iElem++;
        }
      }
    }


    /*--- Compute the number of nodes and write the node coordinates ---*/
    num_Nodes = iDim*jDim*kDim;
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Node coordinates" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NPOIN= " << num_Nodes << std::endl;
    iPoint = 0;
    for (kNode = 0; kNode < kDim; kNode++) {
      for (jNode = 0; jNode < jDim; jNode++) {
        for (iNode = 0; iNode < iDim; iNode++) {
          Mesh_File << iNode*delta_x << "\t";
          Mesh_File << jNode*delta_y << "\t";
          Mesh_File << kNode*delta_z << "\t";
          Mesh_File << iPoint << std::endl;
          iPoint++;
        }
      }
    }



    /*--- Write the header information for the boundary markers ---*/
    Mesh_File << "%" << std::endl;
    Mesh_File << "% Boundary elements" << std::endl;
    Mesh_File << "%" << std::endl;
    Mesh_File << "NMARK= 6" << std::endl;

    /*--- Write the boundary information for each marker ---*/
    Mesh_File << "MARKER_TAG= bottom" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1)*(jDim-1) << std::endl;
    kNode = 0;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (jNode = 0; jNode < jDim-1; jNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + (jNode+1)*jDim + kNode*(iDim*jDim) << "\t";
        Mesh_File << iNode     + (jNode+1)*jDim + kNode*(iDim*jDim) << std::endl;
      }
    }
    Mesh_File << "MARKER_TAG= top" << std::endl;
    Mesh_File << "MARKER_ELEMS= " << (iDim-1)*(jDim-1) << std::endl;
    kNode = kDim-1;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (jNode = 0; jNode < jDim-1; jNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + jNode*jDim     + kNode*(iDim*jDim) << "\t";
        Mesh_File << (iNode+1) + (jNode+1)*jDim + kNode*(iDim*jDim) << "\t";
        Mesh_File << iNode     + (jNode+1)*jDim + kNode*(iDim*jDim) << std::endl;
      }
    }

    Mesh_File << "MARKER_TAG= left" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1)*(kDim-1) << std::endl;
    iNode = 0;
    for (jNode = 0; jNode < jDim-1; jNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode + jNode*jDim     + kNode*(iDim*jDim)     << "\t";
        Mesh_File << iNode + jNode*jDim     + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + kNode*(iDim*jDim)     << std::endl;
      }
    }
    Mesh_File << "MARKER_TAG= right" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (jDim-1)*(kDim-1) << std::endl;
    iNode = iDim-1;
    for (jNode = 0; jNode < jDim-1; jNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode + jNode*jDim     + kNode*(iDim*jDim)     << "\t";
        Mesh_File << iNode + jNode*jDim     + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode + (jNode+1)*jDim + kNode*(iDim*jDim)     << std::endl;
      }
    }

    Mesh_File << "MARKER_TAG= front" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1)*(kDim-1) << std::endl;
    jNode = 0;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode     + jNode*jDim + (kNode+1)*(iDim*jDim) << std::endl;
      }
    }
    Mesh_File << "MARKER_TAG= back" << std::endl;
    Mesh_File << "MARKER_ELEMS= "<< (iDim-1)*(kDim-1) << std::endl;
    jNode = jDim-1;
    for (iNode = 0; iNode < iDim-1; iNode++) {
      for (kNode = 0; kNode < kDim-1; kNode++) {
        Mesh_File << KindBound << "\t";
        Mesh_File << iNode     + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + kNode*(iDim*jDim)     << "\t";
        Mesh_File << (iNode+1) + jNode*jDim + (kNode+1)*(iDim*jDim) << "\t";
        Mesh_File << iNode     + jNode*jDim + (kNode+1)*(iDim*jDim) << std::endl;
      }
    }
    /*--- Close the mesh file and exit ---*/
    Mesh_File.close();
}

void WriteCfgFile(unsigned short nDim, const char* cfg_filename,
                  const char* mesh_filename) {
  std::ofstream cfg_file;

  cfg_file.open(cfg_filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "HYBRID_RANSLES= MODEL_SPLIT" << std::endl;
  cfg_file << "RUNTIME_AVERAGING= POINTWISE" << std::endl;
  cfg_file << "UNSTEADY_SIMULATION= TIME_STEPPING" << std::endl;
  if (nDim == 2)
    cfg_file << "MARKER_FAR= ( lower upper left right )"  << std::endl;
  else
    cfg_file << "MARKER_FAR= ( top bottom back front left right )"  << std::endl;
  cfg_file << "MESH_FILENAME= " << mesh_filename << std::endl;
  cfg_file << "MESH_FORMAT= SU2" << std::endl;

  cfg_file.close();

}

/* ----------------------------------------------------------------------------
 *  Test Fixtures
 * --------------------------------------------------------------------------*/

struct ResolutionFixture {
  ResolutionFixture()
    : machine_eps(std::numeric_limits<su2double>::epsilon()) { }

  ~ResolutionFixture() {
    delete geometry;
    delete config;
  }

  void SetupConfig(unsigned short nDim, const char* mesh_filename) {
    char cfg_filename[100] = "resolution_tensor_test.cfg";
    WriteCfgFile(nDim, cfg_filename, mesh_filename);
    config = new CConfig(cfg_filename, SU2_CFD, 0, 1, 2, VERB_NONE);
    std::remove(cfg_filename);
  }

  void SetupGeometry() {
    // The use of "geometry_aux" is necessary to duplicate a multigrid
    // configuration
    CGeometry *geometry_aux = new CPhysicalGeometry(config, 0, 1);
    geometry = new CPhysicalGeometry(geometry_aux, config);
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
  }

  static unsigned short cfg_file_counter;
  const su2double machine_eps;
  CConfig* config;
  CGeometry* geometry;
};

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

BOOST_AUTO_TEST_SUITE(ResolutionTensorTest);

BOOST_FIXTURE_TEST_CASE(Triangles_Test, ResolutionFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 2;
  WriteTriangleMeshFile("triangle_test.su2");
  SetupConfig(nDim, "triangle_test.su2");
  SetupGeometry();

  unsigned short iPoint;

  geometry->SetResolutionTensor(config);

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    const su2double* const* Mij = geometry->node[iPoint]->GetResolutionTensor();

    // Check that the values of Mij are correct
    BOOST_CHECK_SMALL(Mij[0][0] - 0.5, machine_eps);
    BOOST_CHECK_SMALL(Mij[1][0] - 0.0, machine_eps);
    BOOST_CHECK_SMALL(Mij[0][1] - 0.0, machine_eps);
    BOOST_CHECK_SMALL(Mij[1][1] - 0.5, machine_eps);
  }

  std::remove("triangle_test.su2");
}

BOOST_FIXTURE_TEST_CASE(Quads_Test, ResolutionFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 2;
  WriteQuadMeshFile("quad_test.su2");
  SetupConfig(nDim, "quad_test.su2");
  SetupGeometry();

  unsigned short iPoint;

  geometry->SetResolutionTensor(config);

  for (iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    const su2double* const* Mij = geometry->node[iPoint]->GetResolutionTensor();

    // Check that the values of Mij are correct
    BOOST_CHECK_SMALL(Mij[0][0] - 4.0, machine_eps);
    BOOST_CHECK_SMALL(Mij[1][0] - 0.0, machine_eps);
    BOOST_CHECK_SMALL(Mij[0][1] - 0.0, machine_eps);
    BOOST_CHECK_SMALL(Mij[1][1] - 2.0, machine_eps);
  }
  std::remove("quad_test.su2");
}


BOOST_FIXTURE_TEST_CASE(Hexahedra, ResolutionFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 3;
  WriteHexMeshFile(3,2,1, "hex_mesh.su2");
  SetupConfig(nDim, "hex_mesh.su2");
  SetupGeometry();

  geometry->SetResolutionTensor(config);

  for (unsigned long iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    const su2double* const* Mij = geometry->node[iPoint]->GetResolutionTensor();

    // Build the test info
    std::stringstream msg;
    msg << "Computed array elements:" << std::endl;
    msg << "[[";
    msg << Mij[0][0] << "," << Mij[0][1] << "," << Mij[0][2] << "],[";
    msg << Mij[1][0] << "," << Mij[1][1] << "," << Mij[1][2] << "],[";
    msg << Mij[2][0] << "," << Mij[2][1] << "," << Mij[2][2] << "]]";
    BOOST_TEST_CONTEXT(msg.str()) {
      // Check that the values of Mij are correct
      BOOST_CHECK_CLOSE_FRACTION(Mij[0][0], 3.0, machine_eps);
      BOOST_CHECK_SMALL(Mij[0][1], machine_eps);
      BOOST_CHECK_SMALL(Mij[0][2], machine_eps);
      BOOST_CHECK_SMALL(Mij[1][0], machine_eps);
      BOOST_CHECK_CLOSE_FRACTION(Mij[1][1], 2.0, machine_eps);
      BOOST_CHECK_SMALL(Mij[1][2], machine_eps);
      BOOST_CHECK_SMALL(Mij[2][0], machine_eps);
      BOOST_CHECK_SMALL(Mij[2][1], machine_eps);
      BOOST_CHECK_CLOSE_FRACTION(Mij[2][2], 1.0, machine_eps);
    }
  }

  std::remove("hex_mesh.su2");
}

BOOST_FIXTURE_TEST_CASE(M43_Power, ResolutionFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 3;
  WriteHexMeshFile(3, 2, 1, "hex_mesh_2.cfg");
  SetupConfig(nDim, "hex_mesh_2.cfg");
  SetupGeometry();

  geometry->SetResolutionTensor(config);

  for (unsigned long iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    const su2double* const* M43 = geometry->node[iPoint]->GetResolutionTensor43();

    // Build the test info
    std::stringstream msg;
    msg << "Computed array elements:" << std::endl;
    msg << "[[";
    msg << M43[0][0] << "," << M43[0][1] << "," << M43[0][2] << "],[";
    msg << M43[1][0] << "," << M43[1][1] << "," << M43[1][2] << "],[";
    msg << M43[2][0] << "," << M43[2][1] << "," << M43[2][2] << "]]";
    BOOST_TEST_CONTEXT(msg.str()) {
      // Check that the values of Mij are correct
      BOOST_CHECK_CLOSE_FRACTION(M43[0][0], std::pow(3.0, 4.0/3), machine_eps);
      BOOST_CHECK_SMALL(M43[0][1], machine_eps);
      BOOST_CHECK_SMALL(M43[0][2], machine_eps);
      BOOST_CHECK_SMALL(M43[1][0], machine_eps);
      BOOST_CHECK_CLOSE_FRACTION(M43[1][1], std::pow(2.0, 4.0/3), machine_eps);
      BOOST_CHECK_SMALL(M43[1][2], machine_eps);
      BOOST_CHECK_SMALL(M43[2][0], machine_eps);
      BOOST_CHECK_SMALL(M43[2][1], machine_eps);
      BOOST_CHECK_CLOSE_FRACTION(M43[2][2], 1.0, machine_eps);
    }
  }

  std::remove("hex_mesh_2.su2");
}

BOOST_FIXTURE_TEST_CASE(ResolutionConstantEqualsOneForIsotropic, ResolutionFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 3;
  WriteHexMeshFile(1, 1, 1, "hex_mesh_isotropic.su2");
  SetupConfig(nDim, "hex_mesh_isotropic.su2");
  SetupGeometry();

  geometry->SetResolutionTensor(config);

  for (unsigned long iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    const su2double C_M = geometry->node[iPoint]->GetResolutionCoeff();

    // Check that the values of Mij are correct
    // The constants used for the fit have this precision
    const su2double tolerance = 0.00000000001;
    BOOST_CHECK_SMALL(C_M - 1, tolerance);
  }

  std::remove("hex_mesh_isotropic.su2");
}

BOOST_FIXTURE_TEST_CASE(ResolutionConstantForAnisotropic, ResolutionFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 3;
  WriteHexMeshFile(3, 2, 1, "hex_mesh_aniso.su2");
  SetupConfig(nDim, "hex_mesh_aniso.su2");
  SetupGeometry();

  geometry->SetResolutionTensor(config);

  for (unsigned long iPoint = 0; iPoint < geometry->GetnPointDomain(); iPoint++) {

    const su2double C_M = geometry->node[iPoint]->GetResolutionCoeff();

    // Check value against hand-calculated value
    const su2double correct_value = 1.04805425805;
    // The constants used for the fit have this precision
    const su2double tolerance = 0.00000000001;
    BOOST_CHECK_CLOSE_FRACTION(C_M, correct_value, tolerance);
  }
  std::remove("hex_mesh_aniso.su2");
}

BOOST_AUTO_TEST_SUITE_END();

} // end namespace

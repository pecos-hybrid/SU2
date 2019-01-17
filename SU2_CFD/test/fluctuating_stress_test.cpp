/*!
 * \file fluctuating_stress_test.cpp
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

#define BOOST_TEST_MODULE FluctuatingStress
#include "boost/test/included/unit_test.hpp"

#include "../include/fluctuating_stress.hpp"

#include "../../Common/test/MPI_global_fixture.hpp"

#include <cstdio> // std::remove
#include <fstream>

BOOST_GLOBAL_FIXTURE( MPIGlobalFixture );

void WriteHexMeshFile (su2double delta_x, su2double delta_y, su2double delta_z) {
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
    Mesh_File.open("test.su2", ios::out);
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


void WriteCfgFile(unsigned short nDim, const char* filename) {
  std::ofstream cfg_file;

  cfg_file.open(filename, ios::out);
  cfg_file << "PHYSICAL_PROBLEM= NAVIER_STOKES" << std::endl;
  cfg_file << "HYBRID_RANSLES= MODEL_SPLIT" << std::endl;
  cfg_file << "RUNTIME_AVERAGING= POINTWISE" << std::endl;
  cfg_file << "UNSTEADY_SIMULATION= TIME_STEPPING" << std::endl;
  cfg_file << "KIND_TURB_MODEL= KE" << std::endl;
  if (nDim == 2)
    cfg_file << "MARKER_FAR= ( lower upper left right )"  << std::endl;
  else
    cfg_file << "MARKER_FAR= ( top bottom back front left right )"  << std::endl;
  cfg_file << "MESH_FILENAME= test.su2" << std::endl;
  cfg_file << "MESH_FORMAT= SU2" << std::endl;

  cfg_file.close();

}
/* ----------------------------------------------------------------------------
 *  Test Fixtures
 * --------------------------------------------------------------------------*/

struct FluctuatingStressFixture {
  FluctuatingStressFixture()
    : machine_eps(std::numeric_limits<su2double>::epsilon()) { }

  ~FluctuatingStressFixture() {
    delete geometry;
    delete config;
  }

  void SetupConfig(unsigned short nDim) {
    char cfg_filename[100] = "fluctuating_stress_test.cfg";
    WriteCfgFile(nDim, cfg_filename);
    config = new CConfig(cfg_filename, SU2_CFD, 0, 1, 2, VERB_NONE);
    std::remove(cfg_filename);
  }

  void SetupGeometry() {
    // The use of "geometry_aux" is necessary to imitate a multigrid configuration
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
    geometry->SetResolutionTensor();
  }

  const su2double machine_eps;
  CConfig* config;
  CGeometry* geometry;
};

/* ----------------------------------------------------------------------------
 *  Tests
 * --------------------------------------------------------------------------*/

BOOST_FIXTURE_TEST_CASE(Smoke_Test, FluctuatingStressFixture) {

  // Write out the mesh and configuration files.
  const unsigned short nDim = 3;
  WriteHexMeshFile(1, 1, 1);
  SetupConfig(nDim);
  SetupGeometry();

  CFluctuatingStress* m43 = new CM43Model(nDim, config);

  su2double** eddy_viscosity = new su2double*[nDim];
  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    eddy_viscosity[iDim] = new su2double[nDim];
  unsigned long dummy_point = 0;

  const su2double flow_prim_vars[] = {0.0, 0.0, 0.0, 0.0, 0.0, 1.0};
  const su2double turb_vars[] = {0.0, 1.0, 0.0, 0.0};
  m43->SetPrimitive(flow_prim_vars);
  m43->SetTurbVar(turb_vars);
  const su2double rans_eddy_viscosity = 0.0;

  /*--- Test ---*/

  m43->CalculateEddyViscosity(geometry, config, dummy_point,
                              rans_eddy_viscosity, eddy_viscosity);


  for (unsigned short iDim = 0; iDim < nDim; iDim++)
    delete [] eddy_viscosity[iDim];
  delete [] eddy_viscosity;
  delete m43;
}

/*!
 * \file primal_grid_structure.cpp
 * \brief Main classes for defining the primal grid elements
 * \author F. Palacios
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
 * Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
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

#include "../include/primal_grid_structure.hpp"

unsigned short CPrimalGrid::nDim;

CPrimalGrid::CPrimalGrid(void) {
  
  /*--- Set the default values for the pointers ---*/
  Nodes = NULL;
  Neighbor_Elements = NULL;
  ElementOwnsFace = NULL;
  PeriodIndexNeighbors = NULL;
  Coord_CG = NULL;
  Coord_FaceElems_CG = NULL;
  JacobianFaceIsConstant = NULL;
  GlobalIndex = 0;
  ResolutionTensor = NULL;
  ResolutionValues = NULL;
  ResolutionVectors = NULL;
  
}

CPrimalGrid::~CPrimalGrid() {
  unsigned short iDim;

 if (Nodes != NULL) delete[] Nodes;
 if (Coord_CG != NULL) delete[] Coord_CG;
 if (Neighbor_Elements != NULL) delete[] Neighbor_Elements;
 if (ElementOwnsFace != NULL) delete[] ElementOwnsFace;
 if (PeriodIndexNeighbors != NULL) delete[] PeriodIndexNeighbors;
 if (JacobianFaceIsConstant != NULL) delete[] JacobianFaceIsConstant;
  if (ResolutionValues != NULL) delete[] ResolutionValues;
  if (ResolutionVectors != NULL) {
    for (iDim = 0; iDim < nDim; iDim++) delete [] ResolutionVectors[iDim];
    delete [] ResolutionVectors;
  }
  if (ResolutionTensor != NULL) {
    for (iDim = 0; iDim < nDim; iDim++) delete [] ResolutionTensor[iDim];
    delete [] ResolutionTensor;
  }
}

void CPrimalGrid::SetCoord_CG(su2double **val_coord) {
	unsigned short iDim, iNode, NodeFace, iFace;
	
  AD::StartPreacc();
  AD::SetPreaccIn(val_coord, GetnNodes(), nDim);

	for (iDim = 0; iDim < nDim; iDim++) {
		Coord_CG[iDim] = 0.0;
		for (iNode = 0; iNode < GetnNodes();  iNode++)
			Coord_CG[iDim] += val_coord[iNode][iDim]/su2double(GetnNodes());
	}
	
	for (iFace = 0; iFace < GetnFaces();  iFace++)
		for (iDim = 0; iDim < nDim; iDim++) {
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
			for (iNode = 0; iNode < GetnNodesFace(iFace); iNode++) {
				NodeFace = GetFaces(iFace, iNode);
				Coord_FaceElems_CG[iFace][iDim] += val_coord[NodeFace][iDim]/su2double(GetnNodesFace(iFace));
			}
		}

  AD::SetPreaccOut(Coord_CG, nDim);
  AD::SetPreaccOut(Coord_FaceElems_CG, GetnFaces(), nDim);
  AD::EndPreacc();

}

void CPrimalGrid::GetAllNeighbor_Elements() {
	cout << "( ";
	for (unsigned short iFace = 0; iFace < GetnFaces(); iFace++)
	{
		cout << GetNeighbor_Elements(iFace) << ", ";
	}
	cout << ")"  << endl;
}

void CPrimalGrid::SetResolutionTensor(su2double** val_coord) {
  unsigned short iDim, jDim;
  unsigned short iNode;
  su2double* coord_max = new su2double[nDim];
  su2double* coord_min = new su2double[nDim];

  /*--- Allocate Resolution Tensor ---*/
  ResolutionTensor = new su2double* [nDim];
  ResolutionValues = new su2double[nDim];
  ResolutionVectors = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    ResolutionTensor[iDim] = new su2double [nDim];
    ResolutionVectors[iDim] = new su2double[nDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      ResolutionTensor[iDim][jDim] = 0.0;
    }
  }

  /*--- Initialize the maximum and minimum as the CG ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    coord_max[iDim] = GetCG(iDim);
    coord_min[iDim] = GetCG(iDim);
  }

  /*--- Find the maximum and minimum x,y,z values ---*/
  for (iNode = 0; iNode < GetnNodes(); iNode++) {
    for (iDim = 0; iDim < nDim; iDim++) {
      if (val_coord[iNode][iDim] > coord_max[iDim]) {
        coord_max[iDim] = val_coord[iNode][iDim];
      } else if (val_coord[iNode][iDim] < coord_min[iDim]) {
        coord_min[iDim] = val_coord[iNode][iDim];
      }
    }
  }

  /*--- Set the elements of the resolution tensor ---*/
  for (iDim = 0; iDim < nDim; iDim++) {
    for (jDim = 0; jDim < nDim; jDim++) {
      if (iDim == jDim) ResolutionTensor[iDim][jDim] =
          coord_max[iDim] - coord_min[iDim];
    }
  }

  /*--- Record calculation values in values and vectors ---*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    ResolutionValues[iDim] = coord_max[iDim] - coord_min[iDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      ResolutionVectors[iDim][jDim] = (iDim == jDim);
    }
  }

  delete[] coord_max;
  delete[] coord_min;
}

void CPrimalGrid::GramSchmidt(std::vector<std::vector<su2double> > &w,
                              std::vector<std::vector<su2double> > &v) {
  unsigned short iDim, jDim;
  const unsigned short nDim = w.size();

  vector<su2double> magnitude;
  for (iDim =0; iDim < nDim; iDim++) {
    magnitude.push_back(inline_magnitude(w[iDim]));
  }

  su2double max = magnitude[0], min = magnitude[0];
  unsigned short maxloc = 0, minloc = 0, medianloc=0;
  for (iDim = 1; iDim < nDim; iDim++) {
    if (magnitude[iDim] > max) {
      max = magnitude[iDim]; maxloc = iDim;
    } else if (magnitude[iDim] <= min) {
      min = magnitude[iDim]; minloc = iDim;
    }
  }
  if (nDim == 3) {
    for (iDim = 0; iDim < nDim; iDim++) {
      if (iDim != maxloc && iDim != minloc) {
        medianloc = iDim; break;
      }
    }
  } else {
    medianloc = minloc;
  }

  /*-- Set the largest basis vector to the first input vector --*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    v[maxloc][iDim] = w[maxloc][iDim];
  }

  if (nDim > 1) {
    /*-- Compute the next orthogonal vector --*/
    for (iDim = 0; iDim < nDim; ++iDim) {
      v[medianloc][iDim] = w[medianloc][iDim] -
          inline_dot_prod(w[medianloc],v[maxloc])/
          inline_dot_prod(v[maxloc],v[maxloc])*v[maxloc][iDim];
    }
  }

  if (nDim > 2) {
    /*-- Compute the next orthogonal vector --*/
    for (iDim = 0; iDim < nDim; ++iDim) {
      v[minloc][iDim] = w[minloc][iDim] -
          inline_dot_prod(w[minloc],v[maxloc])/
          inline_dot_prod(v[maxloc],v[maxloc])*v[maxloc][iDim] -
          inline_dot_prod(w[minloc],v[medianloc])/
          inline_dot_prod(v[medianloc],v[medianloc])*v[medianloc][iDim];
    }
  }

  /*-- Normalize results of the Gram-Schmidt --*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    magnitude[iDim] = inline_magnitude(v[iDim]);
    for (jDim = 0; jDim < nDim; ++jDim) {
      v[iDim][jDim] /= magnitude[iDim];
    }
  }
}

void CPrimalGrid::InitializeJacobianConstantFaces(unsigned short val_nFaces) {

  /*--- Allocate the memory for JacobianFaceIsConstant and initialize
        its values to false.     ---*/
  JacobianFaceIsConstant = new bool[val_nFaces];
  for(unsigned short i=0; i<val_nFaces; ++i)
    JacobianFaceIsConstant[i] = false;
}

void CPrimalGrid::InitializeNeighbors(unsigned short val_nFaces) {

  /*--- Allocate the memory for Neighbor_Elements and PeriodIndexNeighbors and
        initialize the arrays to -1 to indicate that no neighbor is present and
        that no periodic transformation is needed to the neighbor. ---*/
  Neighbor_Elements    = new long[val_nFaces];
  ElementOwnsFace      = new bool[val_nFaces];
  PeriodIndexNeighbors = new short[val_nFaces];

  for(unsigned short i=0; i<val_nFaces; ++i) {
    Neighbor_Elements[i]    = -1;
    ElementOwnsFace[i]      =  false;
    PeriodIndexNeighbors[i] = -1;
  }
}

unsigned short CVertexMPI::nFaces = 0;

unsigned short CVertexMPI::nNodes = 1;

unsigned short CVertexMPI::nNeighbor_Elements = 0;

unsigned short CVertexMPI::VTK_Type = 1;

unsigned short CVertexMPI::maxNodesFace = 0;

CVertexMPI::CVertexMPI(unsigned long val_point, unsigned short val_nDim) : CPrimalGrid() {
	unsigned short iDim;
	
	/*--- Allocate CG coordinates ---*/
	nDim = val_nDim;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++) Coord_CG[iDim] = 0.0;
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point;
	
	/*--- By default, no rotation in the solution ---*/
	Rotation_Type = 0;
	
}

CVertexMPI::~CVertexMPI() {
  unsigned short iFaces;
  
    for (iFaces = 0; iFaces < nFaces; iFaces++)
      if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
    if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;

}

void CVertexMPI::Change_Orientation(void) { }

unsigned short CLine::Faces[1][2]={{0,1}};

unsigned short CLine::Neighbor_Nodes[2][1]={{1},{0}};

unsigned short CLine::nNodesFace[1]={2};

unsigned short CLine::nNeighbor_Nodes[2]={1,1};

unsigned short CLine::nFaces = 1;

unsigned short CLine::nNodes = 2;

unsigned short CLine::nNeighbor_Elements = 1;

unsigned short CLine::VTK_Type = 3;

unsigned short CLine::maxNodesFace = 2;

CLine::CLine(unsigned long val_point_0, unsigned long val_point_1,
             unsigned short val_nDim) : CPrimalGrid() {
	unsigned short iDim, iFace;

	/*--- Allocate CG coordinates ---*/
  
	nDim = val_nDim;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
  
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
  
}

CLine::~CLine() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;

}

void CLine::Change_Orientation(void) {
	unsigned long Point_0, Point_1;
  
	Point_0 = Nodes[0];
	Point_1 = Nodes[1];
	Nodes[0] = Point_1;
	Nodes[1] = Point_0;

}

unsigned short CTriangle::Faces[3][2] = {{0,1},{1,2},{2,0}};

unsigned short CTriangle::Neighbor_Nodes[3][2] = {{1,2},{2,0},{0,1}};

unsigned short CTriangle::nNodesFace[3] = {2,2,2};

unsigned short CTriangle::nNeighbor_Nodes[3] = {2,2,2};

unsigned short CTriangle::nFaces = 3;

unsigned short CTriangle::nNodes = 3;

unsigned short CTriangle::nNeighbor_Elements = 3;

unsigned short CTriangle::VTK_Type = 5;

unsigned short CTriangle::maxNodesFace = 2;

CTriangle::CTriangle(unsigned long val_point_0, unsigned long val_point_1,
					 unsigned long val_point_2, unsigned short val_nDim) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = val_nDim;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CTriangle::~CTriangle() {
  unsigned short iFaces;

  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CTriangle::Change_Orientation(void) {
	unsigned long Point_0, Point_2;

	Point_0 = Nodes[0];
	Point_2 = Nodes[2];
	Nodes[0] = Point_2;
	Nodes[2] = Point_0;
  
}

unsigned short CQuadrilateral::Faces[4][2] = {{0,1},{1,2},{2,3},{3,0}};

unsigned short CQuadrilateral::Neighbor_Nodes[4][2] = {{1,3},{2,0},{3,1},{0,2}};

unsigned short CQuadrilateral::nNodesFace[4] = {2,2,2,2};

unsigned short CQuadrilateral::nNeighbor_Nodes[4] = {2,2,2,2};

unsigned short CQuadrilateral::nFaces = 4;

unsigned short CQuadrilateral::nNodes = 4;

unsigned short CQuadrilateral::nNeighbor_Elements = 4;

unsigned short CQuadrilateral::VTK_Type = 9;

unsigned short CQuadrilateral::maxNodesFace = 2;

CQuadrilateral::CQuadrilateral(unsigned long val_point_0, unsigned long val_point_1,
					   unsigned long val_point_2, unsigned long val_point_3, unsigned short val_nDim) 
: CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = val_nDim;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	
	
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CQuadrilateral::~CQuadrilateral() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CQuadrilateral::SetResolutionTensor(su2double **val_coord) {

  unsigned short iDim, jDim, kDim, lDim;
  unsigned short iFace;

  /*-- Allocate ResolutionTensor --*/
  ResolutionTensor = new su2double* [nDim];
  ResolutionValues = new su2double[nDim];
  ResolutionVectors = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    ResolutionTensor[iDim] = new su2double [nDim];
    ResolutionVectors[iDim] = new su2double [nDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      ResolutionTensor[iDim][jDim] = 0.0;
      ResolutionVectors[iDim][jDim] = 0.0;
    }
  }

  /*--- paired_faces is used to sort the faces into matching pairs.
      The code will look for pairs of faces that are mostly opposite, then
      sort them so that the face indices in paired_faces[0] and paired_faces[1]
      match, then paired_faces[2] and paired_faces[3] match, etc. ---*/
  unsigned short* paired_faces = new unsigned short[nFaces];
  vector<vector<su2double> > eigvecs(nDim, vector<su2double>(nDim));

  /*-- Create cell center to face vectors --*/
  vector<vector<su2double> > center2face(nFaces, vector<su2double>(nDim));
  for (iFace = 0; iFace < nFaces; ++iFace) {
    for (iDim = 0; iDim < nDim; ++iDim) {
      center2face[iFace][iDim] = Coord_FaceElems_CG[iFace][iDim] - Coord_CG[iDim];
    }
  }

  /*-- First vector --*/

  /*-- Choose vector 1 as our first vector to pair up --*/
  paired_faces[0] = 0;

  /*-- Find vector mostly parallel to first --*/
  su2double min_dp = 1.0;
  su2double current_dp;
  for (iFace = 1; iFace < nFaces; ++iFace) {
    current_dp = inline_dot_prod(center2face[0],center2face[iFace])
              /(inline_magnitude(center2face[0])*
                  inline_magnitude(center2face[iFace]));
    if (current_dp < min_dp) {
      min_dp = current_dp;
      paired_faces[1] = iFace;
    }
  }

  for (iDim = 0; iDim < nDim; ++iDim) {
    eigvecs[0][iDim] = Coord_FaceElems_CG[0][iDim] -
        Coord_FaceElems_CG[paired_faces[1]][iDim];
  }

  /*-- second vector --*/

  paired_faces[2] = 0;
  paired_faces[3] = 0;
  for (iFace = 1; iFace < nFaces; ++iFace) {
    if (iFace != paired_faces[1]) {
      if (paired_faces[2] == 0) {
        paired_faces[2] = iFace;
      } else {
        paired_faces[3] = iFace;
      }
    }
  }

  for (iDim = 0; iDim < nDim; ++iDim) {
    eigvecs[1][iDim] = Coord_FaceElems_CG[paired_faces[2]][iDim] -
        Coord_FaceElems_CG[paired_faces[3]][iDim];
  }

  /*-- Get magnitudes --*/
  su2double** eigvalues = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; ++iDim) {
    eigvalues[iDim] = new su2double[nDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      eigvalues[iDim][jDim] = 0.0;
    }
    eigvalues[iDim][iDim] = inline_magnitude(eigvecs[iDim]);
    for (jDim = 0; jDim < nDim; ++jDim) {
      eigvecs[iDim][jDim] /= eigvalues[iDim][iDim];
    }
  }

  /*-- Gram-Schmidt Process to make the vectors orthogonal --*/
  vector<vector<su2double> > temp_eigvecs = eigvecs;
  GramSchmidt(temp_eigvecs, eigvecs);

  /*-- Perform matrix multiplication --*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      for (kDim = 0; kDim < nDim; ++kDim) {
        for (lDim = 0; lDim < nDim; ++lDim) {
          ResolutionTensor[iDim][jDim] += eigvecs[kDim][iDim]*
              eigvalues[kDim][lDim]*eigvecs[lDim][jDim];
        }
      }
    }
  }

  /*--- Record calculation values in values and vectors ---*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    ResolutionValues[iDim] = eigvalues[iDim][iDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      ResolutionVectors[iDim][jDim] = eigvecs[iDim][jDim];
    }
  }

  delete [] paired_faces;
  for (iDim = 0; iDim < nDim; ++iDim) {
    delete [] eigvalues[iDim];
  }
  delete [] eigvalues;

}

void CQuadrilateral::Change_Orientation(void) {
	unsigned long Point_1, Point_3;

	Point_1 = Nodes[1];
	Point_3 = Nodes[3];
	Nodes[1] = Point_3;
	Nodes[3] = Point_1;
  
}

unsigned short CTetrahedron::Faces[4][3]={{0,2,1},{0,1,3},{0,3,2},{1,2,3}};

unsigned short CTetrahedron::Neighbor_Nodes[4][3]={{1,2,3},{0,2,3},{0,1,3},{0,1,2}};

unsigned short CTetrahedron::nNodesFace[4]={3,3,3,3};

unsigned short CTetrahedron::nNeighbor_Nodes[4]={3,3,3,3};

unsigned short CTetrahedron::nFaces = 4;

unsigned short CTetrahedron::nNodes = 4;

unsigned short CTetrahedron::nNeighbor_Elements = 4;

unsigned short CTetrahedron::VTK_Type = 10;

unsigned short CTetrahedron::maxNodesFace = 3;

CTetrahedron::CTetrahedron(unsigned long val_point_0, unsigned long val_point_1,
						   unsigned long val_point_2, unsigned long val_point_3) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = 3;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CTetrahedron::~CTetrahedron() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CTetrahedron::Change_Orientation(void) {
	unsigned long Point_0, Point_1;

	Point_0 = Nodes[0];
	Point_1 = Nodes[1];
	Nodes[0] = Point_1;
	Nodes[1] = Point_0;
  
}

unsigned short CHexahedron::Faces[6][4] = {{0,1,5,4},{1,2,6,5},{2,3,7,6},{3,0,4,7},{0,3,2,1},{4,5,6,7}};

unsigned short CHexahedron::Neighbor_Nodes[8][3] = {{1,3,4},{0,2,5},{1,3,6},{0,2,7},{0,5,7},{4,6,1},{2,5,7},{4,3,6}};

unsigned short CHexahedron::nNodesFace[6] = {4,4,4,4,4,4};

unsigned short CHexahedron::nNeighbor_Nodes[8] = {3,3,3,3,3,3,3,3};

unsigned short CHexahedron::nFaces = 6;

unsigned short CHexahedron::nNodes = 8;

unsigned short CHexahedron::nNeighbor_Elements = 6;

unsigned short CHexahedron::VTK_Type = 12;

unsigned short CHexahedron::maxNodesFace = 4;

CHexahedron::CHexahedron(unsigned long val_point_0, unsigned long val_point_1,
						 unsigned long val_point_2, unsigned long val_point_3, 
						 unsigned long val_point_4, unsigned long val_point_5, 
						 unsigned long val_point_6, unsigned long val_point_7) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate center-of-gravity coordinates ---*/
	nDim = 3;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;	Nodes[3] = val_point_3;
	Nodes[4] = val_point_4;	Nodes[5] = val_point_5;
	Nodes[6] = val_point_6;	Nodes[7] = val_point_7;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CHexahedron::~CHexahedron() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
}

void CHexahedron::Change_Orientation(void) {
	unsigned long Point_1, Point_3, Point_5, Point_7;

	Point_1 = Nodes[1];
	Point_3 = Nodes[3];
	Point_5 = Nodes[5];
	Point_7 = Nodes[7];
	Nodes[1] = Point_3;
	Nodes[3] = Point_1;
	Nodes[5] = Point_7;
	Nodes[7] = Point_5;
  
}

void CHexahedron::SetResolutionTensor(su2double** val_coord) {

  unsigned short iDim, jDim, kDim, lDim;
  unsigned short iFace;
  unsigned short* paired_faces;
  /** paired_faces is used to sort the faces into matching pairs.
      The code will look for pairs of faces that are mostly opposite, then
      sort them so that the face indices in paired_faces[0] and paired_faces[1]
      match, then paired_faces[2] and paired_faces[3] match, etc.
   */

  /*-- Allocate ResolutionTensor --*/
  ResolutionTensor = new su2double* [nDim];
  ResolutionValues = new su2double[nDim];
  ResolutionVectors = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; iDim++) {
    ResolutionTensor[iDim] = new su2double [nDim];
    ResolutionVectors[iDim] = new su2double[nDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      ResolutionTensor[iDim][jDim] = 0.0;
    }
  }

  paired_faces = new unsigned short [nFaces];

  /*-- Create cell center to face vectors --*/
  vector<vector<su2double> > center2face(nFaces, vector<su2double>(nDim));
  for (iFace = 0; iFace < nFaces; ++iFace) {
    for (iDim = 0; iDim < nDim; ++iDim) {
      center2face[iFace][iDim] = Coord_FaceElems_CG[iFace][iDim] - Coord_CG[iDim];
    }
  }

  /*-- First vector --*/

  /*-- Choose vector 1 as our first vector to pair up --*/
  paired_faces[0] = 0;
  /*-- Find vector mostly parallel to first --*/
  su2double min_dp = 1.0;
  su2double current_dp;
  for (iFace = 1; iFace < nFaces; ++iFace) {
    current_dp = inline_dot_prod(center2face[0],center2face[iFace])
        /(inline_magnitude(center2face[0])*inline_magnitude(center2face[iFace]));
    if (current_dp < min_dp) {
      min_dp = current_dp;
      paired_faces[1] = iFace;
    }
  }

  /*-- Second vector --*/

  for (iFace = 1; iFace < nFaces; ++iFace) {
    if (iFace != paired_faces[1]) {
      paired_faces[2] = iFace;
      break;
    }
  }

  min_dp = 1.0;
  for (iFace = 1; iFace < nFaces; ++iFace) {
    if (iFace == paired_faces[1]) continue;
    // current_dp = inline_dot_prod(center2face[paired_faces[2]],center2face[iFace])
    //     /(inline_magnitude(center2face[paired_faces[2]])
    //         *inline_magnitude(center2face[1]));
    current_dp = inline_dot_prod(center2face[paired_faces[2]],center2face[iFace])
        /(inline_magnitude(center2face[paired_faces[2]])*inline_magnitude(center2face[iFace]));
    if (current_dp < min_dp) {
      min_dp = current_dp;
      paired_faces[3] = iFace;
    }
  }

  /*-- Third vector --*/

  paired_faces[4] = 0;
  paired_faces[5] = 0;
  for (iFace = 1; iFace < nFaces; ++iFace) {
    if (iFace != paired_faces[1] &&
        iFace != paired_faces[2] &&
        iFace != paired_faces[3]) {
      if (paired_faces[4] == 0) {
        paired_faces[4] = iFace;
      } else {
        paired_faces[5] = iFace;
      }
    }
  }

  /*-- Use paired_faces list to build vectors --*/

  /*-- eigvecs[i][j] is the jth element of the ith eigenvector.
   *   eigvecs[i] is the ith eigenvector ---*/
  vector<vector<su2double> > eigvecs(nDim, vector<su2double>(nDim));
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      eigvecs[jDim][iDim] = Coord_FaceElems_CG[paired_faces[jDim*2]][iDim] -
          Coord_FaceElems_CG[paired_faces[jDim*2+1]][iDim];
    }
  }

  /*-- Get lengths of vectors --*/
  su2double** eigvalues = new su2double*[nDim];
  for (iDim = 0; iDim < nDim; ++iDim) {
    eigvalues[iDim] = new su2double[nDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      eigvalues[iDim][jDim] = 0.0;
    }
    eigvalues[iDim][iDim] = inline_magnitude(eigvecs[iDim]);
  }

  /*-- Gram-Schmidt Process to make the vectors orthogonal --*/
  vector<vector<su2double> > temp_eigvecs = eigvecs;
  GramSchmidt(temp_eigvecs, eigvecs);

  /*--- Change vectors to be positive ---*/
  vector<su2double> x_vector(3);
  vector<su2double> y_vector(3);
  vector<su2double> z_vector(3);
  x_vector[0] = 1.0; y_vector[1] = 1.0; z_vector[2] = 1.0;
  for (iDim = 0; iDim < nDim; iDim++) {
    su2double alignment = inline_dot_prod(eigvecs[iDim],x_vector) +
                          inline_dot_prod(eigvecs[iDim],y_vector) +
                          inline_dot_prod(eigvecs[iDim],z_vector);
    if (alignment < 0) {
      for (jDim = 0; jDim < nDim; jDim++) {
        eigvecs[iDim][jDim] *= -1;
      }
    }
  }

  /*-- Perform matrix multiplication
   * Typically, this is done as A = Q D Q^T, or
   *    Q[i][k] * D[k][m] Q[j][m]
   * where Q is a square matrix whose columns are the eigenvectors and
   * D is the diagonal matrix whose diagonal elements are the
   * corresponding eigenvalues.  But here, Q[i] is the ith eigenvector,
   * rather than Q[:][i].  So the notation is flipped to:
   *    Q[k][i] * D[k][m] Q[m][j]
   * --*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    for (jDim = 0; jDim < nDim; ++jDim) {
      for (kDim = 0; kDim < nDim; ++kDim) {
        for (lDim = 0; lDim < nDim; ++lDim) {
          ResolutionTensor[iDim][jDim] += eigvecs[kDim][iDim]
              *eigvalues[kDim][lDim]*eigvecs[lDim][jDim];
        }
      }
    }
  }

  /*--- Record calculation values in values and vectors ---*/
  for (iDim = 0; iDim < nDim; ++iDim) {
    ResolutionValues[iDim] = eigvalues[iDim][iDim];
    for (jDim = 0; jDim < nDim; ++jDim) {
      ResolutionVectors[iDim][jDim] = eigvecs[iDim][jDim];
    }
  }

  delete [] paired_faces;
  for (iDim = 0; iDim < nDim; ++iDim) {
    delete [] eigvalues[iDim];
  }
  delete [] eigvalues;


}

unsigned short CPrism::Faces[5][4] = {{3,4,1,0},{5,2,1,4},{2,5,3,0},{0,1,2,2},{5,4,3,3}};

unsigned short CPrism::Neighbor_Nodes[6][3] = {{1,2,3},{0,2,4},{1,0,5},{0,4,5},{3,5,1},{4,3,2}};

unsigned short CPrism::nNodesFace[5] = {4,4,4,3,3};

unsigned short CPrism::nNeighbor_Nodes[6] = {3,3,3,3,3,3};

unsigned short CPrism::nFaces = 5;

unsigned short CPrism::nNodes = 6;

unsigned short CPrism::nNeighbor_Elements = 5;

unsigned short CPrism::VTK_Type = 13;

unsigned short CPrism::maxNodesFace = 4;

CPrism::CPrism(unsigned long val_point_0, unsigned long val_point_1, 
			   unsigned long val_point_2, unsigned long val_point_3, 
			   unsigned long val_point_4, unsigned long val_point_5) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/
	nDim = 3;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++) Coord_CG[iDim] = 0.0;
	
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/
	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	Nodes[4] = val_point_4;
	Nodes[5] = val_point_5;
	
	/*--- Allocate and define neighbor elements to a element ---*/
	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CPrism::~CPrism() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CPrism::Change_Orientation(void) {
	unsigned long Point_0, Point_1, Point_3, Point_4;

	Point_0 = Nodes[0];
	Point_1 = Nodes[1];
	Point_3 = Nodes[3];
	Point_4 = Nodes[4];
	Nodes[0] = Point_1;
	Nodes[1] = Point_0;
	Nodes[3] = Point_4;
	Nodes[4] = Point_3;
  
}

unsigned short CPyramid::Faces[5][4] = {{0,3,2,1},{4,3,0,0},{4,0,1,1},{2,4,1,1},{3,4,2,2}};

unsigned short CPyramid::Neighbor_Nodes[5][4] = {{1,3,4,4},{0,2,4,4},{1,3,4,4},{2,0,4,4},{0,1,2,3}};

unsigned short CPyramid::nNodesFace[5] = {4,3,3,3,3};

unsigned short CPyramid::nNeighbor_Nodes[5] = {3,3,3,3,4};

unsigned short CPyramid::nFaces = 5;

unsigned short CPyramid::nNodes = 5;

unsigned short CPyramid::nNeighbor_Elements = 5;

unsigned short CPyramid::VTK_Type = 14;

unsigned short CPyramid::maxNodesFace = 4;

CPyramid::CPyramid(unsigned long val_point_0, unsigned long val_point_1,
				   unsigned long val_point_2, unsigned long val_point_3, 
				   unsigned long val_point_4) : CPrimalGrid() {
	unsigned short iDim, iFace, iNeighbor_Elements;

	/*--- Allocate CG coordinates ---*/

	nDim = 3;
	Coord_CG = new su2double[nDim];
	for (iDim = 0; iDim < nDim; iDim++)
		Coord_CG[iDim] = 0.0;
	Coord_FaceElems_CG = new su2double* [nFaces];
	for (iFace = 0; iFace < nFaces; iFace++) {
		Coord_FaceElems_CG[iFace] = new su2double [nDim];
		for (iDim = 0; iDim < nDim; iDim++)
			Coord_FaceElems_CG[iFace][iDim] = 0.0;
	}
	
	/*--- Allocate and define face structure of the element ---*/

	Nodes = new unsigned long[nNodes];
	Nodes[0] = val_point_0;
	Nodes[1] = val_point_1;
	Nodes[2] = val_point_2;
	Nodes[3] = val_point_3;
	Nodes[4] = val_point_4;
	
	/*--- Allocate and define neighbor elements to a element ---*/

	nNeighbor_Elements = nFaces;
	Neighbor_Elements = new long[nNeighbor_Elements];
	for (iNeighbor_Elements = 0; iNeighbor_Elements<nNeighbor_Elements; iNeighbor_Elements++) {
		Neighbor_Elements[iNeighbor_Elements]=-1;
	}
  
}

CPyramid::~CPyramid() {
  unsigned short iFaces;
  
  for (iFaces = 0; iFaces < nFaces; iFaces++)
    if (Coord_FaceElems_CG[iFaces] != NULL) delete[] Coord_FaceElems_CG[iFaces];
  if (Coord_FaceElems_CG != NULL) delete[] Coord_FaceElems_CG;
  
}

void CPyramid::Change_Orientation(void) {
	unsigned long Point_1, Point_3;

	Point_1 = Nodes[1];
	Point_3 = Nodes[3];
	Nodes[1] = Point_3;
	Nodes[3] = Point_1;

}

CPrimalGridFEM::CPrimalGridFEM(unsigned long  val_elemGlobalID, unsigned short val_VTK_Type,
                               unsigned short val_nPolyGrid,    unsigned short val_nPolySol,
                               unsigned short val_nDOFsGrid,    unsigned short val_nDOFsSol,
                               unsigned long  val_offDOfsSol,   istringstream  &elem_line)
{
  /*--- Store the integer data in the member variables of this object. ---*/
  VTK_Type = val_VTK_Type;
  nDim = (VTK_Type == TRIANGLE || VTK_Type == QUADRILATERAL) ? 2 : 3;

  nPolyGrid = val_nPolyGrid;
  nPolySol  = val_nPolySol;
  nDOFsGrid = val_nDOFsGrid;
  nDOFsSol  = val_nDOFsSol;

  elemIDGlobal        = val_elemGlobalID;
  offsetDOFsSolGlobal = val_offDOfsSol;

  /*--- Allocate the memory for the global nodes of the element to define
        the geometry and read them from elem_line.                        ---*/
  Nodes = new unsigned long[nDOFsGrid];
  for(unsigned short i=0; i<nDOFsGrid; i++)
    elem_line >> Nodes[i];

  /*--- If a linear element is used, the node numbering for non-simplices
        must be adapted. The reason is that compatability with the original
        SU2 format is maintained for linear elements, but for the FEM solver
        the nodes of the elements are stored row-wise.                       ---*/
  if(nPolyGrid == 1){
    switch( VTK_Type ) {

      case QUADRILATERAL:
        swap(Nodes[2], Nodes[3]);
        break;

      case HEXAHEDRON:
        swap(Nodes[2], Nodes[3]);
        swap(Nodes[6], Nodes[7]);
        break;

      case PYRAMID:
        swap(Nodes[2], Nodes[3]);
        break;
    }
  }
}

CPrimalGridFEM::CPrimalGridFEM(unsigned long  val_elemGlobalID, unsigned short val_VTK_Type,
                               unsigned short val_nPolyGrid,    unsigned short val_nPolySol,
                               unsigned short val_nDOFsGrid,    unsigned short val_nDOFsSol,
                               unsigned long  val_offDOfsSol,   const unsigned long *connGrid)
{
  /*--- Store the integer data in the member variables of this object. ---*/
  VTK_Type = val_VTK_Type;
  nDim = (VTK_Type == TRIANGLE || VTK_Type == QUADRILATERAL) ? 2 : 3;

  nPolyGrid = val_nPolyGrid;
  nPolySol  = val_nPolySol;
  nDOFsGrid = val_nDOFsGrid;
  nDOFsSol  = val_nDOFsSol;

  elemIDGlobal        = val_elemGlobalID;
  offsetDOFsSolGlobal = val_offDOfsSol;

  /*--- Allocate the memory for the global nodes of the element to define
        the geometry and copy the data from connGrid. ---*/
  Nodes = new unsigned long[nDOFsGrid];
  for(unsigned short i=0; i<nDOFsGrid; i++)
    Nodes[i] = connGrid[i];
}

CPrimalGridFEM::~CPrimalGridFEM(){}

void CPrimalGridFEM::GetLocalCornerPointsAllFaces(unsigned short elementType,
                                                  unsigned short nPoly,
                                                  unsigned short nDOFs,
                                                  unsigned short &numFaces,
                                                  unsigned short nPointsPerFace[],
                                                  unsigned long  faceConn[6][4]) {

  /*--- Determine the element type and set the face data accordingly.
        The faceConn values are local to the element.                 ---*/

  unsigned short nn2, nn3, nn4;
  switch( elementType ) {
    case TRIANGLE:
      numFaces = 3;
      nPointsPerFace[0] = 2; faceConn[0][0] = 0;        faceConn[0][1] = nPoly;
      nPointsPerFace[1] = 2; faceConn[1][0] = nPoly;    faceConn[1][1] = nDOFs -1;
      nPointsPerFace[2] = 2; faceConn[2][0] = nDOFs -1; faceConn[2][1] = 0;
      break;

    case QUADRILATERAL:
      numFaces = 4; nn2 = nPoly*(nPoly+1);
      nPointsPerFace[0] = 2; faceConn[0][0] = 0;        faceConn[0][1] = nPoly;
      nPointsPerFace[1] = 2; faceConn[1][0] = nPoly;    faceConn[1][1] = nDOFs -1;
      nPointsPerFace[2] = 2; faceConn[2][0] = nDOFs -1; faceConn[2][1] = nn2;
      nPointsPerFace[3] = 2; faceConn[3][0] = nn2;      faceConn[3][1] = 0;
      break;

    case TETRAHEDRON:
      numFaces = 4; nn2 = (nPoly+1)*(nPoly+2)/2 -1; nn3 = nDOFs -1;
      nPointsPerFace[0] = 3; faceConn[0][0] = 0;     faceConn[0][1] = nPoly; faceConn[0][2] = nn2;
      nPointsPerFace[1] = 3; faceConn[1][0] = 0;     faceConn[1][1] = nn3;   faceConn[1][2] = nPoly;
      nPointsPerFace[2] = 3; faceConn[2][0] = 0;     faceConn[2][1] = nn2;   faceConn[2][2] = nn3;
      nPointsPerFace[3] = 3; faceConn[3][0] = nPoly; faceConn[3][1] = nn3;   faceConn[3][2] = nn2;
      break;

    case PYRAMID:
      numFaces = 5; nn2 = (nPoly+1)*(nPoly+1) -1; nn3 = nn2 - nPoly;
      nPointsPerFace[0] = 4; faceConn[0][0] = 0;     faceConn[0][1] = nPoly;    faceConn[0][2] = nn2; faceConn[0][3] = nn3;
      nPointsPerFace[1] = 3; faceConn[1][0] = 0;     faceConn[1][1] = nDOFs -1; faceConn[1][2] = nPoly;
      nPointsPerFace[2] = 3; faceConn[2][0] = nn3;   faceConn[2][1] = nn2;      faceConn[2][2] = nDOFs -1;
      nPointsPerFace[3] = 3; faceConn[3][0] = 0;     faceConn[3][1] = nn3;      faceConn[3][2] = nDOFs -1;
      nPointsPerFace[4] = 3; faceConn[4][0] = nPoly; faceConn[4][1] = nDOFs -1; faceConn[4][2] = nn2;
      break;

    case PRISM:
      numFaces = 5; nn2 = (nPoly+1)*(nPoly+2)/2; nn3 = nPoly*nn2; --nn2;
      nPointsPerFace[0] = 3; faceConn[0][0] = 0;     faceConn[0][1] = nPoly;     faceConn[0][2] = nn2;
      nPointsPerFace[1] = 3; faceConn[1][0] = nn3;   faceConn[1][1] = nn2+nn3;   faceConn[1][2] = nPoly+nn3;
      nPointsPerFace[2] = 4; faceConn[2][0] = 0;     faceConn[2][1] = nn3;       faceConn[2][2] = nPoly+nn3; faceConn[2][3] = nPoly;
      nPointsPerFace[3] = 4; faceConn[3][0] = 0;     faceConn[3][1] = nn2;       faceConn[3][2] = nn2+nn3;   faceConn[3][3] = nn3;
      nPointsPerFace[4] = 4; faceConn[4][0] = nPoly; faceConn[4][1] = nPoly+nn3; faceConn[4][2] = nn2+nn3;   faceConn[4][3] = nn2;
      break;

    case HEXAHEDRON:
      numFaces = 6; nn2 = (nPoly+1)*(nPoly+1); nn4 = nPoly*nn2; --nn2; nn3 = nn2 - nPoly;
      nPointsPerFace[0] = 4; faceConn[0][0] = 0;     faceConn[0][1] = nPoly;     faceConn[0][2] = nn2;       faceConn[0][3] = nn3;
      nPointsPerFace[1] = 4; faceConn[1][0] = nn4;   faceConn[1][1] = nn3+nn4;   faceConn[1][2] = nn2+nn4;   faceConn[1][3] = nPoly+nn4;
      nPointsPerFace[2] = 4; faceConn[2][0] = 0;     faceConn[2][1] = nn4;       faceConn[2][2] = nPoly+nn4; faceConn[2][3] = nPoly;
      nPointsPerFace[3] = 4; faceConn[3][0] = nn3;   faceConn[3][1] = nn2;       faceConn[3][2] = nn2+nn4;   faceConn[3][3] = nn3+nn4;
      nPointsPerFace[4] = 4; faceConn[4][0] = 0;     faceConn[4][1] = nn3;       faceConn[4][2] = nn3+nn4;   faceConn[4][3] = nn4;
      nPointsPerFace[5] = 4; faceConn[5][0] = nPoly; faceConn[5][1] = nPoly+nn4; faceConn[5][2] = nn2+nn4;   faceConn[5][3] = nn2;
      break;
  }
}

void CPrimalGridFEM::GetCornerPointsAllFaces(unsigned short &numFaces,
                                             unsigned short nPointsPerFace[],
                                             unsigned long  faceConn[6][4]) {

  /*--- Get the corner points of the faces local to the element. ---*/

  GetLocalCornerPointsAllFaces(VTK_Type, nPolyGrid, nDOFsGrid,
                               numFaces, nPointsPerFace, faceConn);

  /*--- Convert the local values of faceConn to global values. ---*/

  for(unsigned short i=0; i<numFaces; ++i) {
    for(unsigned short j=0; j<nPointsPerFace[i]; ++j) {
      unsigned long nn = faceConn[i][j];
      faceConn[i][j] = Nodes[nn];
    }
  }

  /*--- Store numFaces in nFaces for later purposes. ---*/

  nFaces = numFaces;
}

CPrimalGridBoundFEM::CPrimalGridBoundFEM(unsigned long         val_elemGlobalID,
                                         unsigned long         val_domainElementID,
                                         unsigned short        val_VTK_Type,
                                         unsigned short        val_nPolyGrid,
                                         unsigned short        val_nDOFsGrid,
                                         vector<unsigned long> &val_nodes)
{
  /*--- Store the integer data in the member variables of this object. ---*/

  VTK_Type = val_VTK_Type;
  nDim = (VTK_Type == LINE) ? 1 : 2;

  nPolyGrid = val_nPolyGrid;
  nDOFsGrid = val_nDOFsGrid;

  boundElemIDGlobal = val_elemGlobalID;
  DomainElement     = val_domainElementID;

  /*--- Allocate the memory for the global nodes of the element to define
        the geometry and copy them from val_nodes.                        ---*/

  Nodes = new unsigned long[nDOFsGrid];
  for(unsigned short i=0; i<nDOFsGrid; i++)
    Nodes[i] = val_nodes[i];

  /*--- For a linear quadrilateral the two last node numbers must be swapped,
        such that the element numbering is consistent with the FEM solver.    ---*/

  if(nPolyGrid == 1 && VTK_Type == QUADRILATERAL) swap(Nodes[2], Nodes[3]);
}

void CPrimalGridBoundFEM::GetLocalCornerPointsFace(unsigned short elementType,
                                                   unsigned short nPoly,
                                                   unsigned short nDOFs,
                                                   unsigned short &nPointsPerFace,
                                                   unsigned long  faceConn[]) {

  /*--- Determine the element type and set the face data accordingly.
        The faceConn values are local to the element.                 ---*/

  switch( elementType ) {
    case LINE:
      nPointsPerFace = 2;
      faceConn[0] = 0; faceConn[1] = nPoly;
      break;

    case TRIANGLE:
      nPointsPerFace = 3;
      faceConn[0] = 0; faceConn[1] = nPoly; faceConn[2] = nDOFs -1;
      break;

    case QUADRILATERAL:
      unsigned short nn2 = nPoly*(nPoly+1);
      nPointsPerFace = 4;
      faceConn[0] = 0; faceConn[1] = nPoly; faceConn[2] = nDOFs -1; faceConn[3] = nn2;
      break;
  }
}

void CPrimalGridBoundFEM::GetCornerPointsAllFaces(unsigned short &nFaces,
                                                  unsigned short nPointsPerFace[],
                                                  unsigned long  faceConn[6][4]) {

  /*--- For compatibility reasons with the base class, nFaces is also an
        argument to this function. However, it is always 1, because the
        boundary element is the face.                                   ---*/
  nFaces = 1;

  /*--- Get the corner points of the face local to the element. ---*/

  unsigned long thisFaceConn[4];
  GetLocalCornerPointsFace(VTK_Type, nPolyGrid, nDOFsGrid,
                           nPointsPerFace[0], thisFaceConn);


  /*--- Convert the local values of thisFaceConn to global values. ---*/

  for(unsigned short j=0; j<nPointsPerFace[0]; ++j)
    faceConn[0][j] = Nodes[thisFaceConn[j]];
}

void CPrimalGridBoundFEM::RemoveMultipleDonorsWallFunctions(void) {

  /* Sort donorElementsWallFunctions in increasing order and remove the
     the double entities. */
  sort(donorElementsWallFunctions.begin(), donorElementsWallFunctions.end());
  vector<unsigned long>::iterator lastEntry;
  lastEntry = unique(donorElementsWallFunctions.begin(),
                     donorElementsWallFunctions.end());
  donorElementsWallFunctions.erase(lastEntry, donorElementsWallFunctions.end());
}

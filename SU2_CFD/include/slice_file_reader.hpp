/*!
 * \file slice_file_reader.hpp
 * \brief Header used to input a slice as a restart file for 3D problems.
 * \author C. Pederson
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

#pragma once

#include "../../Common/include/mpi_structure.hpp"

#include <cassert>
#include <cmath>
#include <string>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <iostream>
#include <set>
#include <stdlib.h>
#include <stdio.h>

#include "variable_structure.hpp"
#include "../../Common/include/geometry_structure.hpp"
#include "../../Common/include/config_structure.hpp"

using namespace std;

class CAbstractFileReader {
 protected:

  int rank, size;
  int* Restart_Vars;
  passivedouble* Restart_Data;

  /**
   * This must be overriden in all base clases.
   */
  virtual void FindClosestPoint(const su2double* coord_i,
                                unsigned long& jPoint_closest,
                                su2double& min_distance) const = 0;
  /**
   * This must be overriden in all base clases.
   */
  virtual void TransformSolution(unsigned long start_index,
                                 unsigned short val_nVar,
                                 const su2double* coord_i,
                                 su2double* solution) const = 0;

 public:

  CAbstractFileReader();
  ~CAbstractFileReader();

  void Read_SliceFile_ASCII(CConfig *config, string val_filename);

  void LoadSolutionFromSlice(const string& restart_filename,
                             CConfig* config,
                             CGeometry* geometry,
                             unsigned short nVar,
                             unsigned short offset,
                             CVariable** node) const;

  const passivedouble* GetRestart_Data() const {
    assert(Restart_Data);
    return Restart_Data;
  }
};

class CFileReader_Cartesian : public CAbstractFileReader {
 protected:
  void FindClosestPoint(const su2double* coord_i,
                        unsigned long& jPoint_closest,
                        su2double& min_distance) const override;
  void TransformSolution(unsigned long start_index,
                         unsigned short val_nVar,
                         const su2double* coord_i,
                         su2double* solution) const override;

 public:
  CFileReader_Cartesian();
};

class CFileReader_Cylindrical : public CAbstractFileReader {
 protected:
  void FindClosestPoint(const su2double* coord_i,
                        unsigned long& jPoint_closest,
                        su2double& min_distance) const override;
  void TransformSolution(unsigned long start_index,
                         unsigned short val_nVar,
                         const su2double* coord_i,
                         su2double* solution) const override;

 public:
  CFileReader_Cylindrical();
};

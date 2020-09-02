/*!
 * \file slice_file_reader.cpp
 * \brief Main subrotuines for reading in a slice as a restart file
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

#include "../include/slice_file_reader.hpp"

CAbstractFileReader::CAbstractFileReader() {
  rank = SU2_MPI::GetRank();
  size = SU2_MPI::GetSize();
};

CAbstractFileReader::~CAbstractFileReader() {
  if (Restart_Vars != NULL) {delete [] Restart_Vars; Restart_Vars = NULL;}
  if (Restart_Data != NULL) {delete [] Restart_Data; Restart_Data = NULL;}
};

void CAbstractFileReader::Read_SliceFile_ASCII(CConfig *config,
                                               string val_filename) {

  /*--- First, check that this is not a binary restart file. ---*/

  char fname[100];
  strcpy(fname, val_filename.c_str());
  int magic_number;

#ifndef HAVE_MPI

  /*--- Serial binary input. ---*/

  FILE *fhw;
  fhw = fopen(fname,"rb");
  size_t ret;

  /*--- Error check for opening the file. ---*/

  if (!fhw) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + fname, CURRENT_FUNCTION);
  }

  /*--- Attempt to read the first int, which should be our magic number. ---*/

  ret = fread(&magic_number, sizeof(int), 1, fhw);
  if (ret != 1) {
    SU2_MPI::Error("Error reading restart file.", CURRENT_FUNCTION);
  }

  /*--- Check that this is an SU2 binary file. SU2 binary files
    have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
        string("SU2 reads/writes binary restart files by default.\n") +
        string("Note that backward compatibility for ASCII restart files is\n") +
        string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
  }

  fclose(fhw);

#else

  /*--- Parallel binary input using MPI I/O. ---*/

  MPI_File fhw;
  int ierr;

  /*--- All ranks open the file using MPI. ---*/

  ierr = MPI_File_open(MPI_COMM_WORLD, fname, MPI_MODE_RDONLY, MPI_INFO_NULL, &fhw);

  /*--- Error check opening the file. ---*/

  if (ierr) {
    SU2_MPI::Error(string("Unable to open SU2 restart file ") + string(fname), CURRENT_FUNCTION);
  }

  /*--- Have the master attempt to read the magic number. ---*/

  if (rank == MASTER_NODE)
    MPI_File_read(fhw, &magic_number, 1, MPI_INT, MPI_STATUS_IGNORE);

  /*--- Broadcast the number of variables to all procs and store clearly. ---*/

  SU2_MPI::Bcast(&magic_number, 1, MPI_INT, MASTER_NODE, MPI_COMM_WORLD);

  /*--- Check that this is an SU2 binary file. SU2 binary files
    have the hex representation of "SU2" as the first int in the file. ---*/

  if (magic_number == 535532) {
    SU2_MPI::Error(string("File ") + string(fname) + string(" is a binary SU2 restart file, expected ASCII.\n") +
        string("SU2 reads/writes binary restart files by default.\n") +
        string("Note that backward compatibility for ASCII restart files is\n") +
        string("possible with the WRT_BINARY_RESTART / READ_BINARY_RESTART options."), CURRENT_FUNCTION);
  }

  MPI_File_close(&fhw);

#endif

  /*--- Open the restart file ---*/

  ifstream restart_file;
  restart_file.open(val_filename.data(), ios::in);

  /*--- In case there is no restart file ---*/

  if (restart_file.fail()) {
    SU2_MPI::Error(string("SU2 ASCII solution file  ") + string(fname) + string(" not found."), CURRENT_FUNCTION);
  }

  /*--- Count the number of lines in the restart file ---*/

  string text_line;
  unsigned long nPoints_slice = 0;
  while (std::getline(restart_file, text_line)) ++nPoints_slice;
  if (nPoints_slice > 1) {
    nPoints_slice -= 1;
  } else {
    SU2_MPI::Error("Did not find any valid lines in the slice file.",
                   CURRENT_FUNCTION);
  }
  restart_file.clear();
  restart_file.seekg(0);

  /*--- Identify the number of fields (and names) in the restart file ---*/

  getline (restart_file, text_line);
  stringstream ss(text_line);
  string Tag;
  config->fields.clear();
  while (ss >> Tag) {
    config->fields.push_back(Tag);
    if (ss.peek() == ',') ss.ignore();
  }

  /*--- Set the number of variables. Don't subtract 1; the slice
   * files don't have any iPoint field. ---*/

  Restart_Vars = new int[5];
  Restart_Vars[1] = (int)config->fields.size();
  Restart_Vars[2] = (int)nPoints_slice;

  /*--- Allocate memory for the restart data. ---*/

  Restart_Data = new passivedouble[Restart_Vars[1]*nPoints_slice];

  /*--- Read all lines in the restart file and extract data. ---*/

  int counter = 0;
  for (unsigned long iPoint = 0; iPoint < nPoints_slice; iPoint++ ) {

    getline (restart_file, text_line);

    istringstream point_line(text_line);

    /*--- We don't worry about whether or not this point lives on the
     * processor. All MPI tasks use the same interpolation data. ---*/

    /*--- Store the solution (starting with node coordinates) --*/

    for (unsigned short iVar = 0; iVar < Restart_Vars[1]; iVar++)
      point_line >> Restart_Data[counter*Restart_Vars[1] + iVar];

    /*--- Increment our local point counter. ---*/

    counter++;

  }
  restart_file.close();
}

void CAbstractFileReader::LoadSolutionFromSlice(const string& restart_filename,
                                                CConfig* config,
                                                CGeometry* geometry,
                                                unsigned short nVar,
                                                unsigned short offset,
                                                CVariable** node) const {
  if (geometry->GetnDim() < 3) {
    SU2_MPI::Error("Slice remapping cannot be used with a 2D geometry.", CURRENT_FUNCTION);
  }

  /*--- The number of variables per row in the restart data ---*/

  const unsigned short nVars_Restart = Restart_Vars[1];

  /*--- Load data from the restart into correct containers. ---*/

  unsigned long iPoint_Global_Local = 0;
  su2double local_max_distance = 0;
  su2double* Solution = new su2double[nVar];
  for (unsigned long iPoint_Global = 0; iPoint_Global < geometry->GetGlobal_nPointDomain(); iPoint_Global++ ) {

    /*--- Retrieve local index. If this node from the restart file lives
     on the current processor, we will load and instantiate the vars. ---*/

    const long iPoint_Local = geometry->GetGlobal_to_Local_Point(iPoint_Global);

    if (iPoint_Local > -1) {

      const su2double* coord_i = geometry->node[iPoint_Local]->GetCoord();

      unsigned long jPoint_closest;
      su2double min_distance;
      FindClosestPoint(coord_i, jPoint_closest, min_distance);
      if (min_distance > local_max_distance) {
        local_max_distance = min_distance;
      }

      /*--- The slice file should only be 2D, and the flow vars
       * should be 3D ---*/
      const unsigned short skipVars = 2 + offset;

      /*--- We need to store this point's data, so jump to the correct
       offset in the buffer of data from the restart file and load it. ---*/

      const long index = jPoint_closest*nVars_Restart + skipVars;
      for (unsigned short iVar = 0; iVar < nVar; iVar++) {
        Solution[iVar] = Restart_Data[index+iVar];
      }
      node[iPoint_Local]->SetSolution(Solution);
      iPoint_Global_Local++;
    }
  }

  delete [] Solution;

  /*--- Detect a wrong solution file ---*/

  unsigned short rbuf_NotMatching = 0, sbuf_NotMatching = 0;
  if (iPoint_Global_Local < geometry->GetnPointDomain()) { sbuf_NotMatching = 1; }

#ifndef HAVE_MPI
  rbuf_NotMatching = sbuf_NotMatching;
#else
  SU2_MPI::Allreduce(&sbuf_NotMatching, &rbuf_NotMatching, 1, MPI_UNSIGNED_SHORT, MPI_SUM, MPI_COMM_WORLD);
#endif
  if (rbuf_NotMatching != 0) {
      SU2_MPI::Error(string("The solution file ") + restart_filename + string(" doesn't match with the mesh file!\n") +
                     string("It could be empty lines at the end of the file."), CURRENT_FUNCTION);
  }

  su2double global_max_distance;
  SU2_MPI::Allreduce(&local_max_distance, &global_max_distance, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);
  if (rank == MASTER_NODE) {
    cout << "Maximum distance in slice remapping: ";
    cout << global_max_distance << endl;
  }
}

CFileReader_Cartesian::CFileReader_Cartesian() : CAbstractFileReader() {}

void CFileReader_Cartesian::FindClosestPoint(const su2double* coord_i,
                      unsigned long& jPoint_closest,
                      su2double& min_distance) const {
  const unsigned short nVars_Restart = Restart_Vars[1];
  min_distance = 1E16;
  jPoint_closest = 0;
  for (unsigned long jPoint = 0; jPoint < Restart_Vars[2]; jPoint++) {
    const long index = jPoint * nVars_Restart;
    const su2double x = Restart_Data[index];
    const su2double y = Restart_Data[index + 1];
    const su2double dist = sqrt(pow(coord_i[0] - x, 2) +
                                pow(coord_i[1] - y, 2));
    if (dist < min_distance) {
      jPoint_closest = jPoint;
      min_distance = dist;
    }
  }
}

CFileReader_Cylindrical::CFileReader_Cylindrical() : CAbstractFileReader() {}

void CFileReader_Cylindrical::FindClosestPoint(const su2double* coord_i,
                      unsigned long& jPoint_closest,
                      su2double& min_distance) const {

  // The axis is hardcoded to be oriented in the x-direction
  // through the point (0, 0, 0)

  const su2double axis_y = 0;
  const su2double axis_z = 0;

  /*--- "i" is the point point on the geometry (3D grid) ---*/

  const su2double x_i = coord_i[0];
  const su2double r_i = sqrt(pow(coord_i[1] - axis_y, 2) +
                             pow(coord_i[2] - axis_z, 2));

  const unsigned short nVars_Restart = Restart_Vars[1];
  min_distance = 1E16;
  jPoint_closest = 0;
  for (unsigned long jPoint = 0; jPoint < Restart_Vars[2]; jPoint++) {

    /*--- "j" is the point in the slice ---*/

    const long index = jPoint * nVars_Restart;
    const su2double x_j = Restart_Data[index];
    const su2double r_j = Restart_Data[index + 1];
    const su2double dist = sqrt(pow(x_i - x_j, 2) +
                                pow(r_i - r_j, 2));
    if (dist < min_distance) {
      jPoint_closest = jPoint;
      min_distance = dist;
    }
  }
}

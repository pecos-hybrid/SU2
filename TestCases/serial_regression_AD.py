#!/usr/bin/env python

## \file serial_regression.py
#  \brief Python script for automated regression testing of SU2 examples
#  \author A. Aranake, A. Campos, T. Economon, T. Lukaczyk, S. Padron
#  \version 6.2.0 "Falcon"
#
# The current SU2 release has been coordinated by the
# SU2 International Developers Society <www.su2devsociety.org>
# with selected contributions from the open-source community.
#
# The main research teams contributing to the current release are:
#  - Prof. Juan J. Alonso's group at Stanford University.
#  - Prof. Piero Colonna's group at Delft University of Technology.
#  - Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#  - Prof. Alberto Guardone's group at Polytechnic University of Milan.
#  - Prof. Rafael Palacios' group at Imperial College London.
#  - Prof. Vincent Terrapon's group at the University of Liege.
#  - Prof. Edwin van der Weide's group at the University of Twente.
#  - Lab. of New Concepts in Aeronautics at Tech. Institute of Aeronautics.
#
# Copyright 2012-2019, Francisco D. Palacios, Thomas D. Economon,
#                      Tim Albring, and the SU2 contributors.
#
# SU2 is free software; you can redistribute it and/or
# modify it under the terms of the GNU Lesser General Public
# License as published by the Free Software Foundation; either
# version 2.1 of the License, or (at your option) any later version.
#
# SU2 is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the GNU
# Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public
# License along with SU2. If not, see <http://www.gnu.org/licenses/>.

# make print(*args) function available in PY2.6+, does'nt work on PY < 2.6
from __future__ import print_function

import sys
from TestCase import TestCase

def main():
    '''This program runs SU2 and ensures that the output matches specified values. 
       This will be used to do checks when code is pushed to github 
       to make sure nothing is broken. '''

    test_list = []

    #######################################################
    ### Disc. adj. compressible RANS                    ###
    #######################################################
    
    # Adjoint turbulent NACA0012 SA
    discadj_rans_naca0012_sa           = TestCase('discadj_rans_naca0012_sa')
    discadj_rans_naca0012_sa.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sa.cfg_file  = "turb_NACA0012_sa.cfg"
    discadj_rans_naca0012_sa.test_iter = 10
    discadj_rans_naca0012_sa.test_vals = [-1.751962, 0.485751, 0.182121, -0.000018] #last 4 columns
    discadj_rans_naca0012_sa.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sa.timeout   = 1600
    discadj_rans_naca0012_sa.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sa)

    # Adjoint turbulent NACA0012 SST
    discadj_rans_naca0012_sst           = TestCase('discadj_rans_naca0012_sst')
    discadj_rans_naca0012_sst.cfg_dir   = "disc_adj_rans/naca0012"
    discadj_rans_naca0012_sst.cfg_file  = "turb_NACA0012_sst.cfg"
    discadj_rans_naca0012_sst.test_iter = 10
    discadj_rans_naca0012_sst.test_vals = [-1.654903, -0.491485, 0.109157, 0.000011] #last 4 columns
    discadj_rans_naca0012_sst.su2_exec  = "SU2_CFD_AD"
    discadj_rans_naca0012_sst.timeout   = 1600
    discadj_rans_naca0012_sst.tol       = 0.00001
    test_list.append(discadj_rans_naca0012_sst)

    #######################################################
    ### Unsteady Disc. adj. compressible RANS           ###
    #######################################################
   
    # Turbulent Cylinder
    discadj_cylinder           = TestCase('unsteady_cylinder')
    discadj_cylinder.cfg_dir   = "disc_adj_rans/cylinder"
    discadj_cylinder.cfg_file  = "cylinder.cfg" 
    discadj_cylinder.test_iter = 9
    discadj_cylinder.test_vals = [3.746904, -1.544886, -0.008345, 0.000014] #last 4 columns
    discadj_cylinder.su2_exec  = "SU2_CFD_AD"
    discadj_cylinder.timeout   = 1600
    discadj_cylinder.tol       = 0.00001
    discadj_cylinder.unsteady  = True
    test_list.append(discadj_cylinder)
    
    ###################################
    ### Coupled FSI Adjoint         ###
    ###################################
   
    # Structural model
    discadj_fsi           = TestCase('discadj_fsi')
    discadj_fsi.cfg_dir   = "disc_adj_fsi"
    discadj_fsi.cfg_file  = "configAD_fsi.cfg" 
    discadj_fsi.test_iter = 3000
    discadj_fsi.test_vals = [0.958848,-0.157183,0.658415,1.302076] #last 4 columns
    discadj_fsi.su2_exec  = "SU2_CFD_AD"
    discadj_fsi.timeout   = 1600
    discadj_fsi.tol       = 0.00001
    test_list.append(discadj_fsi)      

    ######################################
    ### RUN TESTS                      ###
    ######################################  

    pass_list = [ test.run_test() for test in test_list ]
    
    ######################################
    ### RUN PYTHON TESTS               ###
    ######################################
    
    # test discrete_adjoint.py
    discadj_euler_py = TestCase('discadj_euler_py')
    discadj_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    discadj_euler_py.cfg_file  = "inv_NACA0012.cfg"
    discadj_euler_py.test_iter = 10
    discadj_euler_py.su2_exec  = "discrete_adjoint.py"
    discadj_euler_py.timeout   = 1600
    discadj_euler_py.reference_file = "of_grad_cd_disc.dat.ref"
    discadj_euler_py.test_file = "of_grad_cd.dat"
    pass_list.append(discadj_euler_py.run_filediff())
    test_list.append(discadj_euler_py)
    
    # test direct_differentiation.py
    directdiff_euler_py = TestCase('directdiff_euler_py')
    directdiff_euler_py.cfg_dir = "cont_adj_euler/naca0012"
    directdiff_euler_py.cfg_file  = "inv_NACA0012_FD.cfg"
    directdiff_euler_py.test_iter = 10
    directdiff_euler_py.su2_exec  = "direct_differentiation.py"
    directdiff_euler_py.timeout   = 1600
    directdiff_euler_py.reference_file = "of_grad_directdiff.dat.ref"
    directdiff_euler_py.test_file = "DIRECTDIFF/of_grad_directdiff.dat"
    pass_list.append(directdiff_euler_py.run_filediff())
    test_list.append(directdiff_euler_py)

    # Tests summary
    print('==================================================================')
    print('Summary of the serial tests')
    print('python version:', sys.version)
    for i, test in enumerate(test_list):
        if (pass_list[i]):
            print('  passed - %s'%test.tag)
        else:
            print('* FAILED - %s'%test.tag)
    
    if all(pass_list):
        sys.exit(0)
    else:
        sys.exit(1)
    # done

if __name__ == '__main__':
    main()

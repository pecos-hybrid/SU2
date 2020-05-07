#!/usr/bin/env python

## \file generate_ascii_restart.py
#  \brief Python script for generating an ASCII restart file
#  \author C. Pederson
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

# imports
import numpy as np
import argparse
import sys, os
sys.path.append(os.environ['SU2_RUN'])
import SU2

def main(options):

    # load config, start state
    config = SU2.io.Config(options.filename)
    state  = SU2.io.State()

    # prepare config
    config.NUMBER_PART = options.partitions
    config.RESTART_SOL = 'NO'
    config.UNSTEADY_SIMULATION= 'NO'
    config.EXT_ITER = 1
    config.CFL_NUMBER = 1E-3
    config.TIME_DISCRE_FLOW = 'EULER_EXPLICIT'
    config.WRT_BINARY_RESTART = 'NO'
    config.LOW_MEMORY_OUTPUT = 'YES'
    config.WRT_RESIDUALS = 'NO'
    config.WRT_LIMITERS= 'NO'
    config.WRT_RESOLUTION_TENSOR= 'NO'
    config.RESTART_FLOW_FILENAME = options.output

    info = SU2.run.CFD(config)

if __name__ == "__main__":
    # Command Line Options
    short_description="Quickly creates an ASCII restart file."
    parser = argparse.ArgumentParser(description=short_description)
    parser.add_argument("-f", "--file", dest="filename",
                        help="read config from FILE", metavar="FILE")
    parser.add_argument("-n", "--partitions", dest="partitions", default=1,
                        help="number of PARTITIONS", metavar="PARTITIONS")
    parser.add_argument("-o", "--output", dest="output",
                        default="restart_flow_ascii.dat",
                        help="name of the ASCII output file", metavar="RESTART")

    args = parser.parse_args()
    args.partitions = int( args.partitions )
    main(args)

################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD tests
# \author T. A. Oliver
# \version 5.0.0 "Raven"
#
# SU2 Lead Developers: Dr. Francisco Palacios (Francisco.D.Palacios@boeing.com).
#                      Dr. Thomas D. Economon (economon@stanford.edu).
#
# SU2 Developers: Prof. Juan J. Alonso's group at Stanford University.
#                 Prof. Piero Colonna's group at Delft University of Technology.
#                 Prof. Nicolas R. Gauger's group at Kaiserslautern University of Technology.
#                 Prof. Alberto Guardone's group at Polytechnic University of Milan.
#                 Prof. Rafael Palacios' group at Imperial College London.
#                 Prof. Edwin van der Weide's group at the University of Twente.
#                 Prof. Vincent Terrapon's group at the University of Liege.
#
# Copyright (C) 2012-2016 SU2, the open-source CFD code.
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
#
################################################################################

check_PROGRAMS  =
TESTS  =

if BUILD_TESTS

check_PROGRAMS += viscous_proj_flux
check_PROGRAMS += viscous_residual_test
check_PROGRAMS += hybrid_rdelta_test
check_PROGRAMS += load_hybrid_constants_test
if BUILD_LAPACK
check_PROGRAMS += eigensolver_test
check_PROGRAMS += zeta_unit_cell_test
check_PROGRAMS += zeta_transform_test
endif

# Eventually these will likely be required
AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/Common/include
AM_CPPFLAGS += -I$(top_srcdir)/Common/test
AM_CPPFLAGS += -I$(top_srcdir)/SU2_CFD/include
AM_CPPFLAGS += @su2_externals_INCLUDES@ @LAPACK_CXX@ -DBUILD_TESTS @BOOST_CXX@

#AM_LDFLAGS = 
#LIBS = 

LDADD  =
LDADD += $(top_builddir)/SU2_CFD/obj/libSU2Core.a
LDADD += $(top_builddir)/Common/lib/libSU2.a
LDADD += @su2_externals_LIBS@ @LAPACK_LD@
if BUILD_LAPACK
LDFLAGS  += -mkl
endif

# Test sources
viscous_proj_flux_SOURCES              = viscous_proj_flux.cpp
viscous_residual_test_SOURCES          = viscous_residual_test.cpp
hybrid_rdelta_test_SOURCES             = hybrid_rdelta_test.cpp
load_hybrid_constants_test_SOURCES     = load_hybrid_constants_test.cpp
if BUILD_LAPACK
eigensolver_test_SOURCES               = eigensolver_test.cpp
zeta_unit_cell_test_SOURCES            = zeta_unit_cell_test.cpp
zeta_transform_test_SOURCES            = zeta_transform_test.cpp
endif

# Tests to be run
TESTS += viscous_proj_flux
TESTS += viscous_residual_test
TESTS += hybrid_rdelta_test
TESTS += load_hybrid_constants_test
if BUILD_LAPACK
TESTS += eigensolver_test
TESTS += zeta_unit_cell_test
TESTS += zeta_transform_test
endif

# Any of above we expect to fail?
#XFAIL_TESTS =
endif

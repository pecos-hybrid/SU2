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

check_PROGRAMS += load_hybrid_constants_test
if BUILD_LAPACK
check_PROGRAMS += eigensolver_test
check_PROGRAMS += zeta_unit_cell_test
check_PROGRAMS += zeta_transform_test
endif
if HAVE_BOOST
check_PROGRAMS += hybrid_rdelta_test
check_PROGRAMS += fluctuating_stress_test
check_PROGRAMS += forcing_test
check_PROGRAMS += load_heterogeneous_restart
check_PROGRAMS += viscous_ideal_vs_general
check_PROGRAMS += viscous_proj_flux
check_PROGRAMS += viscous_residual_test
check_PROGRAMS += viscous_residual_2d_test
check_PROGRAMS += viscous_model_split
endif

# Eventually these will likely be required
AM_CPPFLAGS =
AM_CPPFLAGS += -I$(top_srcdir)/Common/include
AM_CPPFLAGS += -I$(top_srcdir)/Common/test
AM_CPPFLAGS += -I$(top_srcdir)/SU2_CFD/include
AM_CPPFLAGS += @su2_externals_INCLUDES@ @LAPACK_CXX@ @BOOST_CXX@

#AM_LDFLAGS = 
#LIBS = 

LDADD  =
LDADD += $(top_builddir)/SU2_CFD/obj/libSU2Core.a
LDADD += $(top_builddir)/Common/lib/libSU2.a
#LDADD += @su2_externals_LIBS@ @LAPACK_LD@
LDADD += @su2_externals_LIBS@ @LAPACK_LDALL@
if BUILD_LAPACK
#LDFLAGS  += -mkl
LDFLAGS  += @LAPACK_LDALL@
endif

# Test sources
load_hybrid_constants_test_SOURCES     = load_hybrid_constants_test.cpp
if BUILD_LAPACK
eigensolver_test_SOURCES               = eigensolver_test.cpp
zeta_unit_cell_test_SOURCES            = zeta_unit_cell_test.cpp
zeta_transform_test_SOURCES            = zeta_transform_test.cpp
endif
if HAVE_BOOST
hybrid_rdelta_test_SOURCES             = hybrid_rdelta_test.cpp
fluctuating_stress_test_SOURCES        = fluctuating_stress_test.cpp
forcing_test_SOURCES                   = forcing_test.cpp
load_heterogeneous_restart_SOURCES     = load_heterogeneous_restart.cpp
viscous_ideal_vs_general_SOURCES       = viscous_ideal_vs_general.cpp
viscous_proj_flux_SOURCES              = viscous_proj_flux.cpp
viscous_residual_test_SOURCES          = viscous_residual_test.cpp
viscous_residual_2d_test_SOURCES       = viscous_residual_2d_test.cpp
viscous_model_split_SOURCES            = viscous_model_split.cpp
endif

# Tests to be run
TESTS += load_hybrid_constants_test
if BUILD_LAPACK
TESTS += eigensolver_test
TESTS += zeta_unit_cell_test
TESTS += zeta_transform_test
endif
if HAVE_BOOST
TESTS += hybrid_rdelta_test
TESTS += fluctuating_stress_test
TESTS += forcing_test
TESTS += load_heterogeneous_restart
TESTS += viscous_ideal_vs_general
TESTS += viscous_proj_flux
TESTS += viscous_residual_test
TESTS += viscous_residual_2d_test
TESTS += viscous_model_split
endif

# Any of above we expect to fail?
#XFAIL_TESTS =

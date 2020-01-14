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

if HAVE_BOOST_UTF
check_PROGRAMS += test_driver
endif

# Eventually these will likely be required
AM_CPPFLAGS = -I$(top_srcdir)/Common/include @LAPACK_CXX@ @BOOST_CPPFLAGS@
LDADD = $(top_builddir)/Common/lib/libSU2.a @LAPACK_LD@ @BOOST_UNIT_TEST_FRAMEWORK_LIB@
LDADD += @BOOST_LDFLAGS@ 
if BUILD_LAPACK
LDFLAGS  += -mkl
endif

# always link to built dependencies from ./externals
AM_CPPFLAGS += @su2_externals_INCLUDES@
LDADD += @su2_externals_LIBS@

# Test sources
if HAVE_BOOST_UTF
test_driver_SOURCES = unit_test_config.cpp \
			     main.cpp
if BUILD_LAPACK
test_driver_SOURCES += resolution_tensor_test.cpp
endif
endif

# Tests to be run
if HAVE_BOOST_UTF
TESTS += test_driver
endif

# Any of above we expect to fail?
#XFAIL_TESTS =

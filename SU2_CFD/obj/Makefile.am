################################################################################
#
# \file Makefile.am
# \brief Makefile for SU2_CFD
# \author M. Colonno, T. Economon, F. Palacios
# \version 6.2.0 "Falcon"
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
#
################################################################################

AUTOMAKE_OPTIONS = subdir-objects
ACLOCAL_AMFLAGS = ${ACLOCAL_FLAGS}

bin_PROGRAMS =
noinst_LIBRARIES =

if BUILD_NORMAL
bin_PROGRAMS += ../bin/SU2_CFD
noinst_LIBRARIES+= libSU2Core.a
endif

if BUILD_DIRECTDIFF
bin_PROGRAMS += ../bin/SU2_CFD_DIRECTDIFF
noinst_LIBRARIES+= libSU2Core_DIRECTDIFF.a
endif

if BUILD_REVERSE
bin_PROGRAMS += ../bin/SU2_CFD_AD
noinst_LIBRARIES+= libSU2Core_AD.a
endif

libSU2Core_sources = \
  ../include/data_manufactured_solutions.hpp \
  ../include/definition_structure.hpp \
  ../include/fluid_model.hpp \
  ../include/fluid_model.inl \
  ../include/integration_structure.hpp \
  ../include/integration_structure.inl \
  ../include/driver_structure.hpp \
  ../include/iteration_structure.hpp \
  ../include/numerics_structure.hpp \
  ../include/numerics_direct_mean_hybrid.hpp \
  ../include/numerics_structure.inl \
  ../include/numerics_direct_mean_hybrid.inl \
  ../include/numerics_structure_v2f.hpp \
  ../include/output_structure.hpp \
  ../include/output_structure.inl \
	../include/slice_file_reader.hpp \
  ../include/sgs_model.hpp \
  ../include/sgs_model.inl \
  ../include/solver_structure.hpp \
  ../include/solver_structure.inl \
  ../include/solver_structure_v2f.hpp \
  ../include/solver_structure_v2f.inl \
  ../include/fluctuating_stress.hpp \
  ../include/fluctuating_stress.inl \
  ../include/hybrid_RANS_LES_forcing.hpp \
  ../include/hybrid_RANS_LES_forcing.inl \
  ../include/hybrid_RANS_LES_model.hpp \
  ../include/SU2_CFD.hpp \
  ../include/task_definition.hpp \
  ../include/task_definition.inl \
  ../include/transport_model.hpp \
  ../include/transport_model.inl \
  ../include/variable_structure.hpp \
  ../include/variable_structure.inl \
  ../include/variable_structure_v2f.hpp \
  ../include/variable_structure_v2f.inl \
  ../include/transfer_structure.hpp \
  ../include/transfer_structure.inl \
  ../src/data_manufactured_solutions.cpp \
  ../src/definition_structure.cpp \
  ../src/fluid_model.cpp \
  ../src/fluid_model_pig.cpp \
  ../src/fluid_model_pvdw.cpp \
  ../src/fluid_model_ppr.cpp \
  ../src/fluid_model_inc.cpp \
  ../src/integration_structure.cpp \
  ../src/integration_time.cpp \
  ../src/driver_direct_multizone.cpp \
  ../src/driver_direct_singlezone.cpp \
  ../src/driver_structure.cpp \
  ../src/iteration_structure.cpp \
  ../src/numerics_adjoint_mean.cpp \
  ../src/numerics_adjoint_turbulent.cpp \
  ../src/numerics_direct_heat.cpp \
  ../src/numerics_direct_mean.cpp \
  ../src/numerics_direct_mean_hybrid.cpp \
  ../src/numerics_direct_mean_inc.cpp \
  ../src/numerics_direct_transition.cpp \
  ../src/numerics_direct_turbulent.cpp \
  ../src/numerics_direct_turbulent_v2f.cpp \
  ../src/numerics_direct_elasticity_nonlinear.cpp \
  ../src/numerics_direct_elasticity_linear.cpp \
  ../src/numerics_direct_elasticity.cpp \
  ../src/numerics_structure.cpp \
  ../src/numerics_template.cpp \
  ../src/output_cgns.cpp \
  ../src/output_structure.cpp \
  ../src/output_physics.cpp \
  ../src/output_tecplot.cpp \
  ../src/output_fieldview.cpp \
  ../src/output_su2.cpp \
  ../src/output_paraview.cpp \
  ../src/python_wrapper_structure.cpp \
	../src/slice_file_reader.cpp \
  ../src/solver_adjoint_mean.cpp \
  ../src/solver_adjoint_turbulent.cpp \
  ../src/solver_adjoint_discrete.cpp \
  ../src/solver_adjoint_elasticity.cpp \
  ../src/solver_direct_heat.cpp \
  ../src/solver_direct_mean.cpp \
  ../src/solver_direct_mean_fem.cpp \
  ../src/solver_direct_mean_inc.cpp \
  ../src/solver_direct_transition.cpp \
  ../src/solver_direct_turbulent.cpp \
  ../src/solver_direct_turbulent_v2f.cpp \
  ../src/solver_direct_elasticity.cpp \
  ../src/solver_structure.cpp \
  ../src/solver_template.cpp \
  ../src/fluctuating_stress.cpp \
  ../src/hybrid_RANS_LES_forcing.cpp \
  ../src/hybrid_RANS_LES_model.cpp \
  ../src/transfer_physics.cpp \
  ../src/transfer_structure.cpp \
  ../src/transport_model.cpp \
  ../src/variable_adjoint_mean.cpp \
  ../src/variable_adjoint_turbulent.cpp \
  ../src/variable_adjoint_discrete.cpp \
  ../src/variable_adjoint_elasticity.cpp \
  ../src/variable_direct_heat.cpp \
  ../src/variable_direct_mean.cpp \
  ../src/variable_direct_mean_inc.cpp \
  ../src/variable_direct_transition.cpp \
  ../src/variable_direct_turbulent.cpp \
  ../src/variable_direct_turbulent_v2f.cpp \
  ../src/variable_direct_elasticity.cpp \
  ../src/variable_structure.cpp \
  ../src/variable_template.cpp

su2_cfd_sources = \
  ../include/SU2_CFD.hpp \
  ../src/SU2_CFD.cpp

libSU2Core_cxx_flags = -fPIC
libSU2Core_libadd = 

su2_cfd_cxx_flags =
su2_cfd_ldadd =

# always link to built dependencies from ./externals
su2_cfd_cxx_flags += @su2_externals_INCLUDES@
su2_cfd_ldadd += @su2_externals_LIBS@
su2_cfd_ldadd += @su2_externals_LIBPTHREAD@
libSU2Core_cxx_flags += @su2_externals_INCLUDES@

# if BUILD_MUTATIONPP
su2_cfd_cxx_flags += @MUTATIONPP_CXX@
su2_cfd_ldadd += @MUTATIONPP_LD@
# endif

# if BUILD_LAPACK
libSU2Core_cxx_flags += @LAPACK_CXX@
libSU2Core_libadd += @LAPACK_LD@
su2_cfd_cxx_flags += @LAPACK_CXX@
#su2_cfd_ldadd += @LAPACK_LD@
su2_cfd_ldadd += @LAPACK_LDALL@
# endif

if BUILD_NORMAL
libSU2Core_a_SOURCES = $(libSU2Core_sources)
libSU2Core_a_CXXFLAGS = $(libSU2Core_cxx_flags)
libSU2Core_a_LIBADD = $(libSU2Core_libadd)
___bin_SU2_CFD_SOURCES = $(su2_cfd_sources)
___bin_SU2_CFD_CXXFLAGS = ${su2_cfd_cxx_flags}
___bin_SU2_CFD_LDADD = libSU2Core.a ../../Common/lib/libSU2.a ${su2_cfd_ldadd}
endif

if BUILD_DIRECTDIFF
libSU2Core_DIRECTDIFF_a_SOURCES = $(libSU2Core_sources)
libSU2Core_DIRECTDIFF_a_CXXFLAGS = @DIRECTDIFF_CXX@ $(libSU2Core_cxx_flags)
libSU2Core_DIRECTDIFF_a_LIBADD = @DIRECTDIFF_LIBS@ $(libSU2Core_libadd)
___bin_SU2_CFD_DIRECTDIFF_SOURCES = $(su2_cfd_sources) 
___bin_SU2_CFD_DIRECTDIFF_CXXFLAGS = @DIRECTDIFF_CXX@ ${su2_cfd_cxx_flags}
___bin_SU2_CFD_DIRECTDIFF_LDADD = libSU2Core_DIRECTDIFF.a @DIRECTDIFF_LIBS@ ../../Common/lib/libSU2_DIRECTDIFF.a ${su2_cfd_ldadd}
endif

if BUILD_REVERSE
libSU2Core_AD_a_SOURCES = $(libSU2Core_sources)
libSU2Core_AD_a_CXXFLAGS = @REVERSE_CXX@ $(libSU2Core_cxx_flags)
libSU2Core_AD_a_LIBADD = @REVERSE_LIBS@ $(libSU2Core_libadd)
___bin_SU2_CFD_AD_SOURCES = $(su2_cfd_sources) 
___bin_SU2_CFD_AD_CXXFLAGS = @REVERSE_CXX@ ${su2_cfd_cxx_flags}
___bin_SU2_CFD_AD_LDADD = libSU2Core_AD.a @REVERSE_LIBS@ ../../Common/lib/libSU2_AD.a ${su2_cfd_ldadd}
endif

project('SU2', 'c', 'cpp',
        version: '6.0.1',
	license: 'LGPL2',
	default_options: ['buildtype=release',
                          'warning_level=0',
                          'cpp_std=c++11',
                          'b_ndebug=if-release'])


pymod = import('python')
python = pymod.find_installation()

su2_cpp_args = []
su2_deps     = []

default_warning_flags = ['-Wno-unused-parameter',
                         '-Wno-format-security',
                         '-Wno-deprecated-declarations',
			 '-Wno-non-virtual-dtor']

# external dependencies
mpi_dep = [dependency('mpi', language:'c', required : get_option('with-mpi')),
           dependency('mpi', language:'cpp', required : get_option('with-mpi'))]

if get_option('enable-autodiff') or get_option('enable-directdiff')
  codi_dep = [declare_dependency(include_directories: 'externals/codi/include')]
  codi_rev_args = '-DCODI_REVERSE_TYPE'
  codi_for_args = '-DCODI_FORWARD_TYPE'
endif

# Add cgns library
if get_option('enable-cgns')
  subdir('externals/cgns')
  su2_deps     += cgns_dep
  su2_cpp_args += '-DHAVE_CGNS'
endif

# Add Boost UTF dependency
if get_option('enable-boost-utf')
  boost_dep = dependency('boost', modules : ['unit_test_framework'], static : false)
else
  boost_dep = disabler()
endif

# check if MPI dependencies are found and add them
if mpi_dep[0].found() and mpi_dep[1].found()

  # add MPI dependency
  su2_deps += mpi_dep
  su2_cpp_args += '-DHAVE_MPI'

  # compile metis
  subdir('externals/metis')

  # add metis dependency
  su2_deps     += metis_dep
  su2_cpp_args += '-DHAVE_METIS'

  # compile parmetis
  subdir('externals/parmetis')

  # add parmetis dependency
  su2_deps     += parmetis_dep
  su2_cpp_args += '-DHAVE_PARMETIS'

  # add medi dependency
  if get_option('enable-autodiff') or get_option('enable-directdiff')
    codi_dep += declare_dependency(include_directories: ['externals/medi/include', 'externals/medi/src'])
  endif
endif

if get_option('enable-tecio')
  subdir('externals/tecio')
endif

# blas-type dependencies
if get_option('enable-mkl')

  su2_cpp_args += ['-DHAVE_MKL', '-DHAVE_LAPACK']

  # the following mkl name matches the linked libs in manual dependency
  # see https://software.intel.com/en-us/articles/intel-math-kernel-library-intel-mkl-and-pkg-config-tool
  # should probably switch to iomp if Intel OpenMP is used
  mkl_dep = dependency('mkl-static-lp64-seq', required: false)

  if not mkl_dep.found()
    mkl_root = get_option('mkl_root')
    mkl_dep = declare_dependency(include_directories: mkl_root+'/include', compile_args: '-m64',
      link_args: ['-L'+mkl_root+'/lib/intel64','-lmkl_intel_lp64','-lmkl_sequential','-lmkl_core','-lpthread','-lm','-ldl'])
  endif
  su2_deps += mkl_dep

else

  mkl_dep = disabler()

  if get_option('enable-openblas')

    su2_cpp_args += '-DHAVE_BLAS'
    su2_cpp_args += '-DHAVE_LAPACK'
    blas_dep = dependency(get_option('blas-name'))
    su2_deps += blas_dep

  endif

endif

# compile common library
subdir('Common/src')
# Common unit tests
subdir('Common/test')
# compile SU2_CFD executable
subdir('SU2_CFD/src')
# SU2_CFD unit tests
subdir('SU2_CFD/test')
# compile SU2_DOT executable
subdir('SU2_DOT/src')
# compile SU2_DEF executable
subdir('SU2_DEF/src')
# compile SU2_GEO executable
subdir('SU2_GEO/src')
# compile SU2_MSH executable
subdir('SU2_MSH/src')
# compile SU2_SOL executable
subdir('SU2_SOL/src')
# install python scripts
subdir('SU2_PY')

if get_option('enable-pywrapper')
  subdir('SU2_PY/pySU2')
endif

message('''-------------------------------------------------------------------------
         |    ___ _   _ ___                                                      |
         |   / __| | | |_  )   Release 7.0.0 'Blackbird'                         |
         |   \__ \ |_| |/ /                                                      |
         |   |___/\___//___|   Meson Configuration Summary                       |
         |                                                                       |
         -------------------------------------------------------------------------
         
         Option          Value
         ---------------------
         TecIO:          @2@
         CGNS:           @3@
         AD (reverse):   @4@
         AD (forward):   @5@
         Python Wrapper: @6@
         Intel-MKL:      @7@
         OpenBlas:       @8@
         Boost UTF:      @9@

         Please be sure to add the $SU2_HOME and $SU2_RUN environment variables,
         and update your $PATH (and $PYTHONPATH if applicable) with $SU2_RUN

         Based on the input to this configuration, add these lines to your .bashrc file:

         export SU2_RUN=@0@
         export SU2_HOME=@1@
         export PATH=$PATH:$SU2_RUN
         export PYTHONPATH=$PYTHONPATH:$SU2_RUN

         Use './ninja -C @10@ install' to compile and install SU2
'''.format(get_option('prefix')+'/bin', meson.source_root(), get_option('enable-tecio'), get_option('enable-cgns'),
           get_option('enable-autodiff'), get_option('enable-directdiff'), get_option('enable-pywrapper'), get_option('enable-mkl'),
           get_option('enable-openblas'), get_option('enable-boost-utf'), meson.build_root()))


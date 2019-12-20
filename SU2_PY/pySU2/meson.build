# check for swig
swig = find_program('swig')

# add mpi4py include
if mpi_dep[0].found() and mpi_dep[1].found()
  if get_option('python_path') != ''
    mpi4py_dir = get_option('python_path')
  else
    mpi4py_dir = run_command(python, '-c', 'import mpi4py; print(mpi4py.__path__[0])').stdout().strip()
    message ('Looking for mpi4py in ' + mpi4py_dir)
  endif
  mpi4py_include = mpi4py_dir + '/include'
else
  mpi4py_include = ''
endif

swig_gen = generator(
    swig,
    output: ['@BASENAME@.cxx'],
    arguments: su2_cpp_args +
    [ '-c++', '-python', '-I'+mpi4py_include, '-outdir', meson.current_build_dir(), '-o', './@OUTPUT@', '@INPUT@'],
    depfile: '@BASENAME@.d',
)

swig_sources = [
    'pySU2.i',
]

wrapper_deps = [
    dependency('python3'), 
]

foreach swig_source: swig_sources
    cpp_source = swig_gen.process(swig_source)
    swig_lib = shared_library(
       '_pysu2',
       cpp_source,
       dependencies: [wrapper_deps, common_dep, su2_deps, su2_cfd_dep],
       objects: su2_cfd.extract_all_objects(),
       install: true,
       include_directories : mpi4py_include,
       cpp_args : [default_warning_flags,su2_cpp_args],
       name_prefix : '',
       install_dir: 'bin'
    )
endforeach

meson.add_install_script('install.sh', meson.current_build_dir(), 'pysu2.py')

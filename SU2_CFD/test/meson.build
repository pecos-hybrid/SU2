test_cpp_args = ['-fPIC'] + [default_warning_flags, su2_cpp_args]

if mkl_dep.found()
  eigensolver_test = executable('eigensolver_test',
                   'eigensolver_test.cpp',
                   install : false,
                   dependencies : [su2_cfd_dep, common_dep, su2_deps],
                   cpp_args: test_cpp_args)
  test('eigensolver test', eigensolver_test)
endif

SU2_CFD_boost_tests = files(['convective_blending_test.cpp',
                             'hybrid_rdelta_test.cpp',
                             'slice_test.cpp',
                             'fluctuating_stress_test.cpp',
                             'forcing_test.cpp',
                             'load_heterogeneous_restart.cpp',
                             'viscous_ideal_vs_general.cpp',
                             'viscous_proj_flux.cpp',
                             'viscous_residual_test.cpp',
                             'viscous_residual_2d_test.cpp',
                             'viscous_model_split.cpp',
                             'averaging_timescale_test.cpp',
                             'main.cpp'])

if mkl_dep.found()
  SU2_CFD_boost_tests += files(['resolution_adequacy.cpp'])
endif

SU2_CFD_test_driver = executable('test_driver',
                                 SU2_CFD_boost_tests,
                                 install : false,
                                 dependencies : [su2_cfd_dep, common_dep, su2_deps, boost_dep],
                                 cpp_args: test_cpp_args)
test('SU2_CFD Boost tests', SU2_CFD_test_driver)

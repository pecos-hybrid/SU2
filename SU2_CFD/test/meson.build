test_cpp_args = ['-fPIC'] + [default_warning_flags, su2_cpp_args]
load_hybrid_constants_test = executable('load_hybrid_constants_test',
                 'load_hybrid_constants_test.cpp',
                 install : false,
                 dependencies : [su2_cfd_dep, common_dep, su2_deps],
                 cpp_args: test_cpp_args)
test('load hybrid constants', load_hybrid_constants_test)

boost_tests = files(['convective_blending_test.cpp',
                     'hybrid_rdelta_test.cpp',
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

boost_test_runner = executable('boost_test_runner',
                               boost_tests,
                               install : false,
                               dependencies : [su2_cfd_dep, common_dep, su2_deps, boost_dep],
                               cpp_args: test_cpp_args)
test('boost tests', boost_test_runner)

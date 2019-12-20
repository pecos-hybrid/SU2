test_cpp_args = ['-fPIC'] + [default_warning_flags, su2_cpp_args]

Common_boost_tests = files(['unit_test_config.cpp',
                            'main.cpp'])
if mkl_dep.found()
  Common_boost_tests += files(['resolution_tensor_test.cpp'])
endif

Common_test_driver = executable('test_driver',
                                Common_boost_tests,
                                install : false,
                                dependencies : [common_dep, su2_deps, boost_dep],
                                cpp_args: test_cpp_args)
test('Common Boost tests', Common_test_driver)

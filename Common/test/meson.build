test_cpp_args = ['-fPIC'] + [default_warning_flags, su2_cpp_args]

Common_boost_tests = files(['unit_test_config.cpp',
                            'main.cpp'])

Common_test_driver = executable('test_driver',
                                Common_boost_tests,
                                install : false,
                                dependencies : [common_dep, su2_deps, boost_dep],
                                cpp_args: test_cpp_args)
test('Common Boost tests', Common_test_driver)

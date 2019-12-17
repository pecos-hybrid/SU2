common_include = include_directories('./')
common_src =files(['primal_grid_structure.cpp',
           'config_structure.cpp',
           'dual_grid_structure.cpp',
           'element_structure.cpp',
           'geometry_structure.cpp',
           'ad_structure.cpp',
           'grid_movement_structure.cpp',
           'grid_adaptation_structure.cpp',
           'element_linear.cpp',
           'gauss_structure.cpp',
           'adt_structure.cpp',
           'interpolation_structure.cpp',
           'mpi_structure.cpp',
           'linear_solvers_structure.cpp',
           'linear_solvers_structure_b.cpp',
           'vector_structure.cpp',
           'matrix_structure.cpp'])

if get_option('enable-normal')
  common = static_library('SU2Common', 
                          common_src,
                          install : false, 
	                  dependencies : su2_deps, 
	 		  cpp_args: [default_warning_flags, su2_cpp_args])
  common_dep = declare_dependency(link_with: common,
                                  include_directories : common_include)
endif

if get_option('enable-autodiff')
  
  commonAD = static_library('SU2CommonAD',
                            common_src,
                            install : false, 
			    dependencies : [su2_deps, codi_dep], 
			    cpp_args: [default_warning_flags, su2_cpp_args, codi_rev_args])

  commonAD_dep = declare_dependency(link_with: commonAD, 
                                    include_directories : common_include)

endif

if get_option('enable-directdiff')
  
  commonDD = static_library('SU2CommonDD',
                            common_src,
                            install : false, 
			    dependencies : [su2_deps, codi_dep], 
			    cpp_args: [default_warning_flags, su2_cpp_args, codi_for_args])

  commonDD_dep = declare_dependency(link_with: commonDD, 
                                    include_directories : common_include)

endif

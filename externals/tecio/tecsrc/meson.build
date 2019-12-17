tecio_cpp_flags = ['-fpermissive', '-DUSEENUM', '-DTHREED', '-DMAKEARCHIVE', '-DNO_ASSERTS', '-DNO_THIRD_PARTY_LIBS']

if (build_machine.system() == 'darwin')
  tecio_cpp_flags += ['-DDARWIN', '-DMAC64']
endif

if (build_machine.system() == 'linux')
  tecio_cpp_flags += ['-DLINUX', '-DLINUX64']
endif


tecio_include = include_directories('../', './')
tecio = static_library('tecio',
	   'TranslatedString.cpp',
	   'alloc.cpp',
	   'arrlist.cpp',
	   'auxdata.cpp',
	   'dataio.cpp',
	   'dataio4.cpp',
	   'dataset.cpp',
	   'dataset0.cpp',
	   'datautil.cpp',
	   'filestream.cpp',
	   'geom2.cpp',
	   'q_msg.cpp',
	   'q_unicode.cpp',
	   'set.cpp',
	   'strlist.cpp',
	   'strutil.cpp',
	   'tassert.cpp',
	   'tecxxx.cpp',
           install : false, 
	   cpp_args: [default_warning_flags, tecio_cpp_flags],
	   include_directories: tecio_include)

tecio_dep = declare_dependency(link_with: tecio,
                               include_directories: tecio_include)

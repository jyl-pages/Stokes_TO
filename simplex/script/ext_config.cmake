#ext libs that only contain header files. Just include them.
set(ext_header_lib_lis amgcl autodiff eigen nanoflann nlopt perlin_noise physbam rapidjson)
#ext libs that contain header and source files, in a single folder.
#We need to compile them in the exactly same way as src libs under COMMON mode.
set(ext_uncompiled_lib_lis imgui mma stb tetgen tiny_obj_loader tri2d)

macro(setup_ext_target)
	add_lib_target(ext ${ext_lib_path})
	pass_resources_to_global()
	message("")
endmacro(setup_ext_target)

macro(install_all_ext_libs)
	use_ext_libs("COMPILE" "${ext_uncompiled_lib_lis}")
endmacro(install_all_ext_libs)

#mode is one of LINK/COMPILE
macro(use_ext_libs mode lib_lis)
	set(lib_lis ${lib_lis})
	if(${mode} STREQUAL "LINK")
		set(ENABLE_COMIPLE OFF)
		set(ENABLE_LINK ON)
	elseif(${mode} STREQUAL "COMPILE")
		set(ENABLE_COMIPLE ON)
		set(ENABLE_LINK ON)
	else()
		message(FATAL_ERROR "use_ext_libs unrecognized mode ${mode}")
	endif()
	#ext libs that only consist of header files, just include them.
	#Like: eigen, physbam
	foreach(lib_name ${ext_header_lib_lis})
		if(${lib_name} IN_LIST lib_lis)
			message(STATUS "Use ext lib by including headers: ext/${lib_name}")
			target_add_incs(${simplex_path}/ext/${lib_name})
		endif()
	endforeach()

	#ext libs that we compile our own .lib files.
	#Link these libs, just like linking src libs in COMMON mode
	#The difference between them and src libs, however,
	#is that they will also be linked as .lib in UNIX, instead of source compilation.
	#Like: tiny_obj_loader
	foreach(lib_name ${ext_uncompiled_lib_lis})
		if(${lib_name} IN_LIST lib_lis)
			message(STATUS "Use ext lib by linking compiled .lib files by proj/_install: ext/${lib_name}")
			target_add_incs(${simplex_path}/ext/${lib_name})
			target_add_libs(${ext_lib_path} ${lib_name})
			set(target_lib_build_path ${ext_lib_path}/${lib_name}/cmake-build-${platform})
			target_install_dep(${simplex_path}/ext/${lib_name} ${target_lib_build_path} ${lib_name})
		endif()
	endforeach()

  #openmp
  if("openmp" IN_LIST lib_lis)
    message(STATUS "Use ext lib: OpenMP")
    find_package(OpenMP REQUIRED)
  	if(OPENMP_FOUND)
      message(STATUS "openmp add flags: ${OpenMP_CXX_FLAGS}")
  		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
  		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
    else()
      message(FATAL_ERROR "OpenMP not found")
  	endif()
  endif()

  #opengl
  #This will add support for opengl, glew and glm
  if("opengl" IN_LIST lib_lis)
    message(STATUS "Use ext libs: OpenGL, GLEW and GLM")
    if(WIN32)
			target_add_incs(${simplex_path}/ext/freeglut/include)
			set(freeglut_lib_path ${simplex_path}/ext/freeglut/lib/x64)
			set(freeglut_bin_path ${simplex_path}/ext/freeglut/bin/x64)
  		set(glew_libs ${freeglut_lib_path}/glew32.lib)
  		set(glut_libs debug ${freeglut_lib_path}/freeglutd.lib optimized ${freeglut_lib_path}/freeglut.lib)
			list(APPEND ${current_target}_libs ${glew_libs} ${glut_libs})
			list(APPEND ${current_target}_dlls ${freeglut_bin_path}/freeglut.dll ${freeglut_bin_path}/glew32.dll)
  		#list(APPEND lib_files ${glew_libs} ${glut_libs} PARENT_SCOPE)
  	elseif(UNIX)#freeglut and glew are installed on linux by "sudo apt-get install freeglut3-dev libglew-dev"
  		set(GCC_COVERAGE_COMPILE_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS} -lGL -lglut -lGLU -lGLEW")
  		set(GCC_COVERAGE_LINK_FLAGS "${GCC_COVERAGE_LINK_FLAGS} -lGL -lglut -lGLU -lGLEW")
  		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/lib/x86_64-linux-gnu ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
  		set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
  	endif(WIN32)
		target_add_incs(${simplex_path}/ext/glm)
  endif()

	#ipopt
	if("ipopt" IN_LIST lib_lis)
		if(WIN32)
			message(STATUS "Use ext lib: ipopt")
			add_definitions(-DUSE_IPOPT)
			target_add_incs(${simplex_path}/ext/ipopt/include/coin)
			target_add_libs(${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL IpOptFSS)
			target_add_libs(${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL Ipopt-vc8)
			target_add_libs(${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL Ipopt-vc8)
			list(APPEND ${current_target}_dlls ${simplex_path}/ext/ipopt/lib/x64/ReleaseMKL/Ipopt-vc8.dll)
		elseif(UNIX)
			message(STATUS "Ignore lib ipopt in Linux")
		endif()
	endif()

	#libtorch
	if("libtorch" IN_LIST lib_lis)
		message(STATUS "Use ext lib: libtorch")
		set(CMAKE_PREFIX_PATH ${simplex_path}/ext/libtorch/share/cmake/Torch)
		find_package(Torch REQUIRED)
		list(APPEND ${current_target}_libs ${TORCH_LIBRARIES})
	endif()

	if("cuda" IN_LIST lib_lis)
		message(STATUS "Support CUDA Language")

		set(USE_CUDA ON)
		add_definitions(-DUSE_CUDA)
		set(USE_CPX ON)
		add_definitions(-DUSE_CPX)

		enable_language(CUDA)
		find_package(CUDAToolkit)
		target_add_incs(${simplex_path}/ext/cuda/CUDASamples/v11.1/inc)
		list(APPEND ${current_target}_libs CUDA::cudart CUDA::cublas CUDA::cufft CUDA::cusolver CUDA::curand)
	endif()

endmacro(use_ext_libs)

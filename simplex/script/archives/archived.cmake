macro(add_project)
  add_all_code_files("xproj/${proj_name}" ${proj_src_path} target_files)
  list(APPEND src_files ${target_files})
endmacro(add_project)

macro(add_headers_and_link_libs)
	add_simplex_headers(ext eigen) # Basic matrix like Vector3 here
  add_simplex_headers(ext physbam) # Parse_Args is here

  add_simplex_src_lib_on(USE_COMMON common)
  add_simplex_src_lib_on(USE_GEOMETRY geometry)
  add_simplex_src_lib_on(USE_SOLVER solver)
  add_simplex_src_lib_on(USE_PHYSICS physics)
  add_simplex_src_lib_on(USE_OPTIMIZATION optimization)

  if(USE_VIEWER)	#compiling the opengl viewer needs glm, and need imgui, stb optionally
		add_opengl_headers_and_link_opengl_libs()
		if(WIN32)
			compile_and_link_lib(src viewer ${simplex_lib_path})
		elseif(UNIX)
			add_headers_and_src(src viewer)
		endif(WIN32)
	endif()
  if(USE_GLM)
    add_simplex_headers(ext glm)
  endif()

  if(USE_IMGUI)
    add_simplex_headers(ext imgui)
    add_lib_to_link(${ext_lib_path} imgui)
    add_definitions(-DUSE_IMGUI)
  endif()
  if(USE_STB)
    add_simplex_headers(ext stb)
    add_lib_to_link(${ext_lib_path} stb)
    add_definitions(-DUSE_STB)
  endif()

	if(USE_CUDA)
		add_cuda()
		add_headers_and_src(src parallel)
		add_definitions(-DUSE_CUDA)
	endif(USE_CUDA)

	if(USE_OPENMP)
		add_openmp()
	endif(USE_OPENMP)

	if(USE_TINY_OBJ_LOADER)
		add_simplex_headers(ext tiny_obj_loader)
		add_lib_to_link(${ext_lib_path} tiny_obj_loader)
		add_definitions(-DUSE_TINY_OBJ_LOADER)
	endif()

	if(USE_RAPIDJSON)
		add_simplex_headers(ext rapidjson)
		add_definitions(-DUSE_RAPIDJSON)
	endif()

	if(USE_IPOPT)
		add_headers_and_src(ext ipopt/include/coin)
		add_lib_to_link(${ext_path}/ipopt/lib/x64/ReleaseMKL IpOptFSS)
		add_lib_to_link(${ext_path}/ipopt/lib/x64/ReleaseMKL Ipopt-vc8)
		add_lib_to_link(${ext_path}/ipopt/lib/x64/ReleaseMKL Ipopt-vc10)
		add_definitions(-DUSE_IPOPT)
		add_definitions(-DIPOPT_DLL)
	endif()

	if(USE_NLOPT)
		add_headers_and_src(ext nlopt)
		add_definitions(-DUSE_NLOPT)
	endif()

	if(USE_NANOFLANN)
		#add_headers_and_src(ext nanoflann)
    add_simplex_headers(ext nanoflann)
		add_definitions(-DUSE_NANOFLANN)
	endif()

	if(USE_MMA)
		add_headers_and_src(ext mma)
		add_definitions(-DUSE_MMA)
	endif()

	if(USE_LIBTORCH)
		add_libtorch()
	endif()

	if(USE_TRI2D)
		add_simplex_headers(ext tri2d)
		add_lib_to_link(${ext_lib_path} tri2d)
		add_definitions(-DUSE_TRI2D)
	endif()

	if(USE_TETGEN)
		add_simplex_headers(ext tetgen)
		add_lib_to_link(${ext_lib_path} tetgen)
		add_definitions(-DUSE_TETGEN)
	endif()

	if(USE_AUTODIFF)
		add_simplex_headers(ext autodiff)
		add_definitions(-DUSE_AUTODIFF)
	endif()

	if(USE_AMGCL)
		add_simplex_headers(ext amgcl)

		add_definitions(-DUSE_AMGCL)
		## use boost
		# add_simplex_headers(ext boost)
		## use no boost
		add_definitions(-DAMGCL_NO_BOOST)
	endif()

	if(USE_PERLIN_NOISE)
		add_headers_and_src(ext perlin_noise)
		add_definitions(-DUSE_PERLIN_NOISE)
	endif()

endmacro(add_headers_and_link_libs)

macro(link_project)

	Print_List(ListOption)
	add_headers_and_link_libs()

	if(USE_CUDA)
		add_executable(${proj_name} ${src_files})
		#set_target_properties(${proj_name} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
		#cuda_add_executable(${proj_name} ${src_files} OPTIONS ${cuda_options})						#nvcc dbg: set flag here does not work
		#target_compile_options(${proj_name} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:--std=c++17 >)	 	#nvcc dbg: set flag here does not work either?
	else(USE_CUDA)
		add_executable(${proj_name} ${src_files})
	endif(USE_CUDA)

	if(USE_OPENMP AND OpenMP_CXX_FOUND)
		message("OPENMP found")
		target_link_libraries(${proj_name} PUBLIC OpenMP::OpenMP_CXX)
	endif()

	set_target_properties(${proj_name} PROPERTIES PROJECT_LABEL xproj-${proj_name})
	message(STATUS "Add exe: " ${proj_name})

	list(LENGTH lib_files lib_n)
	if(lib_n GREATER 0)
		target_link_libraries(${proj_name} ${lib_files})
	endif()

	list(LENGTH dep_files dep_n)
	if(dep_n GREATER 0)
		add_dependencies(${proj_name} ${dep_files})
	endif()

	debug_info()
endmacro(link_project)

####CUDA
####TODO: this function needs some update after CUDA11
macro(add_cuda)
	if(USE_CUDA)
		add_definitions(-DUSE_CUDA)
	endif()
	if(WIN32)
		#find_package(CUDA QUIET)	#deprecated
		enable_language(CUDA)
		include(CheckLanguage)
		check_language(CUDA)

		if(CMAKE_CUDA_COMPILER_ID STREQUAL "NVIDIA")
		#if(CUDA_FOUND)
		#TODO: this if section derived from the previous version, need some check here
			message(STATUS "CUDA_PATH: " $ENV{CUDA_PATH})
			message(STATUS "CUDA Found")
			set(cuda_sample_path ${ext_path}/cuda/CUDASamples/v11.1/inc)
			include_directories(${cuda_sample_path})
			#add_simplex_headers(ext cusp)
			add_definitions(-DFORCE_INLINES)
			set(CUDA_COMPUTE_CAPABILITY "61" CACHE STRING "CUDA Compute Capability")

			set(cuda_libs
				$ENV{CUDA_PATH}/lib/x64/cublas.lib
				$ENV{CUDA_PATH}/lib/x64/cuda.lib
				$ENV{CUDA_PATH}/lib/x64/cudadevrt.lib
				$ENV{CUDA_PATH}/lib/x64/cudart.lib
				$ENV{CUDA_PATH}/lib/x64/cudart_static.lib
				$ENV{CUDA_PATH}/lib/x64/cufft.lib
				$ENV{CUDA_PATH}/lib/x64/cufftw.lib
				$ENV{CUDA_PATH}/lib/x64/curand.lib
				$ENV{CUDA_PATH}/lib/x64/cusolver.lib
				$ENV{CUDA_PATH}/lib/x64/cusparse.lib)

			list(APPEND lib_files ${cuda_libs})

			if(NOT CMAKE_CXX_COMPILER_ID STREQUAL "Clang")
				message(STATUS "CUDA Clang setting")
				set(cuda_options "-arch=sm_${CUDA_COMPUTE_CAPABILITY}")
				list(APPEND CUDA_NVCC_FLAGS "-DFORCE_INLINES;")
				list(APPEND CUDA_NVCC_FLAGS "-D__CUDA_NO_HALF_OPERATORS__;")
				#list(APPEND lib_files ${CUDA_cusparse_LIBRARY})
			endif()
		endif()
	elseif(UNIX)
		find_package(CUDAToolkit 11.1 REQUIRED)
		if(CUDAToolkit_FOUND)
			set(CUDA_PATH /usr/local/cuda)
			#include_directories(${CUDA_PATH}/include/cusp)
			add_simplex_headers(ext cusp)
			add_definitions(-DFORCE_INLINES)
			set(cuda_sample_path ${ext_path}/cuda/CUDASamples/v11.1/inc)
			include_directories(${cuda_sample_path})
			include_directories(${CUDA_PATH}/include)
			set(cuda_libs
					${CUDA_PATH}/lib64/libcublas_static.a
					${CUDA_PATH}/lib64/libcublasLt_static.a
					${CUDA_PATH}/lib64/libculibos.a
					${CUDA_PATH}/lib64/libcudadevrt.a
					${CUDA_PATH}/lib64/libcudart_static.a
					${CUDA_PATH}/lib64/libcufft_static.a
					${CUDA_PATH}/lib64/libcufftw_static.a
					${CUDA_PATH}/lib64/libcurand_static.a
				    ${CUDA_PATH}/lib64/libcusolver_static.a
					${CUDA_PATH}/lib64/libmetis_static.a
					${CUDA_PATH}/lib64/libcusparse_static.a)
			list(APPEND lib_files ${cuda_libs})
		endif()
	endif(WIN32)
endmacro(add_cuda)

####OPENMP
macro(add_openmp)
	#set(OpenMP_CXX "${CMAKE_CXX_COMPILER}")
	#set(OpenMP_CXX_FLAGS "-fopenmp=libomp -Wno-unused-command-line-argument")
	#set(OpenMP_CXX_LIB_NAMES "libomp" "libgomp" "libiomp5")
	#set(OpenMP_libomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
	#set(OpenMP_libgomp_LIBRARY ${OpenMP_CXX_LIB_NAMES})
	#set(OpenMP_libiomp5_LIBRARY ${OpenMP_CXX_LIB_NAMES})

	find_package(OpenMP REQUIRED)
	if(OPENMP_FOUND)
		message("OPENMP found")

		set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} ${OpenMP_C_FLAGS}")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${OpenMP_CXX_FLAGS}")
		set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} ${OpenMP_EXE_LINKER_FLAGS}")
	endif()
endmacro(add_openmp)

#compile_and_link_lib: add to src_files, lib_files, and dep_files
#Then cmake will automatically require the lib_files to be compiled
macro(compile_and_link_lib parent_dir lib_name lib_output_dir)
	set(target_lib_src_path ${simplex_path}/${parent_dir}/${lib_name})	#parent_dir = ext or src
  set(target_lib_build_path ${build_path}/${lib_name}/cmake-build-${platform})
	include_directories(${target_lib_src_path})#include .h files only
	message(STATUS "target_lib_build_path: " ${target_lib_build_path})
	add_subdirectory(${target_lib_src_path} ${target_lib_build_path})	#two parameters are src path and build path
  add_lib_to_link(${lib_output_dir} ${lib_name})
	list(APPEND dep_files ${lib_name})
endmacro(compile_and_link_lib)


####called in a subdirectory
macro(compile_ext_lib target_lib_name)
	compile_lib(ext ${target_lib_name} ${ext_lib_path})
endmacro(compile_ext_lib)

macro(add_opengl_headers_and_link_opengl_libs)
	if(WIN32)
		set(freeglut_src_path ${ext_path}/freeglut/include)
		set(freeglut_lib_path ${ext_path}/freeglut/lib/x64)
		include_directories(${freeglut_src_path})
		set(glew_libs ${freeglut_lib_path}/glew32.lib)
		list(APPEND lib_files ${glew_libs})
		set(glut_libs debug ${freeglut_lib_path}/freeglutd.lib optimized ${freeglut_lib_path}/freeglut.lib)
		list(APPEND lib_files ${glut_libs})
	elseif(UNIX)#freeglut and glew are installed on linux by "sudo apt-get install freeglut3-dev libglew-dev"
		set(GCC_COVERAGE_COMPILE_FLAGS "${GCC_COVERAGE_COMPILE_FLAGS} -lGL -lglut -lGLU -lGLEW")
		set(GCC_COVERAGE_LINK_FLAGS "${GCC_COVERAGE_LINK_FLAGS} -lGL -lglut -lGLU -lGLEW")
		set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -L/usr/lib/x86_64-linux-gnu ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
		set(CMAKE_CXX_STANDARD_LIBRARIES "${CMAKE_CXX_FLAGS} ${GCC_COVERAGE_COMPILE_FLAGS} ${GCC_COVERAGE_LINK_FLAGS}")
	endif(WIN32)
endmacro(add_opengl_headers_and_link_opengl_libs)

macro(add_libtorch)
	set(CMAKE_PREFIX_PATH ${ext_path}/libtorch/share/cmake/Torch)
	find_package(Torch REQUIRED)
	list(APPEND lib_files ${TORCH_LIBRARIES})
endmacro(add_libtorch)

####add headers only
macro(add_simplex_headers parent_folder_name header_folder_name)
	include_directories(${simplex_path}/${parent_folder_name}/${header_folder_name})
endmacro(add_simplex_headers)

#Like: add_headers_and_src(src geometry)
macro(add_headers_and_src folder_name target_name)
  add_all_code_files("${folder_name}/${target_name}" ${simplex_path}/${folder_name}/${target_name} target_files)
  list(APPEND src_files ${target_files})
endmacro(add_headers_and_src)

macro(add_rt_headers_and_src rt_path folder_name target_name)
  add_all_code_files("${folder_name}/${target_name}" ${rt_path}/${folder_name}/${target_name} target_files)
  list(APPEND src_files ${target_files})
endmacro(add_rt_headers_and_src)

####this function is called only in simplex/projects/_install/CMakeLists.txt
macro(install_ext_libs)
	#It will call multiple compile_and_link_lib()
	#opengl (needed by imgui)
	add_opengl_headers_and_link_opengl_libs()
	compile_and_link_lib(ext stb ${ext_lib_path})
	add_definitions(-DUSE_STB)
	compile_and_link_lib(ext imgui ${ext_lib_path})
	add_definitions(-DUSE_IMGUI)
	compile_and_link_lib(ext tiny_obj_loader ${ext_lib_path})
	add_definitions(-DUSE_TINY_OBJ_LOADER)
	compile_and_link_lib(ext tri2d ${ext_lib_path})
	add_definitions(-DUSE_TRI2D)
	compile_and_link_lib(ext tetgen ${ext_lib_path})
	add_definitions(-DUSE_TETGEN)
endmacro(install_ext_libs)

####this function is called only in simplex/projects/_install/CMakeLists.txt
macro(install_simplex_libs)
	add_opengl_headers_and_link_opengl_libs()
	compile_and_link_lib(src common ${simplex_lib_path})
	compile_and_link_lib(src geometry ${simplex_lib_path})
	compile_and_link_lib(src solver ${simplex_lib_path})
	compile_and_link_lib(src optimization ${simplex_lib_path})
	compile_and_link_lib(src physics ${simplex_lib_path})
	compile_and_link_lib(src viewer ${simplex_lib_path})
endmacro(install_simplex_libs)

macro(init_project)
	initialize_global_variables()
	set_default_options()
endmacro(init_project)

macro(compile_src_lib target_lib_name)
	add_simplex_headers(ext eigen)
	add_simplex_headers(ext physbam)
	compile_lib(src ${target_lib_name} ${simplex_lib_path})
endmacro(compile_src_lib)

macro(compile_lib lib_src_dir lib_name lib_output)
	message(STATUS "Compile lib: " ${lib_name})
	set(target_lib_src_path ${simplex_path}/${lib_src_dir}/${lib_name})
  include_directories(${target_lib_src_path})
  all_code_files_under(${target_lib_src_path} lib_src_files)
	source_group("${lib_src_dir}/${lib_name}" FILES ${lib_src_files})
	add_library(${lib_name} STATIC ${lib_src_files})

	####add postfix d to debug libs
	set_target_properties(${lib_name} PROPERTIES DEBUG_POSTFIX d)
	####rename target in IDE with path prefix
	set_target_properties(${lib_name} PROPERTIES PROJECT_LABEL ${lib_src_dir}-${lib_name})
	message(STATUS "lib_output: " ${lib_output})

	if(USE_CUDA)
		add_cuda()
		add_headers_and_src(src parallel)
	endif(USE_CUDA)

	if(USE_OPENMP)
		add_openmp()
		#find_package(OpenMP REQUIRED)
		#if(OpenMP_CXX_FOUND)
		#	message(STATUS "link openmp to " ${lib_output})
		#	target_link_libraries(${lib_name} PUBLIC OpenMP::OpenMP_CXX)
		#endif()
	endif()

	####put debug and release libs in the same folder
	set_target_properties(${lib_name} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${lib_output}")
	set_target_properties(${lib_name} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${lib_output}")
	set_target_properties(${lib_name} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${lib_output}")

	unset(target_lib_src_path)
endmacro(compile_lib)

macro(add_headers_and_link_libs)
	####Eigen
	add_simplex_headers(ext eigen)
	####Physbam
	add_simplex_headers(ext physbam)
	if(USE_COMMON)
		if(USE_COMMON_COMPILED)
			add_simplex_headers(src common)
			add_lib_to_link(${simplex_lib_path} common)
		else(USE_COMMON_COMPILED)
			compile_and_link_lib(src common ${simplex_lib_path})
		endif(USE_COMMON_COMPILED)
	endif(USE_COMMON)

	if(USE_GEOMETRY)
		if(WIN32)		#Windows supports linking multiple libs
			if(USE_GEOMETRY_COMPILED)
				add_simplex_headers(src geometry)
				add_lib_to_link(${simplex_lib_path} geometry)
			else(USE_GEOMETRY_COMPILED)
				compile_and_link_lib(src geometry ${simplex_lib_path})
			endif(USE_GEOMETRY_COMPILED)
		elseif(UNIX)	#but Linux doesn't, so src files have to be compiled
			add_headers_and_src(src geometry)
		endif(WIN32)
	endif(USE_GEOMETRY)

	if(USE_SOLVER)
		if(WIN32)
			if(USE_SOLVER_COMPILED)
				add_simplex_headers(src solver)
				add_lib_to_link(${simplex_lib_path} solver)
			else(USE_SOLVER_COMPILED)
				compile_and_link_lib(src solver ${simplex_lib_path})
			endif(USE_SOLVER_COMPILED)
		elseif(UNIX)
			add_headers_and_src(src solver)
		endif(WIN32)
	endif(USE_SOLVER)

	if(USE_PHYSICS)
		if(WIN32)
			if(USE_PHYSICS_COMPILED)
				add_simplex_headers(src physics)
				add_lib_to_link(${simplex_lib_path} physics)
			else(USE_PHYSICS_COMPILED)
				compile_and_link_lib(src physics ${simplex_lib_path})
			endif(USE_PHYSICS_COMPILED)
		elseif(UNIX)
			add_headers_and_src(src physics)
		endif(WIN32)
	endif(USE_PHYSICS)

	if(USE_OPTIMIZATION)
		if(WIN32)
			if(USE_OPTIMIZATION_COMPILED)
				add_simplex_headers(src optimization)
				add_lib_to_link(${simplex_lib_path} optimization)
			else(USE_OPTIMIZATION_COMPILED)
				compile_and_link_lib(src optimization ${simplex_lib_path})
			endif(USE_OPTIMIZATION_COMPILED)
		elseif(UNIX)
			add_headers_and_src(src optimization)
		endif(WIN32)
	endif(USE_OPTIMIZATION)

	if(USE_CUDA)
		add_cuda()
		add_headers_and_src(src parallel)
		add_definitions(-DUSE_CUDA)
	endif(USE_CUDA)

	if(USE_OPENMP)
		add_openmp()
	endif(USE_OPENMP)

	if(USE_VIEWER)	#compiling the opengl viewer needs glm, and need imgui, stb optionally
		add_opengl_headers_and_link_opengl_libs()
		add_simplex_headers(ext glm)
		if(USE_IMGUI)
			add_simplex_headers(ext imgui)
			add_lib_to_link(${ext_lib_path} imgui)
			add_definitions(-DUSE_IMGUI)
		endif()
		if(USE_STB)
			add_simplex_headers(ext stb)
			add_lib_to_link(${ext_lib_path} stb)
			add_definitions(-DUSE_STB)
		endif()
		if(WIN32)
			compile_and_link_lib(src viewer ${simplex_lib_path})
		elseif(UNIX)
			add_headers_and_src(src viewer)
		endif(WIN32)
	endif()

	if(USE_TINY_OBJ_LOADER)
		add_simplex_headers(ext tiny_obj_loader)
		add_lib_to_link(${ext_lib_path} tiny_obj_loader)
		add_definitions(-DUSE_TINY_OBJ_LOADER)
	endif()

	if(USE_RAPIDJSON)
		add_simplex_headers(ext rapidjson)
		add_definitions(-DUSE_RAPIDJSON)
	endif()

	if(USE_IPOPT)
		add_headers_and_src(ext ipopt/include/coin)
		add_lib_to_link(${ext_path}/ipopt/lib/x64/ReleaseMKL IpOptFSS)
		add_lib_to_link(${ext_path}/ipopt/lib/x64/ReleaseMKL Ipopt-vc8)
		add_lib_to_link(${ext_path}/ipopt/lib/x64/ReleaseMKL Ipopt-vc10)
		add_definitions(-DUSE_IPOPT)
		add_definitions(-DIPOPT_DLL)
	endif()

	if(USE_NLOPT)
		add_headers_and_src(ext nlopt)
		add_definitions(-DUSE_NLOPT)
	endif()

	if(USE_NANOFLANN)
		add_headers_and_src(ext nanoflann)
		add_definitions(-DUSE_NANOFLANN)
	endif()

	if(USE_MMA)
		add_headers_and_src(ext mma)
		add_definitions(-DUSE_MMA)
	endif()

	if(USE_LIBTORCH)
		add_libtorch()
	endif()

	if(USE_TRI2D)
		add_simplex_headers(ext tri2d)
		add_lib_to_link(${ext_lib_path} tri2d)
		add_definitions(-DUSE_TRI2D)
	endif()

	if(USE_TETGEN)
		add_simplex_headers(ext tetgen)
		add_lib_to_link(${ext_lib_path} tetgen)
		add_definitions(-DUSE_TETGEN)
	endif()

	if(USE_AUTODIFF)
		add_simplex_headers(ext autodiff)
		add_definitions(-DUSE_AUTODIFF)
	endif()

	if(USE_AMGCL)
		add_simplex_headers(ext amgcl)

		add_definitions(-DUSE_AMGCL)
		## use boost
		# add_simplex_headers(ext boost)
		## use no boost
		add_definitions(-DAMGCL_NO_BOOST)
	endif()

	if(USE_PERLIN_NOISE)
		add_headers_and_src(ext perlin_noise)
		add_definitions(-DUSE_PERLIN_NOISE)
	endif()

endmacro(add_headers_and_link_libs)

####proj_names and proj_paths are set in the CMakelists.txt in simplex/projects/_all
macro(add_and_link_multiple_projects)
	set_compiling_flags()
	add_headers_and_link_libs()

	list(LENGTH proj_names proj_n)
	math(EXPR n "${proj_n}-1")
	if(proj_n GREATER 0)
		foreach(i RANGE ${n})
			list(GET proj_names ${i} pname)
			list(GET proj_paths ${i} ppath)
			file(GLOB_RECURSE ${pname}_cpp ${ppath}/*.cpp)
			file(GLOB_RECURSE ${pname}_h ${ppath}/*.h)
			file(GLOB_RECURSE ${pname}_hpp ${ppath}/*.hpp)
			if(USE_CUDA)
				file(GLOB_RECURSE ${pname}_cu ${ppath}/*.cu)
				file(GLOB_RECURSE ${pname}_cuh ${ppath}/*.cuh)
			elseif(USE_CUDA)
				set(${pname}_cu "")
				set(${pname}_cuh "")
			endif(USE_CUDA)

			source_group("xproj/${pname}" FILES ${${pname}_cpp} ${${ppath}_h} ${${ppath}_cu} ${${ppath}_cuh})
			set(exe_src_files ${${pname}_cpp} ${${ppath}_h} ${${ppath}_hpp} ${${ppath}_cu} ${${ppath}_cuh})

			if(UNIX)	#TODO: fix the problem of redundent compiling
				list(APPEND exe_src_files ${src_files})
			endif(UNIX)

			if(USE_CUDA)
				#cuda_add_executable(${pname} ${exe_src_files} OPTIONS ${cuda_options})
				add_executable(${pname} ${exe_src_files})
			else(USE_CUDA)
				add_executable(${pname} ${exe_src_files})
			endif(USE_CUDA)

			set_target_properties(${pname} PROPERTIES PROJECT_LABEL xproj-${pname})
			message(STATUS "Add exe: " ${pname})

			list(LENGTH lib_files lib_n)
			if(lib_n GREATER 0)
				target_link_libraries(${pname} ${lib_files})
			endif()

			list(LENGTH dep_files dep_n)
			if(dep_n GREATER 0)
				add_dependencies(${pname} ${dep_files})
			endif()
		endforeach()
	endif()
endmacro(add_and_link_multiple_projects)

macro(link_python_project)
	set_compiling_flags()
	add_headers_and_link_libs()

	find_package(PythonLibs REQUIRED)
	# add_subdirectory(${PYTHON_INCLUDE_DIRS})
	add_executable(${proj_name} ${src_files})
	include_directories(${PYTHON_INCLUDE_DIRS})
	message(STATUS "Include python file: " ${PYTHON_INCLUDE_DIRS})
	list(APPEND lib_files ${PYTHON_LIBRARIES})

	set_target_properties(${proj_name} PROPERTIES PROJECT_LABEL xproj-${proj_name})
	message(STATUS "Add python module: " ${proj_name})

	list(LENGTH lib_files lib_n)
	if(lib_n GREATER 0)
		target_link_libraries(${proj_name} ${lib_files})
	endif()

	list(LENGTH dep_files dep_n)
	if(dep_n GREATER 0)
		add_dependencies(${proj_name} ${dep_files})
	endif()

	debug_info()
endmacro(link_python_project)


#set exe path to proj source dir, for pybind11 modulue
macro(set_exe_output_path path)
	set(exe_output_path ${path})
	set(CMAKE_RUNTIME_OUTPUT_DIRECTORY "${exe_output_path}")	#for executable files
	set(CMAKE_LIBRARY_OUTPUT_DIRECTORY "${exe_output_path}")	#for dll files
	message(STATUS "set exe_output_path=" ${exe_output_path})
endmacro(set_exe_output_path)

macro(set_dataset_path path)
	if(NOT dataset_path)
		set(dataset_path ${path})
		add_definitions(-DDATASET_PATH=\"${dataset_path}\")
	endif(NOT dataset_path)
endmacro(set_dataset_path)

macro(pass_flags_to_src)
	if(USE_CUDA)
		add_definitions(-DUSE_CUDA)
	endif()
endmacro(pass_flags_to_src)

function(set_list Listname valname)
  #Note: assume you call set_list(My_List my_val)
  #Listname is a string "Listname"
  #${Listname} is a string "My_List"
  #${${Listname}} is the content of My_List
  foreach(x ${${Listname}})
    #The PARENT_SCOPE cannot be omitted here, otherwise the scope will be only this function.
    set(${x} ${valname} PARENT_SCOPE)
  endforeach()
endfunction()

function(Print_List Listname)
  message(STATUS "List info of ${Listname} :")
  foreach(x ${${Listname}})
    message(STATUS "${x} = ${${x}}")
  endforeach()
endfunction()

#if ${lib_flag} is defined, add simplex lib src/${lib_name} (like: src/common) and link to project.
macro(add_simplex_src_lib_on lib_flag lib_name)
  if(${lib_flag})
    message(STATUS "add simplex lib: src/${lib_name}")
    if(WIN32)		#Windows supports linking multiple libs
        compile_and_link_lib(src ${lib_name} ${simplex_lib_path})
    elseif(UNIX)	#but Linux doesn't, so src files have to be compiled
        add_headers_and_src(src ${lib_name})
    endif(WIN32)
  else()
    message(STATUS "ignore simplex lib: src/${lib_name}")
  endif()
endmacro(add_simplex_src_lib_on)

#update ${lib_files} to link a lib.
macro(add_lib_to_link lib_path lib_name)
	if(WIN32)
		set(libs_to_link debug ${lib_path}/${lib_name}d.lib optimized ${lib_path}/${lib_name}.lib)
	elseif(UNIX)
		set(libs_to_link debug ${lib_path}/lib${lib_name}d.a optimized ${lib_path}/lib${lib_name}.a)
	endif(WIN32)
	list(APPEND lib_files ${libs_to_link})
endmacro(add_lib_to_link)

#If ${var_name} is not defined, set it to ${var_value}
function(new_var var_name var_value)
  if(NOT ${var_name})
    set(${var_name} ${var_value} PARENT_SCOPE)
    message(STATUS "set variable: ${var_name} = ${var_value}")
  else()
    message(STATUS "variable already set: ${var_name} = ${${var_name}}")
  endif()
endfunction()

# macro(set_simplex_root path)
#   new_var_with_def(simplex_path -DROOT_PATH ${path})
# endmacro(set_simplex_root)
#
# macro(set_project_root path)
#   new_var(prjplex_path ${path})
# endmacro(set_project_root)

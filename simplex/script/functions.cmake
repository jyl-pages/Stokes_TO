###############################################################################
####Variables in cmake
###############################################################################
####Names:
####proj_name: set in each project's CMakeLists.txt. The name of project.
####current_target: the name of current compiling target. It can be a project(like fluid_euler), a src lib(like common) or an ext lib(like tri2d)
###############################################################################
####Global configs:
####src_lib_lis: the list of all src lib_lis
####ext_header_lib_lis: the list of all ext libs that only contain headers.
####ext_uncompiled_lib_lis:
###############################################################################
####Persistent compilation configs:
####${current_target}_srcs: source files of current target
####${current_target}_incs: the directories that will be included by it.
####${current_target}_libs: .lib files that will be linked to current target.
####${current_target}_deps: the dependent targets
####${current_target}_dlls: .dll files that needed to be copied along with the .exe
###############################################################################
####Temporal compilation configs:
####target_src_libs: the src libs that are used for the current target.
####target_ext_libs: the ext libs that are used for the current target.
####ENABLE_COMIPLE: ON/OFF. When OFF we will not try to add compile targets by add_subdirectory
####ENABLE_LINK: ON/OFF. When OFF we will not try to link anything.
###############################################################################
####Paths:
####simplex_path=simplex
####prjplex_path=simplex or complex, depending on the project
####data_path=simplex/data
####script_path=simplex/script
####exe_output_path=simplex/bin/${platform}/${proj_name}
####src_lib_path: depends on the compilation mode. Can be under bin/{proj_name}, or simplex/lib/simplex.
####ext_lib_path: simplex/lib/ext.
###############################################################################

#It will include the common_func.cmake in the same directory of this file
cmake_policy(SET CMP0057 NEW)
include(${CMAKE_CURRENT_LIST_DIR}/common_func.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/src_config.cmake)
include(${CMAKE_CURRENT_LIST_DIR}/ext_config.cmake)

#Will be called only once at the beginning of the project.
macro(init_project _simplex_path)
  init_target(${proj_name})
  set_compiling_flags()
	#set platform
	if(WIN32)
		set(platform win)
	elseif(UNIX)
		set(platform unix)
	endif(WIN32)
  #set paths
  new_var_with_def(simplex_path -DROOT_PATH ${_simplex_path})
  set(prjplex_path ${CMAKE_CURRENT_SOURCE_DIR}/../..)
  set(proj_src_path ${PROJECT_SOURCE_DIR})
  add_definitions(-DPROJ_PATH=\"${proj_src_path}\")
  new_var_with_def(data_path -DDATA_PATH ${simplex_path}/data)
  new_var_with_def(script_path -DSCRIPT_PATH ${simplex_path}/script)
  set(exe_output_path ${prjplex_path}/bin/${platform}/${proj_name})
  set(ext_lib_path ${simplex_path}/lib/ext)
  set(src_lib_path ${simplex_path}/lib/simplex)#That will be seen in ext()
endmacro(init_project)

#mode is one of EXE/LIB
macro(init_target _target_name)
  set(current_target ${_target_name})
  if(TARGET ${current_target})
    message(STATUS "target ${current_target} already defined, return")
    return()
  endif()
  message(STATUS "Define target: ${current_target}")
  set(${current_target}_srcs "")
  set(${current_target}_incs ${CMAKE_CURRENT_SOURCE_DIR})
  set(${current_target}_libs "")
  set(${current_target}_deps "")
  set(${current_target}_dlls "")
  set(${current_target}_incs_global ${${current_target}_incs} CACHE INTERNAL "Include directories for ${current_target}")
  set(${current_target}_libs_global ${${current_target}_libs} CACHE INTERNAL "Linked libs for ${current_target}")
  set(${current_target}_dlls_global ${${current_target}_dlls} CACHE INTERNAL "Relied dlls for ${current_target}")
  target_add_srcs("${current_target}" ${CMAKE_CURRENT_SOURCE_DIR})
endmacro(init_target)

####set_compiling_flags needs to be called before any target (lib or exe) is generated
macro(set_compiling_flags)
	set(CMAKE_CXX_STANDARD 17)	#c++17
  message(STATUS "Use c++17")
	if(UNIX)
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wall -std=c++17 -O3")	#c++17
		SET(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Wno-sign-compare")	#turn off sign-compare warning
	endif(UNIX)
	if(WIN32)
		add_definitions(-D_DISABLE_EXTENDED_ALIGNED_STORAGE)	#fix compiling issue for VS2017
    add_definitions(-D_SILENCE_ALL_CXX17_DEPRECATION_WARNINGS)
	endif(WIN32)
endmacro(set_compiling_flags)

macro(ignore_src_modules ign_list)
  set(ign_list ${ign_list})

  #if("parallel" IN_LIST ign_list)
    #set(USE_CUDA OFF)
    #remove_definitions(-DUSE_CUDA)
  #endif()

  if("cpx" IN_LIST ign_list)
    message(STATUS "Disable CPX solver")
    set(USE_CPX OFF)
    remove_definitions(-DUSE_CPX)
  endif()

endmacro(ignore_src_modules)

macro(setup_exe_target)
  add_exe_target(${exe_output_path})
  fast_debug_info()
endmacro(setup_exe_target)

macro(fast_debug_info)
  message(STATUS "Build info for target ${current_target}:")
  message(STATUS "Src files: ${${current_target}_srcs}")
  message(STATUS "Included dirs: ${${current_target}_incs}")
  message(STATUS "Lib files: ${${current_target}_libs}")
  message(STATUS "Dep targets: ${${current_target}_deps}")
  message(STATUS "Dll files: ${${current_target}_dlls}")
endmacro(fast_debug_info)

macro(debug_info)
	message(STATUS "Information for inc,src, lib, and dependencies:")

  list(LENGTH ${current_target}_incs inc_n)
	message(STATUS "inc_n: " ${inc_n})
	math(EXPR n "${inc_n}-1")
	if(inc_n GREATER 0)
		foreach(i RANGE ${n})
			list(GET ${current_target}_incs ${i} f)
			message(STATUS "Inc path " ${i} ": " ${f})
		endforeach()
	endif()

  list(LENGTH ${current_target}_srcs src_n)
	message(STATUS "src_n: " ${src_n})
	math(EXPR n "${src_n}-1")
	if(src_n GREATER 0)
		foreach(i RANGE ${n})
			list(GET ${current_target}_srcs ${i} f)
			message(STATUS "Src file " ${i} ": " ${f})
		endforeach()
	endif()

	list(LENGTH ${current_target}_libs lib_n)
	message(STATUS "lib_n: " ${lib_n})
	math(EXPR n "${lib_n}-1")
	if(lib_n GREATER 0)
		foreach(i RANGE ${n})
			list(GET ${current_target}_libs ${i} f)
			message(STATUS "Lib file " ${i} ": " ${f})
		endforeach()
	endif()

	list(LENGTH ${current_target}_deps dep_n)
	message(STATUS "dep_n: " ${dep_n})
	math(EXPR n "${dep_n}-1")
	if(dep_n GREATER 0)
		foreach(i RANGE ${n})
			list(GET ${current_target}_deps ${i} f)
			message(STATUS "Dep file " ${i} ": " ${f})
		endforeach()
	endif()

endmacro(debug_info)

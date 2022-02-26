#Notes about cmake syntax.
#Suppose you write a function foo(var)
#var is just a string "var"
#${var} is the value of var. Like, if it's integer 1, then ${var}=1
#Further, if ${var} is another thing, say, a list, then ${${var}} is the exact contents of that list.

#If ${var_name} is not defined, set it to ${var_value}, and also set a compilation flag ${def_name} to it.
function(new_var_with_def var_name def_name var_value)
  if(NOT ${var_name})
    set(${var_name} ${var_value} PARENT_SCOPE)
    add_definitions(\"${def_name}\"=\"${var_value}\")
    message(STATUS "set variable with definition: ${var_name} = ${def_name} = ${var_value}")
  else()
    message(STATUS "variable and definition already set: ${var_name} = ${def_name} = ${${var_name}}")
  endif()
endfunction()

#save all code files under ${file_path} to ${all_files}
function(all_code_files_under file_path all_files)
  file(GLOB_RECURSE cpp_files ${file_path}/*.cpp)
	file(GLOB_RECURSE hpp_files ${file_path}/*.hpp)
	file(GLOB_RECURSE h_files ${file_path}/*.h)
  file(GLOB_RECURSE cu_files ${file_path}/*.cu)
	file(GLOB_RECURSE cuh_files ${file_path}/*.cuh)
  set(${all_files} ${cpp_files} ${hpp_files} ${h_files} ${cu_files} ${cuh_files} PARENT_SCOPE)
endfunction(all_code_files_under)

macro(add_exe_target output_path)
  message(STATUS "Compile exe ${current_target} to ${output_path}.")
  set(CMAKE_RUNTIME_OUTPUT_DIRECTORY ${output_path})
  set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${output_path})
  add_executable(${current_target} ${${current_target}_srcs})
  target_setup_all_resources(EXE)
  set_target_properties(${current_target} PROPERTIES PROJECT_LABEL xproj-${current_target})
  if(CMAKE_CUDA_COMPILER)#if cuda is presented
    set_property(TARGET ${current_target} PROPERTY CUDA_ARCHITECTURES 61)
    set_property(TARGET ${current_target} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${current_target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    						--extended-lambda
    						>)
  endif()
  copy_dll_files()
endmacro(add_exe_target)

macro(copy_dll_files)
  list(LENGTH ${current_target}_dlls dll_n)
  if(dll_n GREATER 0)
    foreach(dll_file ${${current_target}_dlls})
      if(WIN32)
        add_custom_command(TARGET ${current_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${exe_output_path}/Release")
        add_custom_command(TARGET ${current_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E make_directory "${exe_output_path}/Debug")
        add_custom_command ( TARGET ${current_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${dll_file}" "${exe_output_path}/Release/")
        add_custom_command ( TARGET ${current_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different "${dll_file}" "${exe_output_path}/Debug/")
      elseif(UNIX)
        add_custom_command ( TARGET ${current_target} POST_BUILD COMMAND ${CMAKE_COMMAND} -E copy_if_different ${dll_file} ${exe_output_path})
      endif()
    endforeach()
  endif()
endmacro()

#like: compile_lib(src common D:/simplex/src/common D:/complex/bin/win/fluid_euler/lib/simplex)
macro(add_lib_target lib_prefix output_path)
  message(STATUS "Compile lib ${lib_prefix}-${current_target} to ${output_path}.")
  add_library(${current_target} STATIC ${${current_target}_srcs})
  target_setup_all_resources(LIB)
  set_target_properties(${current_target} PROPERTIES DEBUG_POSTFIX d)
  set_target_properties(${current_target} PROPERTIES PROJECT_LABEL ${lib_prefix}-${current_target})
  set_target_properties(${current_target} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY "${output_path}")
  set_target_properties(${current_target} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_DEBUG "${output_path}")
  set_target_properties(${current_target} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELEASE "${output_path}")
  set_target_properties(${current_target} PROPERTIES ARCHIVE_OUTPUT_DIRECTORY_RELWITHDEBINFO  "${output_path}")
  if(CMAKE_CUDA_COMPILER)
    set_property(TARGET ${current_target} PROPERTY CUDA_ARCHITECTURES 61)
    set_property(TARGET ${current_target} PROPERTY CUDA_SEPARABLE_COMPILATION ON)
    target_compile_options(${current_target} PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    						--extended-lambda
    						>)
  endif()
endmacro(add_lib_target)

#Set things in ${current_target}_incs, ${current_target}_libs,
#${current_target}_deps actually to the current target
#target_type is "EXE" or "LIB". Only with "EXE", we will link the .lib files.
macro(target_setup_all_resources target_type)
  foreach(inc_path ${${current_target}_incs})
    target_include_directories(${current_target} PUBLIC ${inc_path})
  endforeach()

  list(LENGTH ${current_target}_libs lib_n)
  if(lib_n GREATER 0)
    target_link_libraries(${current_target} PUBLIC ${${current_target}_libs})
  endif()

  list(LENGTH ${current_target}_deps dep_n)
  if(dep_n GREATER 0)
    add_dependencies(${current_target} ${${current_target}_deps})
  endif()
endmacro()

macro(pass_resources_to_global)
  set(${current_target}_incs_global ${${current_target}_incs} CACHE STRING INTERNAL FORCE)
  set(${current_target}_libs_global ${${current_target}_libs} CACHE STRING INTERNAL FORCE)
  set(${current_target}_dlls_global ${${current_target}_dlls} CACHE STRING INTERNAL FORCE)
endmacro()

macro(merge_resources_from_target dep_name)
  if(WIN32)
    list(APPEND ${current_target}_incs ${${dep_name}_incs_global})
    list(APPEND ${current_target}_libs ${${dep_name}_libs_global})
  else()
    #In linux, the referrer goes before the referred
    list(PREPEND ${current_target}_incs ${${dep_name}_incs_global})
    list(PREPEND ${current_target}_libs ${${dep_name}_libs_global})
  endif()
  list(APPEND ${current_target}_dlls ${${dep_name}_dlls_global})
  list(REMOVE_DUPLICATES ${current_target}_incs)
  list(REMOVE_DUPLICATES ${current_target}_dlls)
  #REMOVE_DUPLICATES does not work well with libs, because there're multiple "debug"s and "optimized"s
endmacro()

#Mark header and source files to add
#Like: target_add_srcs(simplex/src/common)
macro(target_add_srcs group_name src_path)
  all_code_files_under(${src_path} _code_files)
  source_group(${group_name} FILES ${_code_files})
  list(APPEND ${current_target}_srcs ${_code_files})
endmacro(target_add_srcs)

#Mark header files to include
#Like: target_add_incs(simplex/src/common)
macro(target_add_incs header_path)
  list(APPEND ${current_target}_incs ${header_path})
endmacro(target_add_incs)

#Mark libs to link
#Like: target_add_libs(simplex/lib/simplex geometry)
macro(target_add_libs lib_path lib_name)
  if(ENABLE_LINK)
    if(WIN32)
  	   set(_libs_to_link debug ${lib_path}/${lib_name}d.lib optimized ${lib_path}/${lib_name}.lib)
       list(APPEND ${current_target}_libs ${_libs_to_link})
    elseif(UNIX)
  	   set(_libs_to_link debug ${lib_path}/lib${lib_name}d.a optimized ${lib_path}/lib${lib_name}.a)
       list(PREPEND ${current_target}_libs ${_libs_to_link})
    endif()
  endif()
endmacro(target_add_libs)

#Mark targets to depend on
#Like: target_add_deps("common;geometry").
#Or: target_add_deps(common)
macro(target_install_dep _dep_target_path _dep_build_path _dep_name)
  if(ENABLE_COMIPLE)
    if(NOT TARGET ${_dep_name})
      add_subdirectory(${_dep_target_path} ${_dep_build_path})
    endif()
    #If CUDA is not set while we're trying to compile "parallel", then no target will be created.
    #Then we should not add dependencies.
    if(TARGET ${_dep_name})
      list(APPEND ${current_target}_deps ${_dep_name})
      merge_resources_from_target(${_dep_name})
      target_add_incs(${simplex_path}/src/${lib_name})
    endif()
  endif()
endmacro(target_install_dep)

#The order here must be carefully set. If lib B depends on lib A, then A must come before B, otherwise A will not be included when compiling B.
set(src_lib_lis common reservoir common_cpx poisson_cpx geometry optimization parallel solver viewer physics)

macro(setup_src_target)
  add_lib_target(src ${src_lib_path})
  fast_debug_info()
  pass_resources_to_global()
  message("")
endmacro(setup_src_target)

#mode: one of SELF/COMMON/COMPILE/LINK
#SELF/COMMON will set src_lib_path, while COMPILE/LINK won't.
#In SELF/COMPILE mode, cmake will try to compile all these src libs.
#In COMMON/LINK mode, cmake will only link existed .lib files.
macro(use_src_libs mode lib_lis)
  set(lib_lis ${lib_lis})
  if(${mode} STREQUAL "SELF")
    set(src_lib_path ${prjplex_path}/bin/${platform}/${proj_name}/lib/simplex)
    set(ENABLE_COMIPLE ON)
    set(ENABLE_LINK ON)
  elseif(${mode} STREQUAL "COMMON")
    set(src_lib_path ${simplex_path}/lib/simplex)
    set(ENABLE_COMIPLE OFF)
    set(ENABLE_LINK ON)
  elseif(${mode} STREQUAL "COMPILE")
    set(ENABLE_COMIPLE ON)
    set(ENABLE_LINK ON)
  elseif(${mode} STREQUAL "LINK")
    set(ENABLE_COMIPLE OFF)
    set(ENABLE_LINK ON)
  elseif(${mode} STREQUAL "HEADERS")
    set(ENABLE_COMIPLE OFF)
    set(ENABLE_LINK OFF)
  else()
    message(FATAL_ERROR "use_src_libs unrecognized mode ${mode}")
  endif()

  #check if the names are correct
  foreach(lib_name ${lib_list})
    if(NOT ${lib_name} IN_LIST src_lib_lis)
      message(FATAL_ERROR "${lib_name} is not a library under simplex/src.")
    endif()
  endforeach()

  foreach(lib_name ${src_lib_lis})
    if(${lib_name} IN_LIST lib_lis)
        set(target_lib_build_path ${src_lib_path}/${lib_name}/cmake-build-${platform})
        target_install_dep(${simplex_path}/src/${lib_name} ${target_lib_build_path} ${lib_name})
        #That's a special treatment for cuda.
        #If "parallel" is used while USE_CUDA is OFF, it will create no target.
        #Then we should not link any libs.
        if(NOT ENABLE_COMIPLE OR TARGET ${lib_name})
          target_add_libs(${src_lib_path} ${lib_name})
        endif()
    endif()
  endforeach()
endmacro(use_src_libs)

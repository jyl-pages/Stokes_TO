cmake_minimum_required(VERSION 2.8.7)
set(proj_name _all)
project(${proj_name})
set(proj_src_path ${PROJECT_SOURCE_DIR})
set(CMAKE_BINARY_DIR ${PROJECT_SOURCE_DIR}/build)
message(STATUS "src: " ${PROJECT_SOURCE_DIR})
message(STATUS "bin: " ${CMAKE_BINARY_DIR})

set(proj_names "")
list(APPEND proj_names 
helloworld
fluid_euler
poisson
levelset
opengl_viewer
rigid
soft_body_fem
)

set(proj_paths "")
foreach(pname ${proj_names})
	list(APPEND proj_paths ${proj_src_path}/../${pname})
endforeach()

include(${PROJECT_SOURCE_DIR}/../../script/functions.cmake)

init_project()
set_all_options_on()
set(USE_CUDA OFF)
add_and_link_multiple_projects()
debug_info()


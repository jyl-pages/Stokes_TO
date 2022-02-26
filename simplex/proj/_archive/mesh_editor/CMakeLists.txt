cmake_minimum_required(VERSION 2.8.7)
set(proj_name mesh_editor)
project(${proj_name})
include(${PROJECT_SOURCE_DIR}/../../script/functions.cmake)
set_simplex_root(../../../simplex)
set_proj_root(../..)


init_project()
set(USE_GEOMETRY OFF)
set(USE_PHYSICS OFF)
set(USE_VIEWER ON)
set(USE_IMGUI ON)
set(USE_STB ON)
set(USE_TINY_OBJ_LOADER ON)
set(USE_RAPIDJSON OFF)

add_project()
link_project()
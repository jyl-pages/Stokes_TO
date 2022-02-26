macro(link_pybind_project_to_python)
	set_compiling_flags()
	add_headers_and_link_libs()

	add_subdirectory(${ext_path}/pybind11 ${ext_path}/pybind11)
	pybind11_add_module(${proj_name} ${src_files})

	set_target_properties(${proj_name} PROPERTIES PROJECT_LABEL xproj-${proj_name})
	message(STATUS "Add pybind module: " ${proj_name})

	list(LENGTH lib_files lib_n)
	if(lib_n GREATER 0)
		target_link_libraries(${proj_name} ${lib_files})
	endif()

	list(LENGTH dep_files dep_n)
	if(dep_n GREATER 0)
		add_dependencies(${proj_name} ${dep_files})
	endif()

	debug_info()
endmacro(link_pybind_project_to_python)

macro(link_pybind_project_from_python)
	set_compiling_flags()
	add_headers_and_link_libs()

	set(PYBIND11_PYTHON_VERSION 3.6.5 CACHE STRING "")
	add_subdirectory(${ext_path}/pybind11 ${ext_path}/pybind11)
	add_executable(${proj_name} ${src_files})

	set_target_properties(${proj_name} PROPERTIES PROJECT_LABEL xproj-${proj_name})
	message(STATUS "Add pybind exe: " ${proj_name})

	target_link_libraries(${proj_name} PRIVATE pybind11::embed)
	list(LENGTH lib_files lib_n)
	if(lib_n GREATER 0)
		target_link_libraries(${proj_name} ${lib_files})
	endif()

	list(LENGTH dep_files dep_n)
	if(dep_n GREATER 0)
		add_dependencies(${proj_name} ${dep_files})
	endif()

	debug_info()
endmacro(link_pybind_project_from_python)

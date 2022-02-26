include(${CMAKE_CURRENT_LIST_DIR}/common_func.cmake)

macro(set_default_options)
	#simplex files
	#ON: USE_COMMON, USE_GEOMETRY, USE_SOLVER, USE_PHYSICS
	#OFF: USE_OPTIMIZATION, USE_VIEWER
	set(ListOptionSource
		USE_COMMON # ON
		USE_GEOMETRY # ON
		USE_SOLVER # ON
		USE_OPTIMIZATION # OFF
		USE_PHYSICS # ON
		USE_VIEWER # OFF
	)
	set_list(ListOptionSource ON)
	set(USE_OPTIMIZATION OFF)
	set(USE_VIEWER OFF)

	#simplex files compiled
	#All default OFF
	set(ListOptionCompiled
		USE_COMMON_COMPILED
		USE_GEOMETRY_COMPILED
		USE_SOLVER_COMPILED
		USE_OPTIMIZATION_COMPILED
		USE_PHYSICS_COMPILED
		USE_VIEWER_COMPILED
	)
	set_list(ListOptionCompiled OFF)

	##Third Parties
	#Part 1: ON with set_all_options_on()
	set(ListOptionThirdCommon
		USE_PHYSBAM # ON
		USE_IMGUI # OFF
		USE_STB # OFF
		USE_TINY_OBJ_LOADER # OFF
		USE_OPENMP # ON
		USE_NANOFLANN # ON for NeighborSearcher.h
		USE_TRI2D # ON
		USE_TETGEN # ON
		USE_AUTODIFF # ON
		USE_CUDA # ON
		USE_AMGCL # ON
		USE_PERLIN_NOISE # ON
	)
	set_list(ListOptionThirdCommon ON)
	set(USE_IMGUI OFF)
	set(USE_STB OFF)
	set(USE_TINY_OBJ_LOADER OFF)

	#Part 2: OFF with set_all_options_on()
	set(ListOptionThirdCustom
		USE_RAPIDJSON # OFF
		USE_IPOPT # OFF
		USE_NLOPT # OFF
		USE_MMA # OFF
		USE_LIBTORCH # OFF
	)
	set_list(ListOptionThirdCustom OFF)

	set(ListOption ${ListOptionSource} ${ListOptionCompiled} ${ListOptionThirdCommon} ${ListOptionThirdCustom})
endmacro(set_default_options)

macro(set_all_options_on)
	set_list(ListOptionSource ON)
	set_list(ListOptionCompiled OFF)
	set_list(ListOptionThirdCommon ON)
	set_list(ListOptionThirdCustom OFF)
endmacro(set_all_options_on)

macro(set_all_options_off)
	set_list(ListOptionSource OFF)
	set_list(ListOptionCompiled OFF)
	set_list(ListOptionThirdCommon OFF)
	set_list(ListOptionThirdCustom OFF)
endmacro(set_all_options_off)

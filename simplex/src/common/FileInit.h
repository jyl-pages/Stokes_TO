//////////////////////////////////////////////////////////////////////////
// Initialize parameters from a file
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FileInit_h__
#define __FileInit_h__
#include "Common.h"
#include "AuxFunc.h"
#include <iostream>
#include <map>

class SimulatorParams {
public:
	std::string name = "";
	SimulatorParams(std::string _name = "") { name = _name; }
	virtual void Read(std::istream& input) {
		std::cerr << "You didn't implement Read() function in some derived class, so it end up here. Please check the code.";
		exit(0);
	}
};

namespace FileInit {

	void Register_Parameter_Lines_From_File(const std::string& file_name, std::map<std::string, std::string>& param_map);
	void Read_Params_With_Map(const std::map<std::string, std::string> param_map, SimulatorParams* param);
	template<typename...Args>
	void Read_Params_With_Map(const std::map<std::string, std::string> param_map, SimulatorParams* param, Args ...rest) {
		Read_Params_With_Map(param_map, param);
		Read_Params_With_Map(param_map, rest...);
	}

	template<typename...Args>
	void Read_Params_From_File(const std::string& file_name, Args ...rest) {
		std::map<std::string, std::string> param_map;
		Register_Parameter_Lines_From_File(file_name, param_map);
		Read_Params_With_Map(param_map, rest...);
	}
};

#endif

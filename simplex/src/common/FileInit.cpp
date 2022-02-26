//////////////////////////////////////////////////////////////////////////
// Functions for FileInit.h. There're some string operation functions here.
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "FileInit.h"
#include <fstream>


void Squash_File_Content(const std::string& file_name, Array<std::string> &tokens) {
	//Ignore all contents commented out by '#'. And split with white spaces
	std::ifstream fin(file_name);
	if (!fin.is_open()) {
		std::cerr << "Squash_File_Content error: file not exitsts: " << file_name << "\n";
		exit(0);
	}
	std::string str;
	tokens.clear();
	while (std::getline(fin, str)) {
		size_t k = str.find('#');
		if (k != std::string::npos) {
			str = str.substr(0, k);
		}
		Array<std::string> line_tokens = AuxFunc::Split_String(str);
		tokens.insert(tokens.end(), line_tokens.begin(), line_tokens.end());
	}
	fin.close();
}

void FileInit::Register_Parameter_Lines_From_File(const std::string& file_name, std::map<std::string, std::string> &param_map) {
	Array<std::string> tokens;
	Squash_File_Content(file_name, tokens);
	param_map.clear();
	std::string now_param_name = "";
	std::string now_param_str = "";
	tokens.push_back("@");//make sure to register the last "@" token
	for (int i = 0; i < tokens.size(); i++) {
		std::string s = tokens[i];
		if (s[0] == '@') {
			if (now_param_name == "") {
				if (now_param_str != "") {
					std::cerr << "FileInit::Register_Parameter_Lines_From_File Error: cannot register a parameter with empty name\"" << now_param_str << "\"\n";
					exit(1);
				}
			}
			else {
				if (param_map.count(now_param_name)) {
					std::cerr << "FileInit::Register_Parameter_Lines_From_File Error: param name \"" << now_param_name << "\" already registered\n";
					exit(1);
				}
				param_map[now_param_name] = now_param_str;
			}
			now_param_name = s.substr(1);//filter out that '@'
			now_param_str = "";
		}
		else {
			if (now_param_str != "") now_param_str += " ";
			now_param_str += s;
		}
	}
}

void FileInit::Read_Params_With_Map(const std::map<std::string, std::string> param_map, SimulatorParams* param){
	auto iter = param_map.find(param->name);
	if (iter == param_map.end()) {
		std::cerr << "FileInit::Read_Params_With_Map error: cannot find param with name=" << param->name;
		exit(1);
	}
	std::istringstream sin(iter->second);
	param->Read(sin);
	if (sin.fail()) {
		std::cerr << "FileInit::Read_Params_With_Map error: " << param->name << " is not fully initialized, missing parameters\n";
		exit(1);
	}
	if (!sin.eof()) {
		std::cerr << "FileInit::Read_Params_With_Map error: params for " << param->name << " are not fully read, there are redundant parameters: " << sin.rdbuf() << "\n";
		exit(1);
	}
}


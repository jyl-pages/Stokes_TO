//#####################################################################
// Opengl viewer
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include "ParseArgs.h"
#include "OpenGLMeshEditor.h"

int main(int argc,char* argv[])
{
    ParseArgs parse_args;
    parse_args.Add_String_Argument("-o","","data path");
    parse_args.Parse(argc,argv);

	std::shared_ptr<OpenGLMeshEditor> viewer;
	viewer.reset(new OpenGLMeshEditor());
	
	std::string file_name=parse_args.Get_String_Value("-o");
	if(file_name!="")viewer->file_name=file_name;
	
	viewer->Initialize();
	viewer->Run();
	viewer->Finish();

	return 0;
}
//#####################################################################
// Copyright 2004-2007, Eran Guendelman, Geoffrey Irving, Andrew Selle.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class PARSE_ARGS  
//##################################################################### 
#ifndef __PARSE_ARGS__
#define __PARSE_ARGS__
#include <cstring>
#include <iostream>
#include "PhysBAMAdaptor.h"
#include "ArgData.h"

class PARSE_ARGS
{
private:
    std::vector<ARG_DATA> arg_data_list;
    std::vector<std::string> extra_arg_list,extra_args_synopsis_list,extra_args_desc_list;
    int num_expected_extra_args;
    std::string program_name;
    bool use_help_option;
    void (*extra_usage_callback)();
public:

    PARSE_ARGS():num_expected_extra_args(0),use_help_option(true),extra_usage_callback(0){}

    void Use_Help_Option(bool use_it)
	{
		use_help_option=use_it;
	}

	void Add_Option_Argument(const std::string& arg_str,const std::string& desc="")
	{
		arg_data_list.push_back(ARG_DATA(arg_str,desc));
	}

	void Add_Integer_Argument(const std::string& arg_str,int default_value,const std::string& val_name="",const std::string& desc="")
	{
		arg_data_list.push_back(ARG_DATA(arg_str,val_name,desc,default_value));
	}

	void Add_Double_Argument(const std::string& arg_str,double default_value,const std::string& val_name="",const std::string& desc="")
	{
		arg_data_list.push_back(ARG_DATA(arg_str,val_name,desc,default_value));
	}

	void Add_Vector_2D_Argument(const std::string& arg_str,const Vector2 &default_value,const std::string& val_name="",const std::string& desc="")
	{
		arg_data_list.push_back(ARG_DATA(arg_str,val_name,desc,default_value));
	}

	void Add_Vector_3D_Argument(const std::string& arg_str,const Vector3 &default_value,const std::string& val_name="",const std::string& desc="")
	{
		arg_data_list.push_back(ARG_DATA(arg_str,val_name,desc,default_value));
	}

	void Add_String_Argument(const std::string& arg_str,const std::string& default_value,const std::string& val_name="",const std::string& desc="")
	{
		arg_data_list.push_back(ARG_DATA(arg_str,val_name,desc,default_value));
	}

	void Set_Extra_Arguments(int num, const std::string& synopsis="",const std::string& desc="") // num=-1 for arbitrary extra arguments
	{
		num_expected_extra_args=(num!=-1)?num_expected_extra_args+1:num;if(synopsis.length())extra_args_synopsis_list.push_back(synopsis);if(desc.length())extra_args_desc_list.push_back(desc);
	}

	void Set_Extra_Usage_Callback(void (*extra_usage_callback_input)())
	{
		extra_usage_callback=extra_usage_callback_input;
	}

	void Parse(int argc,char* argv[])
	{
		program_name=argv[0];int current_arg=1;
		extra_arg_list.clear();
		while(current_arg<argc){
			if(use_help_option && !strcmp(argv[current_arg],"--help"))Print_Usage(true); // print help
			int match=Find_Match(argv[current_arg]);
			if(match==-1){
				if(argv[current_arg][0]=='-') Print_Usage(true);
				else extra_arg_list.push_back(argv[current_arg++]);}
			else if(!arg_data_list[match].Parse_Value(argc,argv,current_arg))
				Print_Usage(true);
			else arg_data_list[match].value_set=true;}
		if(num_expected_extra_args!=-1 && (int)extra_arg_list.size()!=num_expected_extra_args) Print_Usage(true); // didn't get the expected number of extra args
		for(int i=current_arg;i<argc;i++)extra_arg_list.push_back(argv[i]);
	}

	void Parse(int argc,const char* argv[])
	{
		Parse(argc,(char**)argv);
	}

	bool Get_Option_Value(const std::string& arg_str) const
	{
		return arg_data_list[Find_Match(arg_str,ARG_DATA::OPTION)].boolean_value;
	}

	int Get_Integer_Value(const std::string& arg_str) const
	{
		return arg_data_list[Find_Match(arg_str,ARG_DATA::INTEGER)].integer_value;
	}

	double Get_Double_Value(const std::string& arg_str) const
	{
		return arg_data_list[Find_Match(arg_str,ARG_DATA::DOUBLE)].double_value;
	}

	Vector2 Get_Vector_2D_Value(const std::string& arg_str) const
	{
		return arg_data_list[Find_Match(arg_str,ARG_DATA::VECTOR2)].vector_2d_value;
	}

	Vector3 Get_Vector_3D_Value(const std::string& arg_str) const
	{
		return arg_data_list[Find_Match(arg_str,ARG_DATA::VECTOR3)].vector_3d_value;
	}

	const std::string& Get_String_Value(const std::string& arg_str) const
	{
		return arg_data_list[Find_Match(arg_str,ARG_DATA::STRING)].string_value;
	}

	bool Is_Value_Set(const std::string& arg_str) const
	{
		int match=Find_Match(arg_str);
		if(match==-1){std::stringstream ss;ss<<"Argument "<<arg_str<<" undeclared"<<std::endl;}
		return arg_data_list[match].value_set;
	}

	void Override_String_Value(const std::string& arg_str,const std::string& value)
	{
		arg_data_list[Find_Match(arg_str,ARG_DATA::STRING)].string_value=value;
	}

	int Find_Match(const std::string& str) const
	{
		for(int i=0;i<(int)arg_data_list.size();i++)if(arg_data_list[i].str==str)return i;
		return -1;
	}

	int Find_Match(const std::string& str,const ARG_DATA::TYPE& type) const
	{
		int match=Find_Match(str);
		if(match==-1){std::cerr<<"Argument "<<str<<" undeclared"<<std::endl;}
		if(arg_data_list[match].type!=type){std::cerr<<"Type mismatch in Find_Match("<<str<<")"<<std::endl;}
		return match;
	}

	int Num_Extra_Args() const
	{
		return (int)extra_arg_list.size();
	}

	const std::string& Extra_Arg(int i) const
	{
		assert(0<=i && i<(int)extra_arg_list.size());
		return extra_arg_list[i];
	}

	const std::string& Get_Program_Name() const
	{
		return program_name;
	}

	bool Find_And_Remove(const char *str,int& argc,char** argv)
	{
		int i;for(i=0;i<argc;i++)if(!strcmp(str,argv[i]))break;
		if(i<argc){for(;i<argc-1;i++)argv[i]=argv[i+1];argc--;argv[argc]=0;return true;}
		return false;
	}

	double Find_And_Remove_Double(const char *str,int& argc,char** argv)
	{
		double value;
		int i;for(i=0;i<argc;i++)if(!strcmp(str,argv[i]))break;
		if(i+1<argc){value=atof(argv[i+1]);for(;i<argc-2;i++)argv[i]=argv[i+2];argc--;argv[argc]=0;argc--;argv[argc]=0;return value;}
		return 0.0;
	}

	int Find_And_Remove_Integer(const char *str,int& argc,char** argv)
	{
		int value;
		int i;for(i=0;i<argc;i++)if(!strcmp(str,argv[i]))break;
		if(i+1<argc){value=atoi(argv[i+1]);for(;i<argc-2;i++)argv[i]=argv[i+2];argc--;argv[argc]=0;argc--;argv[argc]=0;return value;}
		return 0;
	}

	void Print_Usage(bool do_exit) const
	{
		int i;std::cout<<"Usage: "<<program_name<<" ";
		for(i=0;(int)i<arg_data_list.size();i++){arg_data_list[i].Print_Synopsis();std::cout<<" ";}
		for(i=0;i<(int)extra_args_synopsis_list.size();i++) std::cout<<extra_args_synopsis_list[i]<<" ";std::cout<<std::endl;
		int width=0,len;
		for(i=0;i<(int)arg_data_list.size();i++){len=(int)arg_data_list[i].str.length();if(len>width)width=len;}
		for(i=0;i<(int)extra_args_desc_list.size();i++){len=(int)extra_args_synopsis_list[i].length();if(len>width)width=len;}
		for(i=0;i<(int)arg_data_list.size();i++){arg_data_list[i].Print_Description(width+2);std::cout<<std::endl;}
		for(i=0;i<(int)extra_args_desc_list.size();i++){
			std::cout<<extra_args_synopsis_list[i];for(unsigned j=1;j<=width+2-extra_args_synopsis_list[i].length();j++) std::cout<<" ";std::cout<<extra_args_desc_list[i]<<std::endl;} 
		if(extra_usage_callback)extra_usage_callback();if(do_exit) exit(-1);
	}
};
using ParseArgs=PARSE_ARGS;
#endif

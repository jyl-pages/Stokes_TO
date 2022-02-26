//#####################################################################
// Rapid Json
// Interface functions for SLAX
//#####################################################################
#ifndef __Rapid_Json_Func_h__
#define __Rapid_Json_Func_h__
#include <iostream>
#include <string>
#include "document.h"
#include "writer.h"
#include "stringbuffer.h"

using JsonDoc=rapidjson::Document;
using JsonValue=rapidjson::Value;
using JsonObj=rapidjson::Value::Object;

inline void Parse_Json_From_String(const std::string& str,JsonDoc& doc)
{
	using namespace rapidjson;
	doc.Parse(str.c_str());
}

inline void Parse_Json_From_File(const std::string& file_name,JsonDoc& doc)
{
	std::ifstream f(file_name);std::stringstream ss;ss<<f.rdbuf();
	std::cout<<"Json config file: "<<file_name<<"\n"<<ss.str()<<std::endl;
	Parse_Json_From_String(ss.str(),doc);
}

inline std::string Write_Json_To_String(const JsonDoc& doc)
{
	using namespace rapidjson;
	StringBuffer buffer;
    Writer<StringBuffer> writer(buffer);
    doc.Accept(writer);	
	return std::string(buffer.GetString());
}

inline bool Has_String_Value(const JsonObj& object,const std::string name)
{
	return object.HasMember(name.c_str())&&object.FindMember(name.c_str())->value.IsString();
}

inline std::string Get_String_Value(const JsonObj& object,const std::string name)
{
	auto obj=object.FindMember(name.c_str());
	std::string value=obj->value.GetString();
	return value;
}

inline bool Get_If_Has_String_Value(const JsonObj& object,const std::string name,std::string& value)
{
	if(!Has_String_Value(object,name))return false;
	auto obj=object.FindMember(name.c_str());
	value=obj->value.GetString();return true;
}

inline bool Has_Bool_Value(const JsonObj& object,const std::string name)
{
	return object.HasMember(name.c_str())&&object.FindMember(name.c_str())->value.IsBool();
}

inline bool Get_Bool_Value(const JsonObj& object,const std::string name)
{
	auto obj=object.FindMember(name.c_str());
	bool value=obj->value.GetBool();
	return value;
}

inline bool Has_Float_Value(const JsonObj& object,const std::string name)
{
	return object.HasMember(name.c_str())&&object.FindMember(name.c_str())->value.IsFloat();
}

inline float Get_Float_Value(const JsonObj& object,const std::string name)
{
	auto obj=object.FindMember(name.c_str());
	float value=obj->value.GetFloat();
	return value;
}

inline bool Has_Double_Value(const JsonObj& object,const std::string name)
{
	return object.HasMember(name.c_str())&&object.FindMember(name.c_str())->value.IsDouble();
}

inline double Get_Double_Value(const JsonObj& object,const std::string name)
{
	auto obj=object.FindMember(name.c_str());
	double value=obj->value.GetDouble();
	return value;
}

inline bool Has_Vector4f_Value(const JsonObj& object,const std::string name)
{
	if(!object.HasMember(name.c_str()))return false;
	auto& val=object.FindMember(name.c_str())->value;
	if(!val.IsArray()||val.Size()<4)return false;
	return true;
}

inline Vector4f Get_Vector4f_Value(const JsonObj& object,const std::string name)
{
	auto& val=object.FindMember(name.c_str())->value;
	Vector4f v;for(int i=0;i<(int)val.Size();i++){
		v[i]=val[i].GetFloat();}return v;
}

#endif
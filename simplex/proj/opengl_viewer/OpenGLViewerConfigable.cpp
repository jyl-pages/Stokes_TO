//////////////////////////////////////////////////////////////////////////
// Opengl viewer configable
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "OpenGLViewerConfigable.h"
#include "OpenGLGrid.h"
#include "OpenGLGridField.h"
#include "OpenGLGridVectors.h"
#include "OpenGLGridTensors.h"
#include "OpenGLMesh.h"
#include "OpenGLMeshField.h"
#include "OpenGLMeshVectors.h"
#include "OpenGLMeshTensors.h"
#include "OpenGLParticles.h"
#include "OpenGLScreenObjects.h"
#include "OpenGLViewerTestCase.h"

#ifdef USE_RAPIDJSON
void OpenGLViewerConfiguable::Initialize_Data()
{
	using namespace rapidjson;
	JsonDoc doc;Parse_Json_From_File(config_file_name,doc);
	std::string str2=Write_Json_To_String(doc);

	auto root=doc.FindMember("objects");
	std::cout<<"Parse OpenGL json file "<<config_file_name<<std::endl;
	for(auto iter=root->value.Begin();iter!=root->value.End();++iter){
		const Value::Object& obj=iter->GetObject();
		std::cout<<"object properties"<<std::endl;
		for(auto iter2=obj.MemberBegin();iter2!=obj.MemberEnd();++iter2){
			if(iter2->value.IsInt()){
				std::cout<<iter2->name.GetString()<<" [int]: "<<iter2->value.GetInt()<<std::endl;}
			else if(iter2->value.IsString()){
				std::cout<<iter2->name.GetString()<<" [string]: "<<iter2->value.GetString()<<std::endl;}
			else if(iter2->value.IsBool()){
				std::cout<<iter2->name.GetString()<<" [bool]: "<<iter2->value.GetBool()<<std::endl;}
			else if(iter2->value.IsArray()){
				std::cout<<iter2->name.GetString()<<" [array]: "<<std::endl;;}}
		std::cout<<std::endl;

		std::string type=Get_String_Value(obj,"type");
		if(type=="Grid"){
			std::string name=Get_String_Value(obj,"name");
			auto grid=Add_Object<OpenGLGrid>(name);
			Set_OpenGL_Object_Visibility(obj,grid,'g');}
		else if(type=="Particles"){
			std::string name=Get_String_Value(obj,"name");
			auto p=Add_Object<OpenGLParticles<Particles<3> > >(name);
			Set_OpenGL_Object_Visibility(obj,p,'c');
			Set_OpenGL_Object_Color(obj,p);
			Set_OpenGL_Object_Point_Size(obj,p);
			Set_OpenGL_Object_Particle_Mode(obj,p);
			Set_OpenGL_Object_Particle_Size(obj,p);
			Set_OpenGL_Object_Particle_Velocity_Visibility(obj,p);
			Set_OpenGL_Object_Particle_Force_Visibility(obj,p);}
		else if(type=="SegmentMesh"){
			std::string name=Get_String_Value(obj,"name");
			auto seg=Add_Object<OpenGLSegmentMesh>(name);
			Set_OpenGL_Object_Visibility(obj,seg,'s');
			Set_OpenGL_Object_Line_Width(obj,seg);
			Set_OpenGL_Object_Shading_Mode(obj,seg);
			Set_OpenGL_Object_Color(obj,seg);
			Set_OpenGL_Object_Polygon_Mode(obj,seg);}
		else if(type=="TriangleMesh"){
			std::string name=Get_String_Value(obj,"name");
			auto tri=Add_Object<OpenGLTriangleMesh>(name);
			Set_OpenGL_Object_Visibility(obj,tri,'t');
			Set_OpenGL_Object_Shading_Mode(obj,tri);
			Set_OpenGL_Object_Color(obj,tri);		
			Set_OpenGL_Object_Polygon_Mode(obj,tri);}
		else if(type=="TetrahedronMesh"){
			std::string name=Get_String_Value(obj,"name");
			auto tet=Add_Object<OpenGLTetrahedronMesh>(name);
			Set_OpenGL_Object_Visibility(obj,tet,'t');
			Set_OpenGL_Object_Shading_Mode(obj,tet);
			Set_OpenGL_Object_Color(obj,tet);		
			Set_OpenGL_Object_Polygon_Mode(obj,tet);}

		else if(type=="ScalarField"){
			std::string name=Get_String_Value(obj,"name");
			std::string grid_name=Get_String_Value(obj,"grid_name");
			int c_type=Get_Color_Type(obj);
			int s_type=Get_Storage_Type(obj);
			auto sca=Add_Grid_Scalar_Field(grid_name,{{name,"",(ColorType)c_type,(StoreType)s_type}});
			Set_OpenGL_Object_Visibility(obj,sca,'a');}
		else if(type=="VectorField"){
			std::string name=Get_String_Value(obj,"name");
			std::string grid_name=Get_String_Value(obj,"grid_name");
			int c_type=Get_Color_Type(obj);
			int s_type=Get_Storage_Type(obj);
			auto vec=Add_Grid_Vector_Field(grid_name,{{name,"",(ColorType)c_type,(StoreType)s_type}});
			Set_OpenGL_Object_Visibility(obj,vec,'v');}

	}
}

int OpenGLViewerConfiguable::Get_Color_Type(const JsonObj& obj)
{
	std::string color_type=Get_String_Value(obj,"color_type");
	int c_type=0;
	if(color_type=="Jet"||color_type=="jet")c_type=0;
	else if(color_type=="Hot"||color_type=="hot")c_type=1;
	else if(color_type=="Den"||color_type=="den")c_type=2;
	else if(color_type=="Mat"||color_type=="mat")c_type=3;	
	return c_type;
}

int OpenGLViewerConfiguable::Get_Storage_Type(const JsonObj& obj)
{
	std::string store_type=Get_String_Value(obj,"store_type");
	int s_type=0;
	if(store_type=="Cell"||store_type=="cell")s_type=0;
	else if(store_type=="Node"||store_type=="node")s_type=1;
	return s_type;
}

#endif
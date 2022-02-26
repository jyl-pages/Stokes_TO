//////////////////////////////////////////////////////////////////////////
// OpenGLViewerConfigable
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
////OpenGLViewerConfiguable enables setting up a configurable rendering driver with a Json file
//////////////////////////////////////////////////////////////////////////

#ifndef __OpenGLViewerConfigable_h__
#define __OpenGLViewerConfigable_h__
#include "OpenGLViewer.h"
#ifdef USE_RAPIDJSON
#include "RapidJsonFunc.h"
#endif

class OpenGLViewerConfiguable : public OpenGLViewer
{typedef OpenGLViewer Base;
#ifdef USE_RAPIDJSON
public:
	virtual void Initialize_Data();

protected:
	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Visibility(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj,char default_key)
	{
		if(Has_Bool_Value(obj,"vis")){
			bool vis=Get_Bool_Value(obj,"vis");
			char key=Has_String_Value(obj,"vis_key")?Get_String_Value(obj,"vis_key")[0]:default_key;
			Set_Visibility(opengl_obj,key,vis);}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Color(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{	
		if(Has_Vector4f_Value(obj,"color")){
			Vector4f rgba=Get_Vector4f_Value(obj,"color");
			Set_Color(opengl_obj,OpenGLColor(rgba[0],rgba[1],rgba[2],rgba[3]));}
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Line_Width(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_Float_Value(obj,"line_width")){
			float line_width=Get_Float_Value(obj,"line_width");
			Set_Line_Width(opengl_obj,line_width);}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Point_Size(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_Float_Value(obj,"point_size")){
			float line_width=Get_Float_Value(obj,"point_size");
			Set_Point_Size(opengl_obj,line_width);}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Particle_Mode(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_String_Value(obj,"particle_mode")){
			std::string particle_mode=Get_String_Value(obj,"particle_mode");
			if(particle_mode=="Dot"){Set_Particle_Mode(opengl_obj,ParticleMode::Dot);}
			else if(particle_mode=="Circle"){Set_Particle_Mode(opengl_obj,ParticleMode::Circle);}
			else if(particle_mode=="Sphere"){Set_Particle_Mode(opengl_obj,ParticleMode::Sphere);}}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Particle_Size(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_Float_Value(obj,"particle_size")){
			float particle_size=Get_Float_Value(obj,"particle_size");
			Set_Point_Size(opengl_obj,particle_size);}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Particle_Velocity_Visibility(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_String_Value(obj,"pv_vis_key")){
			bool vis=Has_Bool_Value(obj,"pv_vis")?Get_Bool_Value(obj,"pv_vis"):false;
			char key=Get_String_Value(obj,"pv_vis_key")[0];
			Set_Particle_Velocity_Visibility(opengl_obj,key,vis);}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Particle_Force_Visibility(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_String_Value(obj,"pf_vis_key")){
			bool vis=Has_Bool_Value(obj,"pf_vis")?Get_Bool_Value(obj,"pf_vis"):false;
			char key=Get_String_Value(obj,"pf_vis_key")[0];
			Set_Particle_Force_Visibility(opengl_obj,key,vis);}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Polygon_Mode(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_String_Value(obj,"polygon_mode")){
			std::string polygon_mode=Get_String_Value(obj,"polygon_mode");
			if(polygon_mode=="Fill"){Set_Polygon_Mode(opengl_obj,PolygonMode::Fill);}
			else if(polygon_mode=="Wireframe"){Set_Polygon_Mode(opengl_obj,PolygonMode::Wireframe);}
			else if(polygon_mode=="SurfOnly"){Set_Polygon_Mode(opengl_obj,PolygonMode::SurfOnly);}}		
	}

	template<class T_OPENGL_OBJ> void Set_OpenGL_Object_Shading_Mode(const JsonObj& obj,T_OPENGL_OBJ& opengl_obj)
	{
		if(Has_String_Value(obj,"shading_mode")){
			std::string shading_mode=Get_String_Value(obj,"shading_mode");
			if(shading_mode=="None"){Set_Shading_Mode(opengl_obj,ShadingMode::None);}
			else if(shading_mode=="Lighting"){Set_Shading_Mode(opengl_obj,ShadingMode::Lighting);}
			else if(shading_mode=="TexOnly"){Set_Shading_Mode(opengl_obj,ShadingMode::TexOnly);}
			else if(shading_mode=="TexLighting"){Set_Shading_Mode(opengl_obj,ShadingMode::TexLighting);}
			else if(shading_mode=="Shadow"){Set_Shading_Mode(opengl_obj,ShadingMode::Shadow);}
			else if(shading_mode=="Sprite"){Set_Shading_Mode(opengl_obj,ShadingMode::Sprite);}}		
	}

	int Get_Color_Type(const JsonObj& obj);
	int Get_Storage_Type(const JsonObj& obj);
#endif
};

#endif
//////////////////////////////////////////////////////////////////////////
// Opengl particles
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLTrackerPoints_h__
#define __OpenGLTrackerPoints_h__
#include "Common.h"
#include "Particles.h"
#include "File.h"
#include "ArrayIO.h"
#include "OpenGLObject.h"
#include "OpenGLShaderLibrary.h"

//////////////////////////////////////////////////////////////////////////
////This class is for the visualization of large-scale point clouds

class OpenGLTrackerPoints : public OpenGLObject
{typedef OpenGLObject Base;
public:
	GLfloat point_size=6.f;
	bool use_point_color=false;

	OpenGLTrackerPoints(){color=OpenGLColor::Red();name="tracker_points";}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_vcolor"));
	}

	virtual void Update_Data_To_Render()
	{
		if(!initialized)Initialize();
		if(!data_refreshed)return;

		if(use_point_color){
			GLuint stride_size=8;	
			vtx_size=(int)opengl_vertices.size();
			Set_OpenGL_Vertices();
			Set_OpenGL_Vertex_Attribute(0,4,stride_size,0);		////position
			Set_OpenGL_Vertex_Attribute(1,4,stride_size,4);		////color
		}
		else {
			GLuint stride_size = 4;
			vtx_size = (int)opengl_vertices.size();
			Set_OpenGL_Vertices();
			Set_OpenGL_Vertex_Attribute(0, 4, stride_size, 0);
		}////position
		
		Clear_OpenGL_Arrays();
		Set_Data_Refreshed(false);
	}

	virtual void Refresh(const int frame)
	{
		std::string file_name=output_dir+"/"+std::to_string(frame)+"/"+name;
		if (File::File_Exists(file_name)) {
			BinaryDataIO::Read_Scalar_Array(file_name, opengl_vertices);
			Set_Data_Refreshed();
		}
	}

	virtual void Display() const
    {
    	using namespace OpenGLUbos;using namespace OpenGLTextures;
		if(!visible)return;

		if(use_point_color){
			Update_Polygon_Mode();
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[1];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			glEnable(GL_PROGRAM_POINT_SIZE);
			shader->Set_Uniform("point_size",point_size);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/8);
			shader->End();		
		}else{
			Update_Polygon_Mode();
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			shader->Set_Uniform_Vec4f("color",color.rgba);
			glEnable(GL_PROGRAM_POINT_SIZE);
			shader->Set_Uniform("point_size",point_size);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/4);
			shader->End();}
    }

	virtual void Set_Color(const OpenGLColor& c)
	{color=c;}

	virtual void Set_Point_Size(const GLfloat _point_size)
	{point_size=_point_size;}

	void Use_Point_Color(const bool _use_point_color=true)
	{use_point_color=_use_point_color;}
};


#endif
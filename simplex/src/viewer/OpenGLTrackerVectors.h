//////////////////////////////////////////////////////////////////////////
// Opengl tracker vectors
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLTrackerVectors_h__
#define __OpenGLTrackerVectors_h__
#include "Common.h"
#include "Particles.h"
#include "File.h"
#include "OpenGLObject.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLUbos.h"
#include "OpenGLTexture.h"

//////////////////////////////////////////////////////////////////////////
////This class is for the visualization of large-scale vector clouds

class OpenGLTrackerVectors : public OpenGLObject
{typedef OpenGLObject Base;
public:
	OpenGLTrackerVectors(){color=OpenGLColor::Red();name="tracker_vectors";}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));
	}

	virtual void Update_Data_To_Render()
	{
		if(!initialized)Initialize();
		if(!data_refreshed)return;

		GLuint stride_size=4;	////currently does not support point color
		vtx_size=(int)opengl_vertices.size();
		Set_OpenGL_Vertices();
		int idx=0;{Set_OpenGL_Vertex_Attribute(0,4,stride_size,0);idx++;}	////position
		
		Set_Data_Refreshed(false);
	}

	virtual void Refresh(const int frame)
	{
		Clear_OpenGL_Arrays();
		std::string file_name=output_dir+"/"+std::to_string(frame)+"/"+name;
		if (File::File_Exists(file_name)) {
			std::ifstream input(file_name, std::ios::binary); if (!input)return;
			int n = 0; File::Read_Binary(input, n);
			opengl_vertices.resize(n);
			File::Read_Binary_Array(input, &opengl_vertices[0], n);
			if (scale != (real)1)Rescale();
			Set_Data_Refreshed();
		}
	}

	virtual void Display() const
    {
    	using namespace OpenGLUbos;
		if(!visible)return;
		if (shading_mode == ShadingMode::ColorEncoding) {
			Update_Polygon_Mode();
			{std::shared_ptr<OpenGLShaderProgram> shader = shader_programs[0];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader, "camera");
			int n = (int)opengl_vertices.size() / 8;
			for (int i = 0; i < n; i++) {
				OpenGLColor color_code;
				if (i == color_code_selection) color_code = OpenGLColor::White();
				else color_code = ID_To_Color(i + color_code_offset);
				shader->Set_Uniform_Vec4f("color", color_code.rgba);
				Base::Set_OpenGL_Vertices(opengl_vertices, 8, i * 8);
				
				glBindVertexArray(vao);
				glDrawArrays(GL_LINES, 0, (GLsizei)(2));
			}
			shader->End(); }
		}
		else {
			Update_Polygon_Mode();
			{std::shared_ptr<OpenGLShaderProgram> shader = shader_programs[0];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader, "camera");
			shader->Set_Uniform_Vec4f("color", color.rgba);
			glBindVertexArray(vao);
			glDrawArrays(GL_LINES, 0, (GLsizei)(vtx_size / 4));
			shader->End(); }
		}
    }

	virtual void Set_Color(const OpenGLColor& c)
	{
		color=c;
	}

	void Rescale()
	{
		int n=(int)opengl_vertices.size()/8;
		for(int i=0;i<n;i++){
			opengl_vertices[i*8+4]=opengl_vertices[i*8]+(opengl_vertices[i*8+4]-opengl_vertices[i*8])*(GLfloat)scale;
			opengl_vertices[i*8+5]=opengl_vertices[i*8+1]+(opengl_vertices[i*8+5]-opengl_vertices[i*8+1])*(GLfloat)scale;
			opengl_vertices[i*8+6]=opengl_vertices[i*8+2]+(opengl_vertices[i*8+6]-opengl_vertices[i*8+2])*(GLfloat)scale;}
	}

	virtual int Color_Code_Size(void) {
		return (int)opengl_vertices.size() / 8;
	}

	virtual void Output_Select_Info(int idx) {
		std::cout << name << " select element " << idx<<" with vector: ";
		static GLfloat vec[3];
		for (int axis = 0; axis < 3; axis++) {
			vec[axis] = opengl_vertices[idx * 8 + axis + 4] - opengl_vertices[idx * 8 + axis];
			std::cout << vec[axis] << " ";
		}
		real norm = sqrt(vec[0] * vec[0] + vec[1] * vec[1] + vec[2] * vec[2]);
		std::cout << "with norm = " << norm << std::endl;
	}

	virtual void Output_Debug_Info(std::ostream& out) {
		Base::Output_Debug_Info(out);
		out << ", size: " << Color_Code_Size();
	}
};
#endif
//////////////////////////////////////////////////////////////////////////
// Opengl points
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLPoints_h__
#define __OpenGLPoints_h__
#include "Common.h"
#include "Particles.h"
#include "File.h"
#include "OpenGLObject.h"
#include "OpenGLVectors.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLColorMapper.h"
#include "OpenGLTexture.h"

class OpenGLPoints : public OpenGLObject
{typedef OpenGLObject Base;
public:
    const Array<Vector3>* points=nullptr;
	const Array<real>* colors=nullptr;

	GLfloat point_size=6.f;
	bool use_varying_point_size=false;
	Array<GLfloat> varying_point_size;

	OpenGLPoints(){color=OpenGLColor::Red();name="points";}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));	
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("sprite"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vcolor_zsize_sprite"));	
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vcolor_psize_sprite"));
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vcolor_psize_ptex_sprite"));	
	}

	void Set_Data_Pointers(const Array<Vector3>* _points,const Array<real>* _colors=nullptr){points=_points;colors=_colors;}
	
	virtual void Update_Data_To_Render()
	{
		if(!initialized)Initialize();

		use_vtx_color=(colors!=nullptr&&shading_mode!=ShadingMode::None);
		GLuint stride_size=4+(use_vtx_color?4:0);
		Clear_OpenGL_Arrays();
		GLfloat placeholder=(GLfloat)0;
		for(auto i=0;i<(*points).size();i++){
			OpenGL_Vertex((*points)[i],opengl_vertices);		////position, 3 floats
			if(use_varying_point_size) OpenGL_Vertex(varying_point_size[i],opengl_vertices);
			else OpenGL_Vertex(placeholder,opengl_vertices);	////placeholder, 1 float
			if(use_vtx_color){
				OpenGLColor color=color_mapper->Color((*colors)[i]);
				OpenGL_Color4(color.rgba,opengl_vertices);}}	////color, 4 floats
		
		Set_OpenGL_Vertices();
		int idx=0;{Set_OpenGL_Vertex_Attribute(0,4,stride_size,0);idx++;}	////position
		if(use_vtx_color){Set_OpenGL_Vertex_Attribute(idx,4,stride_size,idx*4);idx++;}	////color
		Clear_OpenGL_Arrays();
	}

	virtual void Display() const
    {
    	using namespace OpenGLUbos;using namespace OpenGLTextures;
		if(!visible)return;
		Update_Polygon_Mode();

		switch(shading_mode){
		case ShadingMode::None:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			shader->Set_Uniform_Vec4f("color",color.rgba);
			glEnable(GL_PROGRAM_POINT_SIZE);
			shader->Set_Uniform("point_size",point_size);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/4);
			shader->End();
		}break;
		case ShadingMode::Sprite:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[1];
			shader->Begin();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SPRITE);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/4);
			shader->End();			
		}break;
		case ShadingMode::ColorSprite:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[2];
			shader->Begin();
			//Enable_Alpha_Blend();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SPRITE);
			shader->Set_Uniform("point_size",point_size);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/8);
			//Disable_Alpha_Blend();
			shader->End();			
		}break;
		case ShadingMode::SizeSprite:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[3];
			shader->Begin();
			Enable_Alpha_Blend();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SPRITE);
			shader->Set_Uniform("point_size",point_size);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/8);
			Disable_Alpha_Blend();
			shader->End();		
		}break;
		case ShadingMode::TexSprite:{
			std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[4];
			shader->Begin();
			Enable_Alpha_Blend();
			Bind_Uniform_Block_To_Ubo(shader,"camera");
			OpenGLTextures::Bind_Texture("sprite2.png",TextureType::Tx2d);
			shader->Set_Uniform("tex2d",0);
			glEnable(GL_VERTEX_PROGRAM_POINT_SIZE);
			glEnable(GL_POINT_SPRITE);
			shader->Set_Uniform("point_size",point_size);
			glBindVertexArray(vao);
			glDrawArrays(GL_POINTS,0,vtx_size/8);
			Disable_Alpha_Blend();
			shader->End();			
		}break;
		}
    }
};
#endif
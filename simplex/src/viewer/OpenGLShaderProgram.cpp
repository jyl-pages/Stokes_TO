//////////////////////////////////////////////////////////////////////////
// Opengl shader program
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "OpenGLShaderProgram.h"

void OpenGLShaderProgram::Initialize(const std::string& vtx_shader_input,const std::string& frg_shader_input)
{vtx_shader=vtx_shader_input;frg_shader=frg_shader_input;compiled=false;
vtx_id=0;frg_id=0;geo_id=0;prg_id=0;
use_geo=false;geo_input_type=GL_POINTS;geo_output_type=GL_TRIANGLE_STRIP;}

void OpenGLShaderProgram::Initialize(const std::string& vtx_shader_input,const std::string& frg_shader_input,
    const std::string& _geo_shader_input,GLenum _geo_input_type,GLenum _geo_output_type,int _max_geo_vtx_output)
{Initialize(vtx_shader_input,frg_shader_input);
use_geo=true;geo_shader=_geo_shader_input;geo_input_type=_geo_input_type;
geo_output_type=_geo_output_type;max_geo_vtx_output=_max_geo_vtx_output;}
	
void OpenGLShaderProgram::Set_Uniform(const std::string& name,GLint value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1i(location,value);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,GLfloat value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1f(location,value);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,Vector2f value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform2f(location,value[0],value[1]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,Vector3f value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform3f(location,value[0],value[1],value[2]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,Vector4f value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform4f(location,value[0],value[1],value[2],value[3]);}

void OpenGLShaderProgram::Set_Uniform(const std::string& name,glm::vec2 value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform2f(location,value[0],value[1]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,glm::vec3 value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform3f(location,value[0],value[1],value[2]);}
void OpenGLShaderProgram::Set_Uniform(const std::string& name,glm::vec4 value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform4f(location,value[0],value[1],value[2],value[3]);}

void OpenGLShaderProgram::Set_Uniform_Array(const std::string& name,GLsizei count,const GLint* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1iv(location,count,value);}
void OpenGLShaderProgram::Set_Uniform_Array(const std::string& name,GLsizei count,const GLfloat* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform1fv(location,count,value);}

void OpenGLShaderProgram::Set_Uniform_Matrix4f(const std::string& name,const GLfloat* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniformMatrix4fv(location,1,GL_FALSE,value);}
void OpenGLShaderProgram::Set_Uniform_Vec4f(const std::string& name,const GLfloat* value)
{GLint location=glGetUniformLocation(prg_id,name.c_str());glUniform4f(location,value[0],value[1],value[2],value[3]);}

void OpenGLShaderProgram::Set_Uniform_Mat(const Material* mat)
{Set_Uniform("mat_amb",mat->mat_amb);Set_Uniform("mat_dif",mat->mat_dif);Set_Uniform("mat_spec",mat->mat_spec);Set_Uniform("mat_shinness",mat->mat_shinness);}

void OpenGLShaderProgram::Bind_Uniform_Block(const std::string& name,const GLuint binding_point)
{GLuint location=glGetUniformBlockIndex(prg_id,name.c_str());glUniformBlockBinding(prg_id,location,binding_point);}

void OpenGLShaderProgram::Bind_Texture2D(const std::string& name,GLuint tex_id,GLint tex_unit)
{GLuint location=glGetUniformLocation(prg_id,name.c_str());glActiveTexture(GL_TEXTURE0+tex_unit);
glBindTexture(GL_TEXTURE_2D,tex_id);glUniform1i(location,tex_unit);}

void OpenGLShaderProgram::Begin(){if(!compiled)Compile();glUseProgram(prg_id);}
void OpenGLShaderProgram::End(){glUseProgram(0);}

bool OpenGLShaderProgram::Compile()
{
	if(compiled)return true;

	vtx_id=glCreateShader(GL_VERTEX_SHADER);
	const char* vtx_shader_string=vtx_shader.c_str();
	GLint vtx_string_length=(GLint)vtx_shader.length()+1;
	glShaderSource(vtx_id,1,&vtx_shader_string,&vtx_string_length);
	glCompileShader(vtx_id);
	GLint vtx_compile_status;
	glGetShaderiv(vtx_id,GL_COMPILE_STATUS,&vtx_compile_status);
	if(vtx_compile_status!=GL_TRUE){
		char log[2048];int log_length;
		glGetShaderInfoLog(vtx_id,2048,(GLsizei*)&log_length,log);
		std::cerr<<"Error: [OpenGLShaderProgram] vertex shader compile error: "<<log<<std::endl;
		glDeleteShader(vtx_id);
		return false;}

	frg_id=glCreateShader(GL_FRAGMENT_SHADER);
	const char* frg_shader_string=frg_shader.c_str();
	GLint frg_string_length=(GLint)frg_shader.length()+1;
	glShaderSource(frg_id,1,&frg_shader_string,&frg_string_length);
	glCompileShader(frg_id);
	GLint frg_compile_status;
	glGetShaderiv(frg_id,GL_COMPILE_STATUS,&frg_compile_status);
	if(frg_compile_status!=GL_TRUE){
		char log[2048];int log_length;
		glGetShaderInfoLog(frg_id,2048,(GLsizei*)&log_length,log);
		std::cerr<<"Error: [OpenGLShaderProgram] fragment shader compile error: "<<log<<std::endl;
		glDeleteShader(frg_id);
		return false;}

	prg_id=glCreateProgram();
	glAttachShader(prg_id,vtx_id);
	glAttachShader(prg_id,frg_id);

	if(use_geo){
		geo_id=glCreateShader(GL_GEOMETRY_SHADER_EXT);
		const char* geo_shader_string=geo_shader.c_str();
		GLint geo_string_length=(GLint)geo_shader.length()+1;
		glShaderSource(geo_id,1,&geo_shader_string,&geo_string_length);
		glCompileShader(geo_id);
		GLint geo_compile_status;
		glGetShaderiv(geo_id,GL_COMPILE_STATUS,&geo_compile_status);
		if(geo_compile_status!=GL_TRUE){
			char log[2048];int log_length;
			glGetShaderInfoLog(geo_id,2048,(GLsizei*)&log_length,log);
			std::cerr<<"Error: [OpenGLShaderProgram] geometry shader compile error: "<<log<<std::endl;
			glDeleteShader(geo_id);
			return false;}

		glAttachShader(prg_id,geo_id);
		glProgramParameteriEXT(prg_id,GL_GEOMETRY_INPUT_TYPE_EXT,geo_input_type);
		glProgramParameteriEXT(prg_id,GL_GEOMETRY_OUTPUT_TYPE_EXT,geo_output_type);
		glProgramParameteriEXT(prg_id,GL_GEOMETRY_VERTICES_OUT_EXT,max_geo_vtx_output);}

	glLinkProgram(prg_id);
	GLint prg_link_status;
	glGetProgramiv(prg_id,GL_LINK_STATUS,&prg_link_status);
	if(prg_link_status!=GL_TRUE){
		char log[2048];int log_length;
		glGetShaderInfoLog(prg_id,2048,(GLsizei*)&log_length,log);
		std::cerr<<"Error: [OpenGLShaderProgram] program link error: "<<log<<std::endl;
		glDeleteProgram(prg_id);
		return false;}

	glDeleteShader(vtx_id);
	glDeleteShader(frg_id);
	if(use_geo)glDeleteShader(geo_id);
	compiled=true;
	return true;
}
//////////////////////////////////////////////////////////////////////////
// Opengl shader library
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLShaderLibrary_h__
#define __OpenGLShaderLibrary_h__
#include <memory>
#include "Hashtable.h"

class OpenGLShaderProgram;

class OpenGLShaderLibrary
{public:
	static OpenGLShaderLibrary* Instance();
	static std::shared_ptr<OpenGLShaderProgram> Get_Shader(const std::string& name);
	static std::shared_ptr<OpenGLShaderProgram> Get_Shader_From_File(const std::string& name,const std::string& vtx_shader_file_name,const std::string& frg_shader_file_name);

	std::shared_ptr<OpenGLShaderProgram> Get(const std::string& name);
	std::shared_ptr<OpenGLShaderProgram> Get_From_File(const std::string& name,const std::string& vtx_shader_file_name,const std::string& frg_shader_file_name);

	void Add_Shader(const std::string& vtx_shader,const std::string& frg_shader,const std::string& name);
	void Add_Shader_From_File(const std::string& vtx_shader_file_name,const std::string& frg_shader_file_name,const std::string& name);

protected:
	Hashtable<std::string,std::shared_ptr<OpenGLShaderProgram> > shader_hashtable;
	Hashtable<std::string,std::string> shader_header_hashtable;

	OpenGLShaderLibrary();
	void Initialize_Shaders();
	void Initialize_Headers();
	std::string Parse(const std::string& shader) const;
};

#endif
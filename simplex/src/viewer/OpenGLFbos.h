//////////////////////////////////////////////////////////////////////////
// Opengl framebuffer object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLFbo_h__
#define __OpenGLFbo_h__
#include <string>
#include <GL/glew.h>
#include "glm.hpp"
#include "Hashtable.h"

namespace OpenGLFbos{
class OpenGLFbo
{public:
	GLuint buffer_index=0;
	std::string name="";
	GLsizei width=0;
	GLsizei height=0;
	GLuint tex_index=0;
	GLuint rbo_index=0;

	virtual void Initialize(const std::string& _name,GLuint _width=0,GLuint _height=0){}
	virtual void Resize(GLuint width,GLuint height){}
	virtual void Resize_To_Window(){}
	virtual void Write_To_File(const std::string& file_name){}
	virtual void Bind_As_Texture(GLuint idx=0){}
	virtual void Bind(){}
	virtual void Unbind(){}
	virtual void Clear(){}
	virtual void Set_Near_And_Far_Plane(float _near,float _far){}
};

class OpenGLFboInstance : public OpenGLFbo
{public:typedef OpenGLFbo Base;
	enum AttachmentType{Color=0,Depth,Position,Stencil} type;
	bool use_linearize_plane=false;
	float near_plane=.01f,far_plane=10.f;	////for depth linearization

	OpenGLFboInstance(AttachmentType _type=Color):Base(),type(_type){}
	virtual void Initialize(const std::string& _name,GLuint _width=0,GLuint _height=0);
	virtual void Resize(GLuint _width,GLuint _height);
	virtual void Clear();
	virtual void Resize_To_Window();
	virtual void Write_To_File(const std::string& file_name);
	virtual void Bind_As_Texture(GLuint idx=0);
	virtual void Set_Near_And_Far_Plane(float _near,float _far);
	virtual void Bind();
	virtual void Unbind();

protected:
	GLuint Generate_Attachment_Texture(const AttachmentType& att_type,GLuint width,GLuint height);
	void Bind_Texture_Color();
	void Bind_Texture_Depth();
	float Linearize_Depth(float depth,float near_plane,float far_plane);
};

class Fbo_Library
{public:
	static Fbo_Library* Instance();
	std::shared_ptr<OpenGLFbo> Get(const std::string& name,const int init_type=0);
protected:
	Hashtable<std::string,std::shared_ptr<OpenGLFbo> > fbo_hashtable;
	std::shared_ptr<OpenGLFbo> Lazy_Initialize_Fbo(const std::string& name,const int type);
};

std::shared_ptr<OpenGLFbo> Get_Fbo(const std::string& name,const int init_type=0);	////0-color,1-depth
std::shared_ptr<OpenGLFbo> Get_Depth_Fbo(const std::string& name);	////0-color,1-depth
OpenGLFboInstance* Get_Fbo_Instance(const std::string& name,const int init_type=0);
void Bind_Fbo(const std::string& name,const int init_type=0);
std::shared_ptr<OpenGLFbo> Get_And_Bind_Fbo(const std::string& name,const int init_type=0);
void Unbind_Fbo();
void Clear_Fbo(const std::string& name,const int init_type=0);
};
#endif
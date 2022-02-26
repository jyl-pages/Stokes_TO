//////////////////////////////////////////////////////////////////////////
// Opengl framebuffer object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "StbImage.h"
#include "OpenGLFbos.h"
#include "OpenGLWindow.h"

namespace OpenGLFbos{
void OpenGLFboInstance::Initialize(const std::string& _name,GLuint _width,GLuint _height)
{
	name=_name;
	glGenFramebuffers(1,&buffer_index);
	Resize(width,height);
}

void OpenGLFboInstance::Resize(GLuint _width,GLuint _height)
{
	if(_width==width&&_height==height||_width==0||_height==0)return;
	width=_width;height=_height;

	switch(type){
	case Color:case Position:Bind_Texture_Color();break;
	case Depth:Bind_Texture_Depth();break;}
}

void OpenGLFboInstance::Clear()
{
	Bind();
	glClearColor(0,0,0,0);
	glClearDepth(1.);
	glClear(GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);
	Unbind();
}

void OpenGLFboInstance::Resize_To_Window()
{
	if(width!=(GLuint)Win_Width()||height!=(GLuint)Win_Height())Resize((GLuint)Win_Width(),(GLuint)Win_Height());
}

void OpenGLFboInstance::Write_To_File(const std::string& file_name)
{
	glBindFramebuffer(GL_FRAMEBUFFER,buffer_index);
	switch(type){
	case Color:{
		int num_pixel=width*height;int num_comp=3;if(num_pixel==0)return;
		GLubyte* pixels=new GLubyte[num_comp*num_pixel];
		GLubyte* pixels_flipped_y=new GLubyte[num_comp*num_pixel];
		glReadPixels(0,0,width,height,GL_RGB,GL_UNSIGNED_BYTE,pixels);
		for(int i=0;i<height;i++){int offset=width*num_comp*(height-i-1);
			std::memcpy(pixels_flipped_y+offset,pixels+width*num_comp*i,width*num_comp);}
		std::stringstream ss;ss<<file_name<<".png";
		Stb::Write_Png(ss.str().c_str(),width,height,num_comp,pixels_flipped_y,0);
		delete pixels;delete pixels_flipped_y;		
	}break;
	case Depth:{
		int num_pixel=width*height;int num_comp=3;if(num_pixel==0)return;
		GLfloat* pixels=new GLfloat[num_pixel];
		GLubyte* pixels_flipped_y=new GLubyte[num_comp*num_pixel];
		glReadPixels(0,0,width,height,GL_DEPTH_COMPONENT,GL_FLOAT,pixels);
		for(int ii=0;ii<height;ii++){int i=height-ii-1;for(int j=0;j<width;j++){
			float depth=pixels[width*ii+j];
			if(use_linearize_plane)depth=Linearize_Depth(depth,near_plane,far_plane)/far_plane;
			int p=(int)(depth*255.);
			for(int k=0;k<num_comp;k++)pixels_flipped_y[width*num_comp*i+j*num_comp+k]=(GLubyte)p;}}
		std::stringstream ss;ss<<file_name<<".png";
		Stb::Write_Png(ss.str().c_str(),width,height,num_comp,pixels_flipped_y,0);
		delete pixels;delete pixels_flipped_y;			
	}break;}
	glBindFramebuffer(GL_FRAMEBUFFER,0);
}

void OpenGLFboInstance::Bind_As_Texture(GLuint idx){glActiveTexture(GL_TEXTURE0+idx);glBindTexture(GL_TEXTURE_2D,tex_index);}
void OpenGLFboInstance::Set_Near_And_Far_Plane(float _near,float _far){near_plane=_near;far_plane=_far;use_linearize_plane=true;}
void OpenGLFboInstance::Bind(){glBindFramebuffer(GL_FRAMEBUFFER,buffer_index);}
void OpenGLFboInstance::Unbind(){glBindFramebuffer(GL_FRAMEBUFFER,0);}

GLuint OpenGLFboInstance::Generate_Attachment_Texture(const AttachmentType& att_type,GLuint width,GLuint height)
{
	GLenum gl_att_type;GLenum gl_att_type_int;
	switch(att_type){
	case Color:{gl_att_type=gl_att_type_int=GL_RGB;}break;
	case Position:{gl_att_type=GL_RGB16F;gl_att_type_int=GL_RGB;}break;
	case Depth:{gl_att_type=gl_att_type_int=GL_DEPTH_COMPONENT;}break;
	case Stencil:{gl_att_type=GL_DEPTH24_STENCIL8;gl_att_type_int=GL_DEPTH_STENCIL;}break;}
	GLuint tex_id;glGenTextures(1,&tex_id);glBindTexture(GL_TEXTURE_2D,tex_id);
	switch(att_type){
	case Color:glTexImage2D(GL_TEXTURE_2D,0,gl_att_type,width,height,0,gl_att_type_int,GL_UNSIGNED_BYTE,0);break;
	case Position:glTexImage2D(GL_TEXTURE_2D,0,gl_att_type,width,height,0,gl_att_type_int,GL_FLOAT,0);break;
	case Depth:glTexImage2D(GL_TEXTURE_2D,0,gl_att_type,width,height,0,gl_att_type,GL_FLOAT,0);break;
	case Stencil:glTexImage2D(GL_TEXTURE_2D,0,gl_att_type,width,height,0,gl_att_type_int,GL_UNSIGNED_INT_24_8,0);break;}
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MIN_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_MAG_FILTER,GL_LINEAR);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_S,GL_CLAMP_TO_EDGE);
	glTexParameteri(GL_TEXTURE_2D,GL_TEXTURE_WRAP_T,GL_CLAMP_TO_EDGE);
	glBindTexture(GL_TEXTURE_2D,0);
	return tex_id;
}

void OpenGLFboInstance::Bind_Texture_Color()
{
	glBindFramebuffer(GL_FRAMEBUFFER,buffer_index);
	if(tex_index!=0)glDeleteTextures(1,&tex_index);
	tex_index=Generate_Attachment_Texture(type,width,height);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_COLOR_ATTACHMENT0,GL_TEXTURE_2D,tex_index,0);
	if(rbo_index!=0)glDeleteRenderbuffers(1,&rbo_index);
	glGenRenderbuffers(1,&rbo_index);glBindRenderbuffer(GL_RENDERBUFFER,rbo_index);
	glRenderbufferStorage(GL_RENDERBUFFER,GL_DEPTH24_STENCIL8,width,height);
	glBindRenderbuffer(GL_RENDERBUFFER,0);
	glFramebufferRenderbuffer(GL_FRAMEBUFFER,GL_DEPTH_STENCIL_ATTACHMENT,GL_RENDERBUFFER,rbo_index);
	if(glCheckFramebufferStatus(GL_FRAMEBUFFER)!=GL_FRAMEBUFFER_COMPLETE)std::cerr<<"Error: [OpenGLFboInstance] framebuffer not complete"<<std::endl;
	glBindFramebuffer(GL_FRAMEBUFFER,0);	
}

void OpenGLFboInstance::Bind_Texture_Depth()
{
	glBindFramebuffer(GL_FRAMEBUFFER,buffer_index);
	if(tex_index!=0)glDeleteTextures(1,&tex_index);
	tex_index=Generate_Attachment_Texture(Depth,width,height);
	glFramebufferTexture2D(GL_FRAMEBUFFER,GL_DEPTH_ATTACHMENT,GL_TEXTURE_2D,tex_index,0);
	glBindFramebuffer(GL_FRAMEBUFFER,0);		
}

float OpenGLFboInstance::Linearize_Depth(float depth,float near_plane,float far_plane)
{float z=depth*2.f-1.f;/*Back to NDC*/return (2.f*near_plane*far_plane)/(far_plane+near_plane-z*(far_plane-near_plane));}

Fbo_Library* Fbo_Library::Instance(){static Fbo_Library instance;return &instance;}

std::shared_ptr<OpenGLFbo> Fbo_Library::Get(const std::string& name,const int init_type)
{
	auto search=fbo_hashtable.find(name);
	if(search!=fbo_hashtable.end())return search->second;
	else return Lazy_Initialize_Fbo(name,init_type);
}

std::shared_ptr<OpenGLFbo> Fbo_Library::Lazy_Initialize_Fbo(const std::string& name,const int type)
{
	std::shared_ptr<OpenGLFbo> ptr=nullptr;if(name=="")return ptr;
	OpenGLFboInstance* fbo=new OpenGLFboInstance((OpenGLFboInstance::AttachmentType)type);
	fbo->Initialize(name);ptr.reset(fbo);
	fbo_hashtable.insert(std::make_pair(fbo->name,ptr));return ptr;
}

std::shared_ptr<OpenGLFbo> Get_Fbo(const std::string& name,const int init_type)	////0-color,1-depth
{return Fbo_Library::Instance()->Get(name,init_type);}

std::shared_ptr<OpenGLFbo> Get_Depth_Fbo(const std::string& name)	////0-color,1-depth
{return Fbo_Library::Instance()->Get(name,/*depth*/1);}

OpenGLFboInstance* Get_Fbo_Instance(const std::string& name,const int init_type)
{return dynamic_cast<OpenGLFboInstance*>(Get_Fbo(name,init_type).get());}

void Bind_Fbo(const std::string& name,const int init_type)
{auto fbo=Get_Fbo(name,init_type);if(fbo==nullptr)return;glBindFramebuffer(GL_FRAMEBUFFER,fbo->buffer_index);}

std::shared_ptr<OpenGLFbo> Get_And_Bind_Fbo(const std::string& name,const int init_type)
{auto fbo=Get_Fbo(name,init_type);if(fbo!=nullptr)glBindFramebuffer(GL_FRAMEBUFFER,fbo->buffer_index);return fbo;}

void Unbind_Fbo(){glBindFramebuffer(GL_FRAMEBUFFER,0);}

void Clear_Fbo(const std::string& name,const int init_type)
{
	auto fbo=Get_Fbo("depth",init_type);
	if(fbo==nullptr||fbo->width==0||fbo->height==0)return;
	fbo->Clear();
}
};
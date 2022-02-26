//////////////////////////////////////////////////////////////////////////
// Opengl uniform buffer object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "OpenGLUbos.h"

namespace OpenGLUbos{

//////////////////////////////////////////////////////////////////////////
////OpenGL UBO classes and library

void Bind_Shader_Ubo_Headers(Hashtable<std::string,std::string>& shader_header_hashtable)
{
	shader_header_hashtable.insert(std::make_pair("camera",camera));
	//shader_header_hashtable.insert(std::make_pair("lighting",lighting));
	shader_header_hashtable.insert(std::make_pair("lights",lights));
}

template<class T_UBO> void OpenGLUboInstance<T_UBO>::Initialize(const std::string& _uniform_block_name,GLuint _binding_point,GLuint _ubo,GLint _block_offset)
{
	name=_uniform_block_name;binding_point=_binding_point;
	if(_ubo!=0){buffer_index=_ubo;block_offset=_block_offset;}
	else{
		buffer_size=sizeof(T_UBO);
		glGenBuffers(1,&buffer_index);glBindBuffer(GL_UNIFORM_BUFFER,buffer_index);
		glBufferData(GL_UNIFORM_BUFFER,buffer_size,0,GL_STATIC_DRAW);
		glBindBuffer(GL_UNIFORM_BUFFER,0);	
		block_offset=0;block_size=buffer_size;}
	Bind_Block(binding_point,block_offset,block_size);
	std::cout<<"[OpenGLUboInstance] Initialized buffer_index: "<<buffer_index<<", name: "<<name<<", buffer_size: "<<buffer_size<<", binding_point: "<<binding_point<<
		", block_offset: "<<block_offset<<", block_size: "<<block_size<<std::endl;
}

template<class T_UBO> void OpenGLUboInstance<T_UBO>::Set_Block_Attributes(){Set_Block_Attributes_Helper(object,block_offset);}

template<class T_UBO> void OpenGLUboInstance<T_UBO>::Bind_Block(){Bind_Block(binding_point,block_offset,block_size);}

template<class T_UBO> void OpenGLUboInstance<T_UBO>::Bind_Block(GLuint binding_point,GLint block_offset,size_type block_size)
{glBindBufferRange(GL_UNIFORM_BUFFER,binding_point,buffer_index,block_offset,block_size);}

template<class T_UBO> void OpenGLUboInstance<T_UBO>::Set_Block_Attributes_Helper(const Camera& obj,GLint block_offset)
{
	GLint offset=block_offset;
	offset=Set_Block_Attribute(obj.projection,offset);
	offset=Set_Block_Attribute(obj.view,offset);
	offset=Set_Block_Attribute(obj.pvm,offset);
	offset=Set_Block_Attribute(obj.ortho,offset);
	offset=Set_Block_Attribute(obj.position,offset);
}

template<class T_UBO> void OpenGLUboInstance<T_UBO>::Set_Block_Attributes_Helper(const Lights& obj,GLint block_offset)
{
	GLint offset=block_offset;
	offset=Set_Block_Attribute(obj.amb,offset);
	offset=Set_Block_Attribute(obj.lt_att,offset);
	for(int i=0;i<2;i++){
		offset=Set_Block_Attribute(obj.lt[i].att,offset);
		offset=Set_Block_Attribute(obj.lt[i].pos,offset);
		offset=Set_Block_Attribute(obj.lt[i].dir,offset);
		offset=Set_Block_Attribute(obj.lt[i].amb,offset);
		offset=Set_Block_Attribute(obj.lt[i].dif,offset);
		offset=Set_Block_Attribute(obj.lt[i].spec,offset);
		offset=Set_Block_Attribute(obj.lt[i].atten,offset);
		offset=Set_Block_Attribute(obj.lt[i].r,offset);}
}

template class OpenGLUboInstance<Camera>;
template class OpenGLUboInstance<Lights>;

//////////////////////////////////////////////////////////////////////////
////OpenGL UBO classes and library

std::shared_ptr<OpenGLUbo> Ubo_Library::Get(const std::string& name)
{
	auto search=ubo_hashtable.find(name);
	if(search!=ubo_hashtable.end())return search->second;
	else return std::shared_ptr<OpenGLUbo>(nullptr);
}

GLuint Ubo_Library::Get_Binding_Point(const std::string& name)
{std::shared_ptr<OpenGLUbo> ptr=Get(name);if(ptr!=nullptr)return ptr->binding_point;else return GL_INVALID_INDEX;}

Ubo_Library::Ubo_Library(){Initialize_Ubos();}

void Ubo_Library::Initialize_Ubos()
{using namespace OpenGLUbos;
	std::cout<<"Initialize ubo library"<<std::endl;
	int binding_point=0;
	////camera
	{OpenGLUboInstance<Camera>* ubo=new OpenGLUboInstance<Camera>();
	ubo->Initialize("camera",binding_point++);;ubo->Set_Block_Attributes();
	ubo_hashtable.insert(std::make_pair(ubo->name,std::shared_ptr<OpenGLUbo>(ubo)));}
	////lights
	{OpenGLUboInstance<Lights>* ubo=new OpenGLUboInstance<Lights>();
	ubo->Initialize("lights",binding_point++);ubo->Set_Block_Attributes();
	ubo_hashtable.insert(std::make_pair(ubo->name,std::shared_ptr<OpenGLUbo>(ubo)));}
}

//////////////////////////////////////////////////////////////////////////
////UBO functions

void Initialize_Ubos()
{Ubo_Library::Instance()->Get("");}

std::shared_ptr<OpenGLUbo> Get_Ubo(const std::string& name)
{return Ubo_Library::Instance()->Get(name);}

GLuint Get_Ubo_Binding_Point(const std::string& name)
{return Ubo_Library::Instance()->Get_Binding_Point(name);}

bool Bind_Uniform_Block_To_Ubo(std::shared_ptr<OpenGLShaderProgram>& shader,const std::string& ubo_name)	////assuming uniform block name=ubo name
{GLuint binding_point=Get_Ubo_Binding_Point(ubo_name);if(binding_point==GL_INVALID_INDEX)return false;
shader->Bind_Uniform_Block(ubo_name,binding_point);return true;}

////Camera
OpenGLUboInstance<Camera>* Get_Camera_Ubo()
{
	std::shared_ptr<OpenGLUbo> ubo=Get_Ubo("camera");
	OpenGLUboInstance<Camera>* camera_ubo=dynamic_cast<OpenGLUboInstance<Camera>* >(ubo.get());return camera_ubo;	
}

Camera* Get_Camera()
{
	auto* camera_ubo=Get_Camera_Ubo();
	if(camera_ubo!=nullptr)return &camera_ubo->object;return nullptr;
}

Vector3 Get_Camera_Pos()
{
	auto camera=Get_Camera_Ubo();
	return Vector3(camera->object.position[0],camera->object.position[1],camera->object.position[2]);
}

////Lighting
OpenGLUboInstance<Lights>* Get_Lights_Ubo()
{
	std::shared_ptr<OpenGLUbo> ubo=Get_Ubo("lights");
	OpenGLUboInstance<Lights>* lights_ubo=dynamic_cast<OpenGLUboInstance<Lights>* >(ubo.get());return lights_ubo;
}

Lights* Get_Lights()
{
	auto* lights_ubo=Get_Lights_Ubo();
	if(lights_ubo!=nullptr)return &lights_ubo->object;return nullptr;
}

void Update_Lights_Ubo()
{auto* lights_ubo=Get_Lights_Ubo();lights_ubo->Set_Block_Attributes();}

Light* Get_Light(const int i)
{
	Lights* lights=Get_Lights();if(lights==nullptr)return nullptr;
	Light* light=lights->Get(i);return light;
}

void Clear_Lights()
{
	Lights* lights=Get_Lights();if(lights==nullptr)return;
	int& i=lights->Light_Num();i=0;
}

Light* Add_Directional_Light(const glm::vec3& dir)
{
	Lights* lights=Get_Lights();if(lights==nullptr)return nullptr;

	int& i=lights->Light_Num();Light* lt=&lights->lt[i];i++;
	lt->Initialize();
	lt->Set_Directional();
	lt->dir=glm::vec4(dir,1.f);
	lt->pos=glm::vec4(-dir*2.f,1.f);

	Update_Lights_Ubo();return lt;
}

Light* Add_Point_Light(const glm::vec3& pos)
{
	Lights* lights=Get_Lights();if(lights==nullptr)return nullptr;

	int& i=lights->Light_Num();Light* lt=&lights->lt[i];i++;
	lt->Initialize();
	lt->Set_Point();
	lt->pos=glm::vec4(pos,1.f);
	
	Update_Lights_Ubo();return lt;
}

Light* Add_Spot_Light(const glm::vec3& pos,glm::vec3& dir)
{
	Lights* lights=Get_Lights();if(lights==nullptr)return nullptr;

	int& i=lights->Light_Num();Light* lt=&lights->lt[i];i++;
	lt->Initialize();
	lt->Set_Spot();
	lt->pos=glm::vec4(pos,1.f);
	lt->dir=glm::vec4(dir,1.f);

	Update_Lights_Ubo();return lt;
}

void Set_Ambient(const glm::vec4& amb)
{
	Lights* lights=Get_Lights();if(lights==nullptr)return;
	lights->amb=amb;
}
};

//////////////////////////////////////////////////////////////////////////
// Opengl shader library
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "File.h"
#include "OpenGLShaders.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLShaderProgram.h"

using namespace OpenGLShaders;

OpenGLShaderLibrary* OpenGLShaderLibrary::Instance(){static OpenGLShaderLibrary instance;return &instance;}

std::shared_ptr<OpenGLShaderProgram> OpenGLShaderLibrary::Get(const std::string& name)
{
	auto search=shader_hashtable.find(name);
	if(search!=shader_hashtable.end())return search->second;
	else return std::shared_ptr<OpenGLShaderProgram>(nullptr);
}

std::shared_ptr<OpenGLShaderProgram> OpenGLShaderLibrary::Get_From_File(const std::string& name, 
	const std::string& vtx_shader_file_name,const std::string& frg_shader_file_name)
{
	if(shader_hashtable.find(name)==shader_hashtable.end())
		Add_Shader_From_File(vtx_shader_file_name,frg_shader_file_name,name);
	return Get(name);
}

OpenGLShaderLibrary::OpenGLShaderLibrary()
{Initialize_Shaders();}

void OpenGLShaderLibrary::Initialize_Shaders()
{
	Initialize_Headers();
	Add_Shader(vpos_vtx_shader,ucolor_frg_shader,"vpos");
	Add_Shader(vpos_model_vtx_shader,ucolor_frg_shader,"vpos_model");
	Add_Shader(vpos_model_vnormal_vfpos_vtx_shader,vnormal_vfpos_lt_frg_shader,"vpos_model_vnormal_lt");
	Add_Shader(vpos_model_vnormal_vfpos_vtx_shader,vnormal_vfpos_dl_fast_frg_shader,"vpos_model_vnormal_dl_fast");
	Add_Shader(vpos_model_vnormal_model_vfpos_vtx_shader,vnormal_vfpos_dl_fast_frg_shader,"vpos_model_vnormal_dl_fast2");

	Add_Shader(vcolor_vtx_shader,vcolor_frg_shader,"vcolor");
	Add_Shader(psize_vtx_shader,ucolor_frg_shader,"psize_ucolor");
	Add_Shader(psize_vcolor_vtx_shader,vcolor_frg_shader,"psize_vcolor");
	Add_Shader(vcolor_ortho_vtx_shader,vcolor_frg_shader,"ortho_vcolor");
	Add_Shader(vtex_vtx_shader,vtex_frg_shader,"vtex");
	Add_Shader(vnormal_vfpos_vtx_shader,vnormal_vfpos_dl_frg_shader,"vnormal_dl");
	Add_Shader(vnormal_vfpos_vtx_shader,vnormal_vfpos_pl_frg_shader,"vnormal_pl");
	Add_Shader(vnormal_vfpos_vtx_shader,vnormal_vfpos_lt_frg_shader,"vnormal_lt");
	Add_Shader(vnormal_vfpos_vtex_vtx_shader,vnormal_vfpos_vtex_lt_frg_shader,"vnormal_vtex_lt");
	Add_Shader(vnormal_vfpos_vtx_shader,vnormal_vfpos_lt_env_frg_shader,"vnormal_vtex_lt_env");
	Add_Shader(skybox_vtx_shader,skybox_frg_shader,"skybox");
	Add_Shader(sprite_vtx_shader,sprite_frg_shader,"sprite");
	Add_Shader(sprite_vcolor_vtx_shader,sprite_vcolor_frg_shader,"vcolor_zsize_sprite");
	Add_Shader(sprite_vcolor_vsize_vtx_shader,sprite_vcolor_frg_shader,"vcolor_psize_sprite");
	Add_Shader(sprite_vcolor_vsize_vtx_shader,sprite_tex_vcolor_frg_shader,"vcolor_psize_ptex_sprite");
	Add_Shader(vfpos_vtx_shader,vfpos_frg_shader,"aabb_ray_depth");
	Add_Shader(vfpos_vtx_shader,vol_frg_shader,"vol_ray_casting");
	Add_Shader(vfpos_vtx_shader,vec_vol_frg_shader,"vec_vol_ray_casting");
	Add_Shader(shadow_vtx_shader,none_frg_shader,"sd_depth");
	Add_Shader(vnormal_vfpos_vsdpos_vtx_shader,vnormal_vfpos_lt_sd_frg_shader,"sd_lt");
	Add_Shader(vclip_vtx_shader,ucolor_frg_shader,"ucolor_bk");
	Add_Shader(vclip_vfpos_vtx_shader,gcolor_frg_shader,"gcolor_bk");
	Add_Shader(vclip_vtex_vtx_shader,vtex_frg_shader,"vtex_bk");
}

void OpenGLShaderLibrary::Initialize_Headers()
{
	shader_header_hashtable.insert(std::make_pair("version",version));
	shader_header_hashtable.insert(std::make_pair("material",material));
	shader_header_hashtable.insert(std::make_pair("phong_dl_func",phong_dl_func));
	shader_header_hashtable.insert(std::make_pair("phong_pl_func",phong_pl_func));
	shader_header_hashtable.insert(std::make_pair("phong_sl_func",phong_sl_func));
	shader_header_hashtable.insert(std::make_pair("shadow_func",shadow_func));
	shader_header_hashtable.insert(std::make_pair("phong_dl_fast_func",phong_dl_fast_func));

	OpenGLUbos::Bind_Shader_Ubo_Headers(shader_header_hashtable);
}

std::string OpenGLShaderLibrary::Parse(const std::string& shader) const
{
	std::string s=shader;
	std::replace(s.begin(),s.end(),'~','#');	////replace ~ with #, fix for linux compile
	std::string name="#include";
	size_type p1=s.find(name);
	while(p1!=std::string::npos){
		size_type p2=s.find(' ',p1);
		size_type p3=s.find(';',p1);
		if(p2==std::string::npos||p3==std::string::npos)break;
		size_type n_var=p3-p2-1;
		std::string var=s.substr(p2+1,n_var);
		auto hash_pair=shader_header_hashtable.find(var);
		if(hash_pair==shader_header_hashtable.end())break;
		const std::string& replace=hash_pair->second;
		size_type n_replace=p3-p1+1;
		s.replace(p1,n_replace,replace);
		p1=s.find(name);
	}
	std::replace(s.begin(),s.end(),'~','#');	////replace ~ with #, fix for linux compile
	return s;
}

void OpenGLShaderLibrary::Add_Shader(const std::string& vtx_shader,const std::string& frg_shader,const std::string& name)
{
	std::shared_ptr<OpenGLShaderProgram> shader=std::make_shared<OpenGLShaderProgram>();
	shader->Initialize(Parse(vtx_shader),Parse(frg_shader));shader->name=name;
	shader_hashtable.insert(std::make_pair(shader->name,shader));	
}

void OpenGLShaderLibrary::Add_Shader_From_File(const std::string& vtx_shader_file_name,const std::string& frg_shader_file_name,const std::string& name)
{
	std::string vtx_shader;File::Read_Text_To_String(vtx_shader_file_name,vtx_shader);
	std::string frg_shader;File::Read_Text_To_String(frg_shader_file_name,frg_shader);
	std::cout<<"vtx_file: "<<vtx_shader_file_name<<std::endl;
	std::cout<<"frg_file: "<<frg_shader_file_name<<std::endl;

	std::cout<<"vtx:\n"<<vtx_shader<<std::endl;
	std::cout<<"frag:\n"<<frg_shader<<std::endl;
	Add_Shader(vtx_shader,frg_shader,name);
}

std::shared_ptr<OpenGLShaderProgram> OpenGLShaderLibrary::Get_Shader(const std::string& name)
{return OpenGLShaderLibrary::Instance()->Get(name);}

std::shared_ptr<OpenGLShaderProgram> OpenGLShaderLibrary::Get_Shader_From_File(const std::string& name, 
	const std::string& vtx_shader_file_name,const std::string& frg_shader_file_name)
{return OpenGLShaderLibrary::Instance()->Get_From_File(name,vtx_shader_file_name,frg_shader_file_name);}


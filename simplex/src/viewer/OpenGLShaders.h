//////////////////////////////////////////////////////////////////////////
// Opengl shaders
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLShaders_h__
#define __OpenGLShaders_h__
#include <string>
#include "Common.h"
#include "Hashtable.h"
#include "OpenGLUbos.h"

namespace OpenGLShaders{
#define To_String(S) #S

//////////////////////////////////////////////////////////////////////////
////predefined uniform blocks and functions
const std::string version=To_String(
~version 430 core\n
);

const std::string material=To_String(
uniform vec4 mat_amb=vec4(1.f);
uniform vec4 mat_dif=vec4(1.f,1.f,1.f,1.f);
uniform vec4 mat_spec=vec4(1.f);
uniform vec4 mat_shinness=vec4(32.f,0.f,0.f,0.f);
);

const std::string phong_dl_func=To_String(
vec3 phong_dl(int i,vec3 norm)
{
	vec4 amb=lt[i].amb*mat_amb;
	vec3 light_dir=lt[i].dir.xyz;
    float dif_coef=max(dot(norm,-light_dir),0.);
    vec4 dif=dif_coef*lt[i].dif*mat_dif;
	vec4 color=amb+dif;return color.rgb;
}
);

const std::string phong_pl_func=To_String(
vec3 phong_pl(int i,vec3 pos,vec3 norm)
{
	vec4 amb=lt[i].amb*mat_amb;
	vec3 light_dir=lt[i].pos.xyz-pos;float dis=length(light_dir);light_dir=light_dir/dis;
    float dif_coef=max(dot(norm,light_dir),0.f);
    vec4 dif=dif_coef*mat_dif;
	vec3 view_dir=normalize(position.xyz-pos);
	vec3 half_dir=normalize(light_dir+view_dir);
	float spec_coef=pow(max(dot(norm,half_dir),0.f),mat_shinness[0]);
	vec4 spec=spec_coef*lt[i].spec*mat_spec;
	
	vec4 color=amb+dif+spec;
	float atten_coef=1.f/(lt[i].atten[0]+lt[i].atten[1]*dis+lt[i].atten[2]*dis*dis);
	color*=atten_coef;
	
	return color.rgb;
}
);

const std::string phong_sl_func=To_String(
vec3 phong_sl(int i,vec3 pos,vec3 norm)
{
	vec4 amb=lt[i].amb*mat_amb;
	vec3 light_dir=lt[i].pos.xyz-pos;float dis=length(light_dir);light_dir=light_dir/dis;
	float theta=dot(light_dir,-lt[i].dir.xyz);
	float spot_coef=clamp((theta-lt[i].r[1])/lt[i].r[2],0.f,1.f);

    float dif_coef=max(dot(norm,light_dir),0.f);
    vec4 dif=dif_coef*mat_dif;
	vec3 view_dir=normalize(position.xyz-pos);
	vec3 half_dir=normalize(light_dir+view_dir);
	float spec_coef=pow(max(dot(norm,half_dir),0.f),mat_shinness[0]);
	vec4 spec=spec_coef*lt[i].spec*mat_spec;

	vec4 color=amb+(dif+spec)*spot_coef;
	float atten_coef=1.f/(lt[i].atten[0]+lt[i].atten[1]*dis+lt[i].atten[2]*dis*dis);
	color*=atten_coef;

	return color.rgb;
}
);

const std::string shadow_func=To_String(
float shadow(vec4 shadow_pos,vec3 normal,vec3 light_dir)
{
	vec3 proj_coord=shadow_pos.xyz/shadow_pos.w;
	proj_coord=proj_coord*.5f+.5f;
	
	float shadow=0.f;float dp=proj_coord.z;float step=1.f/512.f;
	float bias=max(.05f*(1.f-dot(normal,light_dir)),.005f);
	for(int i=-1;i<=1;i++)for(int j=-1;j<=1;j++){
		vec2 coord=proj_coord.xy+vec2(i,j)*step;
		float dp0=texture(shadow_map,coord).r;
		shadow+=dp>dp0+bias?0.2f:1.f;}shadow/=9.f;
	return shadow;
	
	/*float dp0=texture(shadow_map,proj_coord.xy).r;
	float dp=proj_coord.z;
	float bias=max(.0001f*(1.f-dot(normal,light_dir)),.00001f);
	float coef=dp>dp0+bias?0.f:1.f;return coef;*/
}
);

const std::string phong_dl_fast_func=To_String(
vec3 phong_dl_fast(vec3 norm)
{
    float dif_coef=abs(dot(norm,vec3(1.f,1.f,1.f)));
    vec4 dif=dif_coef*vec4(.5f)*mat_dif+vec4(.1f);
	return dif.rgb;
}
);

//////////////////////////////////////////////////////////////////////////
////vtx shader
const std::string vpos_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
}														
);

const std::string vpos_model_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
uniform mat4 model=mat4(1.0f);
void main()												
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
}														
);

const std::string vpos_model_vnormal_vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
	vtx_normal=vec3(normal);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}														
);

const std::string vpos_model_vnormal_model_vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
	vtx_normal=vec3(pvm*model*vec4(normal.xyz,1.f));
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}														
);

const std::string vcolor_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec4 vtx_color;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_color=v_color;
}														
);

const std::string vclip_vtx_shader=To_String(
~include version;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
void main()
{
	gl_Position=vec4(pos.xyz,1.f);
}
);

const std::string vclip_vfpos_vtx_shader=To_String(
~include version;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec3 vtx_frg_pos;
void main()
{
	gl_Position=vec4(pos.xyz,1.f);
	vtx_frg_pos=pos.xyz;
}
);

const std::string vclip_vtex_vtx_shader=To_String(
~include version;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 tex;
out vec2 vtx_frg_tex;
void main()
{
	gl_Position=vec4(pos.xyz,1.f);
	vtx_frg_tex=tex.xy;
}
);

const std::string vnormal_vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_normal=vec3(normal);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}
);

const std::string vtex_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 tex;
out vec2 vtx_frg_tex;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_frg_tex=tex.xy;
}
);

const std::string vnormal_vfpos_vtex_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
layout (location=2) in vec4 tex;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
out vec2 vtx_frg_tex;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_normal=vec3(normal);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
	vtx_frg_tex=tex.xy;
}
);

const std::string psize_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
uniform float point_size=1.f;
void main()												
{
	gl_PointSize=point_size;
	gl_Position=pvm*vec4(pos.xyz,1.f);
}														
);

const std::string psize_vcolor_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
uniform float point_size=1.f;
out vec4 vtx_color;
void main()												
{
	gl_PointSize=point_size;
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_color=v_color;
}														
);

const std::string vcolor_ortho_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec4 vtx_color;
void main()												
{
	gl_Position=ortho*model*vec4(pos.xy,1.f,1.f);
	vtx_color=v_color;
}														
);

const std::string skybox_vtx_shader=To_String(
~include version;
~include camera;
layout (location=0) in vec4 pos;
out vec3 vtx_frg_tex;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_frg_tex=pos.xyz;
}														
);

const std::string sprite_vtx_shader=To_String(
~include version;
~include camera;
uniform float point_size=1.f;
layout (location=0) in vec4 pos;

void main()												
{
	vec4 p=pvm*vec4(pos.xyz,1.f);
	gl_PointSize=(1.f-p.z/p.w)*point_size;
	gl_Position=p;
}														
);

const std::string sprite_vcolor_vtx_shader=To_String(
~include version;
~include camera;
uniform float point_size=1.f;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec4 vtx_color;
void main()												
{
	vec4 p=pvm*vec4(pos.xyz,1.f);
	gl_PointSize=(1.f-p.z/p.w)*point_size;
	gl_Position=p;
	vtx_color=v_color;
}														
);

const std::string sprite_vcolor_vsize_vtx_shader=To_String(
~include version;
~include camera;
uniform float point_size=1.f;
layout (location=0) in vec4 pos;
layout (location=1) in vec4 v_color;
out vec4 vtx_color;
void main()												
{
	gl_PointSize=pos.w*point_size;
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_color=v_color;
}														
);

const std::string vfpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 model=mat4(1.0f);
layout (location=0) in vec4 pos;
out vec3 vtx_frg_pos;
void main()												
{
	gl_Position=pvm*vec4(pos.xyz,1.f);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
}
);

const std::string shadow_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 shadow_pv;
uniform mat4 model=mat4(1.f);
layout (location=0) in vec4 pos;
void main()
{
    gl_Position=shadow_pv*model*vec4(pos.xyz,1.0);
}
);

const std::string vnormal_vfpos_vsdpos_vtx_shader=To_String(
~include version;
~include camera;
uniform mat4 shadow_pv;
uniform mat4 model=mat4(1.f);
layout (location=0) in vec4 pos;
layout (location=1) in vec4 normal;
out vec3 vtx_normal;
out vec3 vtx_frg_pos;
out vec4 vtx_shadow_pos;
void main()
{
	gl_Position=pvm*model*vec4(pos.xyz,1.f);
	vtx_normal=vec3(normal);
	vtx_frg_pos=vec3(model*vec4(pos.xyz,1.f));
    vtx_shadow_pos=shadow_pv*model*vec4(pos.xyz,1.f);
}
);

//////////////////////////////////////////////////////////////////////////
////frg shader
const std::string none_frg_shader=To_String(
~include version;
void main(){}
);

const std::string ucolor_frg_shader=To_String(
~include version;
uniform vec4 color=vec4(1.f,1.f,0.f,1.f);
out vec4 frag_color;
void main()								
{										
    frag_color=color;	  
}										
);

const std::string gcolor_frg_shader=To_String(
~include version;
uniform vec4 color_0=vec4(1.f,1.f,1.f,1.f);
uniform vec4 color_1=vec4(.8f,.9f,.8f,1.f);
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()								
{ 
	float m=abs(vtx_frg_pos.x);
	vec3 c=mix(color_0.rgb,color_1.rgb,m*m);
	frag_color=vec4(c,1.f);
}										
);

const std::string vcolor_frg_shader=To_String(
~include version;
in vec4 vtx_color;
out vec4 frag_color;
void main()								
{										
    frag_color=vtx_color;
}										
);

const std::string vtex_frg_shader=To_String(
~include version;
uniform sampler2D tex2d;
in vec2 vtx_frg_tex;
out vec4 frag_color;
void main()								
{										
    frag_color=texture(tex2d,vtx_frg_tex);	  
	//frag_color=vec4(vtx_frg_tex.x,vtx_frg_tex.y,0.f,1.f);
}										
);

const std::string vnormal_vfpos_dl_frg_shader=To_String(
~include version;
~include material;
~include lights;
~include phong_dl_func;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
    vec3 norm=normalize(vtx_normal);
	vec3 color=phong_dl(0,norm);
	frag_color=vec4(color,1.f);
}
);

const std::string vnormal_vfpos_dl_fast_frg_shader=To_String(
~include version;
~include material;
~include lights;
~include phong_dl_fast_func;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
    vec3 norm=normalize(vtx_normal);
	vec3 color=phong_dl_fast(norm);
	if(gl_FrontFacing) frag_color=vec4(color,1.f);
	else frag_color=vec4(1.f-color.r,1.f-color.g,1.f-color.b,1.f);
}
);

const std::string vnormal_vfpos_pl_frg_shader=To_String(
~include version;
~include material;
~include lights;
~include phong_dl_func;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
    vec3 norm=normalize(vtx_normal);
	vec3 color=phong_pl(0,norm);
	frag_color=vec4(color,1.f);
}
);

const std::string vnormal_vfpos_lt_frg_shader=To_String(
~include version;
~include material;
~include camera;
~include lights;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
~include phong_dl_func;
~include phong_pl_func;
~include phong_sl_func;
void main()
{
    vec3 normal=normalize(vtx_normal);
	vec3 color=mat_amb.rgb*amb.rgb;
	for(int i=0;i<lt_att[0];i++){
		vec3 c0=vec3(0.f);
		switch(lt[i].att[0]){
		case 0:{c0=phong_dl(i,normal);}break;
		case 1:{c0=phong_pl(i,vtx_frg_pos,normal);}break;
		case 2:{c0=phong_sl(i,vtx_frg_pos,normal);}break;}
		color+=c0;}
	if(gl_FrontFacing) frag_color=vec4(color,1.f);
	else frag_color=vec4(1.f-color.r,1.f-color.g,1.f-color.b,1.f);
}
);

const std::string vnormal_vfpos_lt_env_frg_shader=To_String(
~include version;
~include material;
~include camera;
~include lights;
uniform samplerCube tex3d;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
out vec4 frag_color;
~include phong_dl_func;
~include phong_pl_func;
~include phong_sl_func;
void main()
{
    vec3 normal=normalize(vtx_normal);
	vec3 color=mat_amb.rgb*amb.rgb;
	for(int i=0;i<lt_att[0];i++){
		vec3 c0=vec3(0.f);
		switch(lt[i].att[0]){
		case 0:{c0=phong_dl(i,normal);}break;
		case 1:{c0=phong_pl(i,vtx_frg_pos,normal);}break;
		case 2:{c0=phong_sl(i,vtx_frg_pos,normal);}break;}
		color+=c0;}
	
	vec3 view_dir=normalize(vtx_frg_pos-position.xyz);
	vec3 refl_dir=reflect(view_dir,normal);
	vec3 env_color=texture(tex3d,refl_dir).rgb;
	vec3 color2=min(mat_amb.rgb*env_color+color,1.f);
	
	frag_color=vec4(color2,1.f);
}
);

const std::string vnormal_vfpos_vtex_lt_frg_shader=To_String(
~include version;
~include material;
~include camera;
~include lights;
uniform sampler2D tex2d;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
in vec2 vtx_frg_tex;
out vec4 frag_color;
~include phong_dl_func;
~include phong_pl_func;
~include phong_sl_func;
void main()
{
    vec3 normal=normalize(vtx_normal);
	vec3 color=mat_amb.rgb*amb.rgb;
	for(int i=0;i<lt_att[0];i++){
		vec3 c0=vec3(0.f);
		switch(lt[i].att[0]){
		case 0:{c0=phong_dl(i,normal);}break;
		case 1:{c0=phong_pl(i,vtx_frg_pos,normal);}break;
		case 2:{c0=phong_sl(i,vtx_frg_pos,normal);}break;}
		color+=c0;}
	vec4 tex_color=texture(tex2d,vtx_frg_tex);
	frag_color=vec4(color,1.f)*tex_color;
}
);

const std::string vnormal_vfpos_lt_sd_frg_shader=To_String(
~include version;
~include material;
~include camera;
~include lights;
uniform sampler2D shadow_map;
~include shadow_func;
in vec3 vtx_normal;
in vec3 vtx_frg_pos;
in vec4 vtx_shadow_pos;
out vec4 frag_color;
~include phong_dl_func;
~include phong_pl_func;
~include phong_sl_func;
void main()
{
    vec3 normal=normalize(vtx_normal);
	vec3 color=mat_amb.rgb*amb.rgb;
	for(int i=0;i<lt_att[0];i++){
		vec3 c0=vec3(0.f);
		switch(lt[i].att[0]){
		case 0:{c0=phong_dl(i,normal);}break;
		case 1:{c0=phong_pl(i,vtx_frg_pos,normal);}break;
		case 2:{c0=phong_sl(i,vtx_frg_pos,normal);}break;}
		float s=1.f;
		if(lt[i].att[1]!=0){
			vec3 light_dir=lt[i].att[0]==0?-lt[i].dir.xyz:normalize(lt[i].pos.xyz-vtx_frg_pos);
			s=shadow(vtx_shadow_pos,normal,light_dir);}
		color+=c0*s;}
	frag_color=vec4(color,1.f);
}
);

const std::string skybox_frg_shader=To_String(
~include version;
uniform samplerCube tex3d;
in vec3 vtx_frg_tex;
out vec4 frag_color;
void main()								
{										
    frag_color=texture(tex3d,vtx_frg_tex);	  
}										
);

const std::string sprite_frg_shader=To_String(
~include version;
out vec4 frag_color;
void main()								
{									
	vec2 r=gl_PointCoord-vec2(.5f);
	float r2=dot(r,r);
	if(r2>.25f)discard;
    frag_color=vec4(1.f-4.f*r2,1.f-4.f*r2,0.f,1.f);	  
	//frag_color=vec4(0.f,1.f,0.f,1.f);
}	
);

const std::string sprite_vcolor_frg_shader=To_String(
~include version;
in vec4 vtx_color;
out vec4 frag_color;
void main()								
{									
	vec2 r=gl_PointCoord-vec2(.5f);
	float z=sqrt(.25f-dot(r,r));
	vec3 n=vec3(r.x,r.y,z);
	vec3 light_dir=vec3(1.f,1.f,1.f);
	float dif=max(0.f,dot(light_dir,n));

	vec3 light_dir2=vec3(-1.f,-2.f,0.f);
	float dif2=2.f*pow(max(0.f,dot(light_dir2,n)),1);

	vec3 light_dir3=vec3(0.f,-1.f,1.f);
	float dif3=1.f*pow(max(0.f,dot(light_dir3,n)),2);

	float r2=dot(r,r);
	if(r2>.25f)discard;
	/*float c=1.25f-4.f*r2;
	if(c<.2f)c=.2f;*/
    //frag_color=vec4(/*c**/vtx_color.xyz,.4f)*(dif+dif2+dif3);
	frag_color=vec4(1.f,1.f,1.f,1.f);
}
);

const std::string sprite_tex_vcolor_frg_shader=To_String(
~include version;
uniform sampler2D tex2d;
in vec4 vtx_color;
out vec4 frag_color;
void main()								
{									
	vec2 r=gl_PointCoord-vec2(.5f);
	float r2=dot(r,r);
	if(r2>.25f)discard;
	float alpha=1.f-r2*4.f;
	vec4 tex_color=texture(tex2d,gl_PointCoord);
	frag_color=vec4(tex_color.rgb*vtx_color.rgb,alpha);
}
);

const std::string vfpos_frg_shader=To_String(
~include version;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
	frag_color=vec4(vtx_frg_pos,1.f);
}
);

const std::string vol_frg_shader=To_String(
~include version;
uniform vec2 screen_size;
uniform float dx;
uniform vec3 org;
uniform vec3 one_over_length;
uniform sampler2D tex2d;
uniform sampler3D tex3d;
uniform sampler1D tex1d;
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
	vec2 screen_coord=vec2(gl_FragCoord.x-.5f,gl_FragCoord.y-.5f);
	vec2 tex_coord=screen_coord/screen_size;
	vec4 end=texture(tex2d,tex_coord);
	vec3 ray=end.xyz-vtx_frg_pos;
	int n=int(ray.length()/dx)+1;
	vec4 ac=vec4(0.f,0.f,0.f,0.f);
	for(int i=0;i<=n;i++){
		vec3 pos=vtx_frg_pos+float(i)/float(n)*ray;
		vec3 coord=(pos-org)*one_over_length;
		float val=texture(tex3d,coord).r;
		vec4 c=texture(tex1d,val);
		c.rgb=(1.f-c.rgb);
		c.a=pow(c.a*.1f,1.f);
		ac.rgb=(1-ac.a)*c.rgb*c.a+ac.rgb;
		ac.a=(1-ac.a)*c.a+ac.a;
		if(ac.a>.95f)break;}
	ac.r=min(ac.r,1.f);ac.g=min(ac.g,1.f);ac.b=min(ac.b,1.f);
	frag_color=ac;
}
);

const std::string vec_vol_frg_shader=To_String(
~include version;
uniform vec2 screen_size;
uniform float dx;
uniform int tex_dim;
uniform vec3 org;
uniform vec3 one_over_length;
uniform vec4 pl_pos;
uniform sampler2D tex2d;
uniform sampler3D tex3d;
uniform sampler1D tex1d[4];
in vec3 vtx_frg_pos;
out vec4 frag_color;
void main()
{
	vec2 screen_coord=vec2(gl_FragCoord.x-.5f,gl_FragCoord.y-.5f);
	vec2 tex_coord=screen_coord/screen_size;
	vec4 end=texture(tex2d,tex_coord);
	vec3 ray=end.xyz-vtx_frg_pos;
	int n=int(ray.length()/dx)+1;
	vec4 ac=vec4(0.f,0.f,0.f,0.f);
	for(int i=0;i<=n;i++){
		vec3 pos=vtx_frg_pos+float(i)/float(n)*ray;
		vec3 coord=(pos-org)*one_over_length;
		vec4 val=texture(tex3d,coord);
		vec4 c=vec4(0.f,0.f,0.f,0.f);
		for(int j=0;j<tex_dim;j++){
			vec4 c0=texture(tex1d[j],val[j]);
			c.rgb+=c0.rgb*c0.a;c.a+=c0.a;}
		c.rgb/=(c.a+.000001f);c.a*=.25f;

		ac.rgb=(1-ac.a)*c.rgb*c.a+ac.rgb;
		ac.a=(1-ac.a)*c.a+ac.a;
		if(ac.a>.95f)break;}
	ac.r=min(ac.r,1.f);ac.g=min(ac.g,1.f);ac.b=min(ac.b,1.f);
	frag_color=ac;

	/*simple pl*/
	/*vec3 coord=(vtx_frg_pos-org)*one_over_length;
	vec3 normal=vec3(0.f);
	float eps=.001f;
	for(int i=0;i<3;i++){
		if(coord[i]<eps){normal[i]=-1.f;}
		else if(coord[i]>1.f-eps)normal[i]=1.f;}
	normal=normalize(normal);
	vec3 lt_dir=normalize(pl_pos.xyz-vtx_frg_pos);
	float dif_coef=max(dot(normal,lt_dir),0.f);
	float amb=.3;
	float coef=dif_coef+amb;
	vec4 dif=coef*ac;
	frag_color=dif;*/
}
);
}
#endif
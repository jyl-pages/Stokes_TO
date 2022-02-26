//////////////////////////////////////////////////////////////////////////
// Opengl color mapper
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "OpenGLColorMapper.h"

//////////////////////////////////////////////////////////////////////////
////OpenGLColorMapper
//////////////////////////////////////////////////////////////////////////

OpenGLColorMapper::OpenGLColorMapper()
{	
	density_color_max=OpenGLColor::White();
	Resize(5);
	Set_Color_Mode((real)0,(real)1,ColorType::Jet,0);
	Set_Color_Mode((real)0,(real)1,ColorType::Hot,1);
	Set_Color_Mode((real)0,(real)1,ColorType::Den,2);
	Set_Color_Mode((real)0,(real)1,ColorType::Mat,3);
	Set_Color_Mode((real)0,(real)1,ColorType::MC,4);
	current_mode_index=0;

	multi_comp_colors.clear();
	multi_comp_colors.push_back(OpenGLColor(.0f,1.f,1.f,1.f));
	multi_comp_colors.push_back(OpenGLColor(1.f,1.f,.0f,1.f));
	multi_comp_colors.push_back(OpenGLColor(1.f,.0f,.0f,1.f));
	multi_comp_colors.push_back(OpenGLColor(.0f,.0f,1.f,1.f));
	multi_comp_colors.push_back(OpenGLColor(1.f,.0f,1.f,1.f));
	multi_comp_colors.push_back(OpenGLColor(.0f,1.f,.0f,1.f));
}

void OpenGLColorMapper::Set_Matlab_Jet(real v_min,real v_max,const int i)
{Resize(i+1);Valid_Min_And_Max(v_min,v_max);Set_Color(OpenGLColorRamp::Matlab_Jet((real)v_min,(real)v_max),ColorType::Jet,i);}

void OpenGLColorMapper::Set_Matlab_Hot(real v_min,real v_max,const int i)
{Resize(i+1);Valid_Min_And_Max(v_min,v_max);Set_Color(OpenGLColorRamp::Matlab_Hot((real)v_min,(real)v_max),ColorType::Hot,i);}

void OpenGLColorMapper::Set_Density_Base(real v_min,real v_max,const int i,const ColorType type)	/*type==ColorType::Den,::Mat,::MC*/
{
	OpenGLColor color_min=OpenGLColor::Transparent();
	//OpenGLColor color_min=OpenGLColor(1.f,1.f,1.f,0.f);

	if(type==ColorType::Mat)color_min=OpenGLColor(0.f,0.f,0.f,1.f);

	OpenGLColor color_max=density_color_max;
	//OpenGLColor color_max=OpenGLColor(0.f,0.f,0.f,1.f);
	if(type==ColorType::Mat)color_max=OpenGLColor(1.f,1.f,1.f,1.f);
	else if(type==ColorType::MC)color_max=OpenGLColor(1.f,0.f,0.f,0.f);

	Valid_Min_And_Max(v_min,v_max);OpenGLColorRamp* c=new OpenGLColorRamp();
	c->Add_Color(v_min-(real)1e10,color_min);c->Add_Color(v_min+(real)1e-10,color_min);c->Add_Color(v_max-(real)1e-10,color_max);
	c->Add_Color(v_max,color_max);c->Add_Color(v_max+(real)1e10,color_max);
	Resize(i+1);Set_Color(c,type,i);
}

void OpenGLColorMapper::Set_Density(real v_min,real v_max,const int i/*=0*/)
{Set_Density_Base(v_min,v_max,i,ColorType::Den);}

void OpenGLColorMapper::Set_Material(real v_min,real v_max,const int i/*=0*/)
{Set_Density_Base(v_min,v_max,i,ColorType::Mat);}

void OpenGLColorMapper::Set_Material_Multi_Component(real v_min,real v_max,const int i/*=0*/)
{Set_Density_Base(v_min,v_max,i,ColorType::MC);}

void OpenGLColorMapper::Set_Color_Mode(real v_min,real v_max,const ColorType type,const int i)
{
	current_mode_index=0;
	real epsilon=(v_max-v_min)*(real).1;
	range[i]=Vector2(v_min-epsilon,v_max+epsilon);
	switch(type){
	case ColorType::Jet:{Set_Matlab_Jet(v_min-epsilon,v_max+epsilon,i);}break;
	case ColorType::Hot:{Set_Matlab_Hot(v_min,v_max,i);}break;
	case ColorType::Den:{Set_Density((real).0,(real)1,i);}break;
	case ColorType::Mat:{Set_Material(v_min,v_max,i);}break;
	case ColorType::MC:{Set_Material_Multi_Component((real)0,(real)1,i);}break;default:break;}
}
	
OpenGLColor OpenGLColorMapper::Color(const real v,const int i) const 
{OpenGLColor c=colors[i]->Lookup(v);/*if(alpha[i]<1.f)c.rgba[3]=alpha[i];*/return c;}

void OpenGLColorMapper::Resize(const int size)
{if((int)colors.size()<size){colors.resize(size);types.resize(size);alpha.resize(size);range.resize(size);}}

OpenGLColor OpenGLColorMapper::Multi_Comp_Color(const Array<real>& v) const
{OpenGLColor c=OpenGLColor(0.f,0.f,0.f,0.f);for(int i=0;i<(int)v.size();i++){float c_i=Color(v[i],4).rgba[0];c=c+c_i*multi_comp_colors[i]*(float)one_third;}
c.rgba[3]=alpha[0];/*c.rgba[3]=(float)one_third*(c.rgba[0]+c.rgba[1]+c.rgba[2]);*/return c;}

void OpenGLColorMapper::Set_Alpha(const GLfloat a,const int i){alpha[i]=a;}
	
void OpenGLColorMapper::Set_Alpha_All(const GLfloat _a){for(auto& a:alpha)a=_a;}

void OpenGLColorMapper::Valid_Min_And_Max(real& v_min,real& v_max){if(v_min>=v_max){v_max=v_min+(real)1e-3;}}

void OpenGLColorMapper::Set_Color(OpenGLColorRamp* color_ramp,ColorType type,const int i,const GLfloat a)
{colors[i].reset(color_ramp);types[i]=type;alpha[i]=a;}

//////////////////////////////////////////////////////////////////////////
//// OpenGLColorTransfer
//////////////////////////////////////////////////////////////////////////

OpenGLColorTransfer::OpenGLColorTransfer(const Array<ArrayF<GLfloat,5> >& _maps){Initialize(_maps);}

void OpenGLColorTransfer::Initialize(const Array<ArrayF<GLfloat,5> >& _maps)
{
	maps=_maps;for(auto& m:maps){color_ramp.Add_Color(m[0],OpenGLColor(m[1],m[2],m[3],m[4]));}
	int size=256;transfer.resize(size*4);std::fill(transfer.begin(),transfer.end(),(ushort)0);
	for(auto i=0;i<size;i++){real c=(real)i/(real)(size-1);OpenGLColor color=color_ramp.Lookup(c);
		for(int j=0;j<4;j++){transfer[i*4+j]=(ushort)((real)0xffff*color.rgba[j]);}}
}
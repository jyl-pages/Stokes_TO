//////////////////////////////////////////////////////////////////////////
// Opengl color mapper
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLColor_h__
#define __OpenGLColor_h__
#include "Constants.h"
#include "OpenGLColorRamp.h"
#include "OpenGLCommon.h"

class OpenGLColorMapper
{public:
	////color modes
	Array<std::shared_ptr<OpenGLColorRamp> > colors;
	Array<ColorType> types;
	Array<GLfloat> alpha;
	Array<Vector2> range;

	int current_mode_index;
	OpenGLColor density_color_max;
	Array<OpenGLColor> multi_comp_colors;

	OpenGLColorMapper();

	void Set_Matlab_Jet(real v_min,real v_max,const int i=0);
	void Set_Matlab_Hot(real v_min,real v_max,const int i=0);
	void Set_Density(real v_min,real v_max,const int i=0);
	void Set_Material(real v_min,real v_max,const int i=0);
	void Set_Material_Multi_Component(real v_min,real v_max,const int i=0);

	void Set_Density_Base(real v_min=(real)1e-5,real v_max=(real)1,const int i=0,const ColorType type=ColorType::Den);	/*type==Density,material,MultiComp*/
	void Set_Color_Mode(real v_min,real v_max,const ColorType type,const int i=0);
	OpenGLColor Color(const real v,const int i=0) const;
	void Resize(const int size);
	OpenGLColor Multi_Comp_Color(const Array<real>& v) const;
	void Set_Alpha(const GLfloat a,const int i=0);
	void Set_Alpha_All(const GLfloat _a);

	template<class T_VAL> void Set_Matlab_Jet(const Array<T_VAL>& v,const int i=0)
	{T_VAL v_min,v_max;Min_And_Max(v,v_min,v_max);Set_Matlab_Jet((real)v_min,(real)v_max,i);}

	template<class T_VAL> void Set_Matlab_Hot(const Array<T_VAL>& v,const int i=0)
	{T_VAL v_min,v_max;Min_And_Max(v,v_min,v_max);Set_Matlab_Hot((real)v_min,(real)v_max,i);}

	template<class T_VAL> void Set_Density(const Array<T_VAL>& v,const int i=0)
	{T_VAL v_min,v_max;Min_And_Max(v,v_min,v_max);Set_Density(v_min,v_max,i);}

	template<class T_VAL> void Set_Color_Mode(const Array<T_VAL>& v,const ColorType type,const int i=0)
	{T_VAL v_min,v_max;Min_And_Max(v,v_min,v_max);Set_Color_Mode(v_min,v_max,type,i);}

	template<class T_VAL> void Reset_Color_Mode(const Array<T_VAL>& v,const int i=0)
	{T_VAL v_min,v_max;Min_And_Max(v,v_min,v_max);Reset_Color_Mode((real)v_min,(real)v_max,i);}
	void Reset_Color_Mode(real v_min,real v_max,const int i=0)
	{Set_Color_Mode(v_min,v_max,types[i],i);}

	template<class T_VAL> void Update(const Array<T_VAL>& v,const int i=0,const ColorType type=ColorType::Den)
	{
		switch(types[i]){
		case ColorType::Jet:Set_Matlab_Jet(v,i);break;
		case ColorType::Hot:Set_Matlab_Hot(v,i);break;
		case ColorType::Den:Set_Density(v,i);break;
		case ColorType::Mat:Set_Material(v,i);break;
		case ColorType::MC:Set_Material_Multi_Component(v,i);break;}
	}

protected:
	void Valid_Min_And_Max(real& v_min,real& v_max);
	void Set_Color(OpenGLColorRamp* color_ramp,ColorType type,const int i=0,const GLfloat a=(GLfloat)0);

	template<class T_VAL> void Offset(T_VAL& v_min,T_VAL& v_max,const T_VAL offset){v_min-=offset;v_max+=offset;}
};

class OpenGLColorTransfer
{
public:
	Array<ArrayF<GLfloat,5> > maps;
	OpenGLColorRamp color_ramp;
	Array<ushort> transfer;

	OpenGLColorTransfer(){}
	OpenGLColorTransfer(const Array<ArrayF<GLfloat,5> >& _maps);
	void Initialize(const Array<ArrayF<GLfloat,5> >& _maps);
};
#endif
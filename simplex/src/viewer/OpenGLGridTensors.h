 //////////////////////////////////////////////////////////////////////////
// Opengl grid tensors
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLGridTensors_h__
#define __OpenGLGridTensors_h__
#include "OpenGLGridField.h"

template<class T> class OpenGLGridTensors : public OpenGLGridField<Matrix<T,3>,1>
{typedef OpenGLGridField<Matrix<T,3>,1> Base;
public:
	using Base::use_subgrid_view;using Base::subgrid_cell_start;using Base::subgrid_cell_end;	
	using Base::name;using Base::field_draw_type;using Base::field_val;using Base::mac_grid;
	using Base::Update_Data_To_Render_Grid_Lines;using Base::Display_Lines;
	OpenGLGridTensors(){name="grid_tensors";field_draw_type=FieldDrawType::Line;}
	
	virtual void Update_Data_To_Render()
	{
		if(use_subgrid_view)Update_Data_To_Render_Grid_Lines(mac_grid,field_val,&subgrid_cell_start,&subgrid_cell_end);
		else Update_Data_To_Render_Grid_Lines(mac_grid,field_val);
	}

	virtual void Display() const {Display_Lines();}
};
#endif
//////////////////////////////////////////////////////////////////////////
// Opengl grid height field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLGridHeightField_h__
#define __OpenGLGridHeightField_h__
#include "OpenGLGridField.h"

////TODO: cell color
template<class T> class OpenGLGridHeightField : public OpenGLGridField<T,4>
{typedef OpenGLGridField<T,4> Base;
public:
	using Base::name;using Base::polygon_mode;using Base::field_draw_type;using Base::line_norm;using Base::color;
	using Base::mac_grid;using Base::scale;using Base::field_val;using Base::data_idx;using Base::visible;using Base::dim;
	using Base::opengl_vertices;using Base::opengl_elements;using Base::shader_programs;using Base::vao;using Base::ele_size;
	using Base::Set_OpenGL_Vertex_Attribute;using Base::Set_OpenGL_Vertices;using Base::Set_OpenGL_Elements;using Base::Update_Polygon_Mode;
	using Base::Update_Data_To_Render_Pre;using Base::Update_Data_To_Render_Post;using Base::Scalar;
	using Base::use_alpha_blend; using Base::data; using Base::color_mapper; using Base::Update_Scalar_Range; using Base::Enable_Alpha_Blend; using Base::Disable_Alpha_Blend;

	OpenGLGridHeightField():Base()
	{name="grid_height_field";polygon_mode=PolygonMode::Wireframe;field_draw_type=FieldDrawType::Line;line_norm=(real)1;color=OpenGLColor::Blue();}

	virtual void Update_Color_Mapper()
	{
		if(use_alpha_blend){
			T v_min=(T)0,v_max=(T)0;Update_Scalar_Range(field_val[data_idx]->array,v_min,v_max);
			if(data[data_idx].color_type==ColorType::Den){
				color_mapper->density_color_max=OpenGLColor(data[data_idx].color);}
			color_mapper->Reset_Color_Mode(v_min,v_max,(int)data[data_idx].color_type);}
		////No alpha update
	}

	virtual void Update_Data_To_Render()
	{
		if(!Check_Valid_Height_Field())return;
		if(!Update_Data_To_Render_Pre())return;

		int axis=2;
		iterate_node_d(iter,mac_grid.grid,3){
			const Vector3i& node=iter.Coord();Vector3 pos=mac_grid.grid.Node(node);
			OpenGLColor color=color_mapper->Color(Scalar((*field_val[data_idx])(node)),(int)data[data_idx].color_type);
			pos[axis]=scale*Scalar((*field_val[data_idx])(node));
			OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats

		Base::template Update_OpenGL_Elements_For_Node_Data<2>();

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();

		Update_Data_To_Render_Post();
	}

    virtual void Display() const
    {
		if(!visible)return;
		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		shader->Set_Uniform_Vec4f("color",color.rgba);
		glBindVertexArray(vao);

		if(use_alpha_blend)Enable_Alpha_Blend();
		glEnable(GL_PRIMITIVE_RESTART);
		glPrimitiveRestartIndex((GLuint)UINT_MAX);
		glDrawElements(GL_TRIANGLE_STRIP,ele_size,GL_UNSIGNED_INT,0);
		glDisable(GL_PRIMITIVE_RESTART);
		if(use_alpha_blend)Disable_Alpha_Blend();
		shader->End();}
    }

protected:
	bool Check_Valid_Height_Field() const {bool valid=/*dim==2&&*/mac_grid.grid.node_counts[2]==1;
	if(!valid)std::cout<<"Invalid dim for OpenGLGridHeightField"<<std::endl;return valid;}
};
#endif

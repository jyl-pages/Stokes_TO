//////////////////////////////////////////////////////////////////////////
// Opengl grid vectors
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLGridField_h__
#define __OpenGLGridField_h__
#include "OpenGLField.h"

////TODO: cell color
template<class T_VAL,int data_size>
class OpenGLGridField : public OpenGLField<Field<T_VAL,3>,data_size>
{typedef OpenGLField<Field<T_VAL,3>,data_size> Base;
public:
	using Base::name;using Base::visible;using Base::data;using Base::data_idx;using Base::dim;using Base::color_mapper;
	using Base::field_val;using Base::opengl_vertices;using Base::shader_programs;using Base::color;using Base::vao;using Base::ele_size;
	using Base::output_dir;using Base::data_refreshed;using Base::opengl_elements;using Base::verbose;
	using Base::Update_Data_To_Render_Pre;using Base::Update_Data_To_Render_Post;using Base::Set_Data_Refreshed;
	using Base::Set_OpenGL_Vertices;using Base::Set_OpenGL_Vertex_Attribute;using Base::Set_OpenGL_Elements;
	using Base::Update_Polygon_Mode;using Base::Read_Grid_Field;using Base::Object_File_Name;using Base::Scalar;
	using Base::Update_Scalar_Range;using Base::Enable_Alpha_Blend;using Base::Disable_Alpha_Blend;
	using Base::use_const_color_map; using Base::color_map_inited;
	using Base::Update_Data_To_Render_Grid_Lines;using Base::Display_Lines;
	MacGrid<3> mac_grid;
	bool use_alpha_blend=true;
	bool use_subgrid_view=false;
	Vector3i subgrid_cell_start=Vector3i::Zero();
	Vector3i subgrid_cell_end=Vector3i::Zero();

	OpenGLGridField():Base(){name="grid_field";}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;
		switch(data[data_idx].store_type){
		case StoreType::Node:Update_Node_Data_To_Render();break;
		case StoreType::Cell:
			{if(dim==2)Update_Cell_Data_To_Render<2>();
			else Update_Cell_Data_To_Render<3>();}break;}
		Update_Data_To_Render_Post();
	}

	virtual void Update_Color_Mapper()
	{
		if(use_alpha_blend){
			if(data[data_idx].color_type==ColorType::Den){
				//color_mapper->density_color_max=OpenGLColor(data[data_idx].color);	
				////multiple data_idx is being deprecated, let's assume we always want to use single piece of data
				color_mapper->density_color_max=OpenGLColor(color);}

			if(!use_const_color_map||!color_map_inited){
				real v_min=(real)0,v_max=(real)0;Update_Scalar_Range(field_val[data_idx]->array,v_min,v_max);
				color_mapper->Reset_Color_Mode(v_min,v_max,(int)data[data_idx].color_type);
				color_map_inited=true;}}
		////No alpha update
	}

	void Update_Node_Data_To_Render()
	{
		if(dim==2){
			if(!use_subgrid_view)Create_Node_Data<2>();
			else Create_Node_Data_With_Subgrid<2>();
			Update_OpenGL_Elements_For_Node_Data<2>();}
		else{
			if(!use_subgrid_view)Create_Node_Data<3>();
			else Create_Node_Data_With_Subgrid<3>();
			Update_OpenGL_Elements_For_Node_Data<3>();}
		
		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();	
	}

	template<int dim> void Update_Cell_Data_To_Render()
	{
		if(!use_subgrid_view)Create_Cell_Data<dim>();
		else Create_Cell_Data_With_Subgrid<dim>();

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		Set_OpenGL_Elements();		
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
	
	////Grid helper
	Grid<3> Grid_2D(const Grid<3>& grid)	////input grid, cell_counts[2]=0
	{Grid<3> g3=grid;g3.cell_counts[2]=1;return g3;}

	virtual void Refresh(const int frame)
	{
		if(frame==0){std::string file_name=Object_File_Name(output_dir,frame,name);
			if(File::File_Exists(file_name)){
				Grid<3> lat;File::Read_Binary_From_File(file_name,lat);mac_grid.Initialize(lat);
				if(mac_grid.grid.cell_counts[2]==0){dim=2;mac_grid.grid=Grid_2D(mac_grid.grid);}
				Set_Data_Refreshed();
				if(verbose)std::cout<<"Read file "<<file_name<<std::endl;}}

		for(size_type i=0;i<data.size();i++){std::string file_name=Object_File_Name(output_dir,frame,data[i].name);
			if(File::File_Exists(file_name)){
				Read_Grid_Field(file_name,mac_grid.grid.cell_counts,mac_grid.grid.node_counts,field_val[i],data[i].store_type);
				Set_Data_Refreshed();}}
	}

	void Use_Subgrid_View(bool use=true){use_subgrid_view=use;}

	void Set_Subgrid(const Vector3i& start,const Vector3i& end)
	{subgrid_cell_start=start;subgrid_cell_end=end;}

	void Set_Subgrid_Slice(const int axis,const int start,const int end=-1)
	{
		Use_Subgrid_View();
		subgrid_cell_start=Vector3i::Zero()+Vector3i::Unit(axis)*start;
		subgrid_cell_end=mac_grid.grid.cell_counts-Vector3i::Ones();
		subgrid_cell_end[axis]=(end==-1)?subgrid_cell_start[axis]:end;
	}

protected:
	void Update_OpenGL_Elements(ArrayF2P<int,3>& nodes)
	{OpenGL_Cube(nodes,opengl_elements);/*cube, 18 uint*/}

	void Update_OpenGL_Elements(ArrayF2P<int,2>& nodes)
	{OpenGL_Quad(nodes,opengl_elements);/*cube, 5 uint*/}

	template<int dim> void Update_OpenGL_Elements_For_Node_Data()
	{
		iterate_cell_d(iter,mac_grid.grid,3){const Vector3i& cell=iter.Coord();
			ArrayF2P<int,dim> nodes;for(int i=0;i<(int)nodes.size();i++){nodes[i]=mac_grid.grid.Node_Index(mac_grid.grid.Cell_Incident_Node(cell,dim==2?(i*2):i));}
			Update_OpenGL_Elements(nodes);/*cube, 18 uint*/}
	}

	template<int dim> void Create_Node_Data()
	{
		iterate_node_d(iter,mac_grid.grid,3){
			const Vector3i& node=iter.Coord();Vector3 pos=mac_grid.grid.Node(node);
			OpenGLColor color=color_mapper->Color(Scalar((*field_val[data_idx])(node)),(int)data[data_idx].color_type);
			OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats	
	}

	template<int dim> void Create_Node_Data_With_Subgrid()
	{
		Vector3i node_start=subgrid_cell_start;
		Vector3i node_end=subgrid_cell_end+Vector3i::Ones();
		iterate_range_d(iter,node_start,node_end,3){
			const Vector3i& node=iter.Coord();Vector3 pos=mac_grid.grid.Node(node);
			OpenGLColor color=color_mapper->Color(Scalar((*field_val[data_idx])(node)),(int)data[data_idx].color_type);
			OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats		
	}

	template<int dim> void Create_Cell_Data()
	{
		int c=0;
		iterate_cell_d(iter,mac_grid.grid,3){const Vector3i cell=iter.Coord();
			OpenGLColor color=color_mapper->Color(Scalar((*field_val[data_idx])(cell)),(int)data[data_idx].color_type);
			ArrayF2P<int,dim> nodes;for(int i=0;i<nodes.size();i++){nodes[i]=c++;
				Vector3i inc_node=mac_grid.grid.Cell_Incident_Node(cell,dim==2?(i*2):i);
				Vector3 pos=mac_grid.grid.Node(inc_node);
				OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats	
			Update_OpenGL_Elements(nodes);}	
	}

	template<int dim> void Create_Cell_Data_With_Subgrid()
	{
		int c=0;
		iterate_range_d(iter,subgrid_cell_start,subgrid_cell_end,3){const Vector3i cell=iter.Coord();
			OpenGLColor color=color_mapper->Color(Scalar((*field_val[data_idx])(cell)),(int)data[data_idx].color_type);
			ArrayF2P<int,dim> nodes;for(int i=0;i<nodes.size();i++){nodes[i]=c++;
				Vector3i inc_node=mac_grid.grid.Cell_Incident_Node(cell,dim==2?(i*2):i);
				Vector3 pos=mac_grid.grid.Node(inc_node);
				OpenGL_Vertex4_And_Color4(pos,color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats	
			Update_OpenGL_Elements(nodes);}		
	}
};
#endif

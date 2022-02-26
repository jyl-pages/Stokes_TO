//////////////////////////////////////////////////////////////////////////
// Opengl grid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLGrid_h__
#define __OpenGLGrid_h__
#include "MacGrid.h"
#include "OpenGLObject.h"
#include "OpenGLShaderProgram.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLPrimitives.h"
#include "OpenGLUbos.h"

class OpenGLGrid : public OpenGLObject
{typedef OpenGLObject Base;
public:
	MacGrid<3> mac_grid;
	bool use_subgrid_view=false;
	Vector3i subgrid_cell_start=Vector3i::Zero();
	Vector3i subgrid_cell_end=Vector3i::Zero();

	OpenGLGrid(){color=OpenGLColor::Gray(.5f);name="grid";}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));
	}

	virtual void Update_Data_To_Render()
	{
		if(!Update_Data_To_Render_Pre())return;

		if(!use_subgrid_view)Create_Grid_Lines();
		else Create_Grid_Lines_With_Subgrid();

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,3);

		Update_Data_To_Render_Post();
	}

    virtual void Display() const
    {
    	using namespace OpenGLUbos;
		if(!visible)return;
		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		Bind_Uniform_Block_To_Ubo(shader,"camera");
		shader->Set_Uniform_Vec4f("color",color.rgba);
		glBindVertexArray(vao);
		glDrawArrays(GL_LINES,0,vtx_size/3);
		shader->End();}
    }
	
	virtual void Refresh(const int frame)
	{
		std::string file_name=output_dir+"/"+std::to_string(frame)+"/"+name;
		if(File::File_Exists(file_name)){
			File::Read_Binary_From_File(file_name,mac_grid.grid);mac_grid.Initialize(mac_grid.grid);Set_Data_Refreshed();
			if(verbose)std::cout<<"Read file "<<file_name<<std::endl;}
	}

	void Create_Grid_Lines()
	{
		Grid<3>& grid=mac_grid.grid;
		////x lines
		if(grid.node_counts[0]>1){for(int i=0;i<grid.node_counts[1];i++)for(int j=0;j<grid.node_counts[2];j++){
			Vector3 start=grid.Node(Vector3i(0,i,j));
			Vector3 end=grid.Node(Vector3i(grid.node_counts[0]-1,i,j));
			OpenGL_Vertex(start,opengl_vertices);
			OpenGL_Vertex(end,opengl_vertices);}}
		////y lines
		if(grid.node_counts[1]>1){for(int i=0;i<grid.node_counts[2];i++)for(int j=0;j<grid.node_counts[0];j++){
			Vector3 start=grid.Node(Vector3i(j,0,i));
			Vector3 end=grid.Node(Vector3i(j,grid.node_counts[1]-1,i));
			OpenGL_Vertex(start,opengl_vertices);
			OpenGL_Vertex(end,opengl_vertices);}}
		////z lines
		if(grid.node_counts[2]>1){for(int i=0;i<grid.node_counts[0];i++)for(int j=0;j<grid.node_counts[1];j++){
			Vector3 start=grid.Node(Vector3i(i,j,0));
			Vector3 end=grid.Node(Vector3i(i,j,grid.node_counts[2]-1));
			OpenGL_Vertex(start,opengl_vertices);
			OpenGL_Vertex(end,opengl_vertices);}}	
	}

	void Create_Grid_Lines_With_Subgrid()
	{
		Grid<3>& grid=mac_grid.grid;
		////x lines
		if(subgrid_cell_end[0]-subgrid_cell_start[0]>=0){
			for(int i=subgrid_cell_start[1];i<=subgrid_cell_end[1]+1;i++)
				for(int j=subgrid_cell_start[2];j<=subgrid_cell_end[2]+1;j++){
					Vector3 start=grid.Node(Vector3i(subgrid_cell_start[0],i,j));
					Vector3 end=grid.Node(Vector3i(subgrid_cell_end[0]+1,i,j));
					OpenGL_Vertex(start,opengl_vertices);
					OpenGL_Vertex(end,opengl_vertices);}}
		////y lines
		if(subgrid_cell_end[1]-subgrid_cell_start[1]>=0){
			for(int i=subgrid_cell_start[2];i<=subgrid_cell_end[2]+1;i++)
				for(int j=subgrid_cell_start[0];j<=subgrid_cell_end[0]+1;j++){
					Vector3 start=grid.Node(Vector3i(j,subgrid_cell_start[1],i));
					Vector3 end=grid.Node(Vector3i(j,subgrid_cell_end[1]+1,i));
					OpenGL_Vertex(start,opengl_vertices);
					OpenGL_Vertex(end,opengl_vertices);}}
		////z lines
		if(subgrid_cell_end[2]-subgrid_cell_start[2]>=0){
			for(int i=subgrid_cell_start[0];i<=subgrid_cell_end[0]+1;i++)
				for(int j=subgrid_cell_start[1];j<=subgrid_cell_end[1]+1;j++){
					Vector3 start=grid.Node(Vector3i(i,j,subgrid_cell_start[2]));
					Vector3 end=grid.Node(Vector3i(i,j,subgrid_cell_end[2]+1));
					OpenGL_Vertex(start,opengl_vertices);
					OpenGL_Vertex(end,opengl_vertices);}}		
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
};
#endif
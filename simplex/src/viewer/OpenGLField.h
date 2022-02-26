//////////////////////////////////////////////////////////////////////////
// Opengl field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLField_h__
#define __OpenGLField_h__
#include "File.h"
#include "Field.h"
#include "AuxFunc.h"
#include "MacGrid.h"
#include "Interpolation.h"
#include "OpenGLObject.h"
#include "OpenGLMesh.h"
#include "OpenGLColorMapper.h"

template<class T_FIELD,int data_size=4>
class OpenGLField : public OpenGLObject
{typedef OpenGLObject Base;
public:
	ArrayF<std::shared_ptr<T_FIELD>,data_size> field_val;
	FieldDrawType field_draw_type;
	TensorFieldDrawType tensor_field_draw_type;
	bool draw_arrow=true;
	bool use_const_color_map=false;
	bool color_map_inited=false;

	OpenGLField()
	{color=OpenGLColor::Red(.5f);name="field";polygon_mode=PolygonMode::Fill;field_draw_type=FieldDrawType::Color;
	tensor_field_draw_type=TensorFieldDrawType::Eigenvector;line_norm=(real).2;
	for(int i=0;i<data_size;i++){field_val[i]=std::make_shared<T_FIELD>();}}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("vcolor"));
	}

	virtual void Update_Color_Mapper()
	{
		if(data[data_idx].color_type==ColorType::Den){
			color_mapper->density_color_max=OpenGLColor(data[data_idx].color);}

		if(!use_const_color_map||!color_map_inited){
			real v_min=(real)0,v_max=(real)0;Update_Scalar_Range(field_val[data_idx]->array,v_min,v_max);
			color_mapper->Reset_Color_Mode(v_min,v_max,(int)data[data_idx].color_type);
			color_map_inited=true;}

		color_mapper->Set_Alpha_All(alpha);
	}
	
	////Read a node field; interpolate if the data is a cell field
	template<class T_VAL> void Read_Grid_Field(const std::string& file_name,const Vector3i& cell_counts,const Vector3i& node_counts,
		std::shared_ptr<Field<T_VAL,3> >& val,const StoreType store_type)
	{
		std::shared_ptr<Field<T_VAL,3> > read_field=std::make_shared<Field<T_VAL,3> >();
		File::Read_Binary_From_File(file_name,*read_field);

		StoreType read_type=(StoreType)Field_Type<3>(cell_counts,node_counts,read_field->counts);
		if(read_type==StoreType::None){
			std::cerr<<"Error: [OpenGLField] grid field size do not match, read_counts: "<<read_field->counts.transpose()
				<<", cell_counts: "<<cell_counts.transpose()<<", node_counts: "<<node_counts.transpose()<<std::endl;return;}

		switch(store_type){
		case StoreType::Node:{
			switch(read_type){
			case StoreType::Node:{val=read_field;}break;
			case StoreType::Cell:{val->Resize(node_counts);Interpolate_Grid_Cells_To_Nodes(cell_counts,node_counts,*read_field,*val);}break;}
		}break;
		case StoreType::Cell:{
			switch(read_type){
			case StoreType::Node:{val->Resize(cell_counts);Interpolate_Grid_Nodes_To_Cells(cell_counts,node_counts,*read_field,*val);}break;
			case StoreType::Cell:{val=read_field;}break;}
		}break;}

		if(verbose)std::cout<<"Read grid field: read type: "<<(int)read_type<<", store type: "<<(int)store_type<<", file: "<<file_name<<std::endl;
	}

	////read node field, if element field, do interpolation
	template<class T_VAL,int e> void Read_Mesh_Field(const std::string& file_name,const int element_num,const int vertex_num,
		const Array<Vector<int,e> >& elements,std::shared_ptr<Field<T_VAL,1> >& val,const StoreType store_type)
	{
		std::shared_ptr<Field<T_VAL,1> > read_field=std::make_shared<Field<T_VAL,1> >();
		File::Read_Binary_From_File(file_name,*read_field);

		Vector1i cell_counts=Vec1i(element_num);Vector1i node_counts=Vec1i(vertex_num);
		StoreType read_type=(StoreType)Field_Type<1>(cell_counts,node_counts,read_field->counts);
		if(read_type==StoreType::None){std::cerr<<"Error: [OpenGLField] Mesh field size do not match, read_counts: "<<read_field->counts.transpose()
			<<", cell_counts: "<<cell_counts.transpose()<<", node_counts: "<<node_counts.transpose()<<std::endl;
			return;}

		switch(store_type){
		case StoreType::Node:{
			switch(read_type){
			case StoreType::Node:{val=read_field;}break;
			case StoreType::Cell:{val->Resize(node_counts);Interpolate_Mesh_Cells_To_Nodes(cell_counts[0],node_counts[0],elements,*read_field,*val);}break;}
		}break;
		case StoreType::Cell:{
			switch(read_type){
			case StoreType::Node:{val->Resize(cell_counts);Interpolate_Mesh_Nodes_To_Cells(cell_counts[0],node_counts[0],elements,*read_field,*val);}break;
			case StoreType::Cell:{val=read_field;}break;}
		}break;}

		if(verbose)
			std::cout<<"Read mesh field: read type: "<<(int)read_type<<", store type: "<<(int)store_type<<", file: "<<file_name<<std::endl;
	}

	virtual void Normalize_Data()
	{if(field_draw_type==FieldDrawType::Line){scale=Normalized_Scale(field_val[0]->array,line_norm);Update_Data_To_Render();}}
	
	////interactive data
	template<class T_VAL,int d> void Copy_Field(const Field<T_VAL,d>& read_field,std::shared_ptr<T_FIELD>& val)
	{
		Vector3i counts_3d;Dim_Conversion(read_field.counts,counts_3d,1);
		val->Resize(counts_3d);for(auto i=0;i<val->array.size();i++)val->array[i]=(real)read_field.array[i];
		Set_Data_Refreshed();
	}

protected:
	////field helper functions
	template<int d> StoreType Field_Type(const Vector<int,d>& cell_counts,const Vector<int,d>& node_counts,const Vector<int,d>& field_counts)
	{
		if(field_counts==cell_counts){return StoreType::Cell;}
		else if(field_counts==node_counts){return StoreType::Node;}
		else return StoreType::None;			
	}

	Grid<3> Dim_Consistent_Grid(const Vector3i& cell_counts,const Vector3i& node_counts)
	{Grid<3> lat(cell_counts);lat.node_counts=node_counts;return lat;}

	template<class T_VAL> void Interpolate_Grid_Cells_To_Nodes(const Vector3i& cell_counts,const Vector3i& node_counts,const Field<T_VAL,3>& read_field,Field<T_VAL,3>& write_field)
	{
		Grid<3> intp_grid=Dim_Consistent_Grid(cell_counts,node_counts);
		Interpolation<3> intp(intp_grid);intp.Interpolate_Cells_To_Nodes(read_field,write_field);	
	}

	template<class T_VAL> void Interpolate_Grid_Nodes_To_Cells(const Vector3i& cell_counts,const Vector3i& node_counts,const Field<T_VAL,3>& read_field,Field<T_VAL,3>& write_field)
	{
		Grid<3> intp_grid=Dim_Consistent_Grid(cell_counts,node_counts);
		Interpolation<3> intp(intp_grid);intp.Interpolate_Nodes_To_Cells(read_field,write_field);	
	}

	template<class T_VAL,int e> void Interpolate_Mesh_Cells_To_Nodes(const int& nc,const int& nv,const Array<Vector<int,e> >& elements,const Field<T_VAL,1>& read_field,Field<T_VAL,1>& write_field)
	{
		const Array<T_VAL>& cell_values=read_field.array;
		Array<T_VAL>& node_values=write_field.array;AuxFunc::Fill(node_values,Zero<T_VAL>());
		Array<real> node_weights;node_weights.resize(nv);AuxFunc::Fill(node_weights,(real)0);
		for(size_type i=0;i<cell_values.size();i++){const Vector<int,e>& vtx=elements[i];
			for(int j=0;j<e;j++){const int v=vtx[j];real c=(real)1;
				node_values[v]=node_values[v]+cell_values[i]*c;node_weights[v]+=c;}}
		for(size_type i=0;i<node_values.size();i++)if(node_weights[i]!=(real)0)node_values[i]=node_values[i]/node_weights[i];
	}

	template<class T_VAL,int e> void Interpolate_Mesh_Nodes_To_Cells(const int& nc,const int& nv,const Array<Vector<int,e> >& elements,const Field<T_VAL,1>& read_field,Field<T_VAL,1>& write_field)
	{
		const Array<T_VAL>& node_values=read_field.array;
		Array<T_VAL>& cell_values=write_field.array;
		for(size_type i=0;i<cell_values.size();i++){const Vector<int,e>& vtx=elements[i];
			cell_values[i]=Zero<T_VAL>();	
			for(int j=0;j<e;j++){cell_values[i]=cell_values[i]+node_values[vtx[j]];}
			cell_values[i]=cell_values[i]/(real)e;}
	}

	void Update_Vertex_Lines_Helper(const Vector3& val,const Vector3& pos)
	{
		OpenGLColor color=color_mapper->Color(Scalar(val));
		ArrayF<Vector3,2> p;
		p[0]=pos;p[1]=p[0]+val*scale;
		for(int i=0;i<2;i++){
			OpenGL_Vertex4_And_Color4(p[i],color.rgba,opengl_vertices);}		////position, 4 floats; color, 4 floats
		if(draw_arrow){
			Vector3 y=Vector3::Unit(1);
			Vector3 axis=y.cross(val).normalized();
			real length=val.norm()*scale;
			real cone_r=std::min(length*(real).04,(real).04);
			Vector3 cone_vec=axis.cross(val).normalized()*cone_r;
			ArrayF<Vector3,2> q;
			q[0]=pos+val*scale*(real).8+cone_vec;
			q[1]=q[0]-(real)2*cone_vec;
			for(int i=0;i<2;i++){
				OpenGL_Vertex4_And_Color4(q[i],color.rgba,opengl_vertices);
				OpenGL_Vertex4_And_Color4(p[1],color.rgba,opengl_vertices);}}
	}

	void Update_Vertex_Lines_Helper(const Matrix3& val,const Vector3& pos)
	{
		ArrayF<Vector3,2> p;OpenGLColor color;
		for(int j=0;j<3;j++){
			switch(tensor_field_draw_type){
			case TensorFieldDrawType::Eigenvector:{
				Sym_Mat3_Eig eig(val);
				real lambda=eig.eigenvalues()[j];
				Vector3 half_v=eig.eigenvectors().col(j)*lambda*scale*(real).5;
				p[0]=pos-half_v;p[1]=pos+half_v;
				color=color_mapper->Color(abs(lambda));
			}break;
			case TensorFieldDrawType::Frame:{
				Vector3 v=val.col(j);
				p[0]=pos;p[1]=pos+v*scale;
				if(j==0) color=OpenGLColor::Red();
				else if(j==1) color=OpenGLColor::Green();
				else color=OpenGLColor::Blue();
			}break;}

			for(int i=0;i<2;i++){
				OpenGL_Vertex4_And_Color4(p[i],color.rgba,opengl_vertices);}}		////position, 4 floats; color, 4 floats
	}

	virtual void Display_Lines() const
    {
		if(!visible)return;
		Update_Polygon_Mode();
		{std::shared_ptr<OpenGLShaderProgram> shader=shader_programs[0];
		shader->Begin();
		glBindVertexArray(vao);
		glDrawArrays(GL_LINES,0,(GLsizei)(vtx_size/8));
		shader->End();}
    }

	template<class T_FIELD2> void Update_Data_To_Render_Grid_Lines(const MacGrid<3>& mac_grid,const T_FIELD2& field_val,const Vector3i* subgrid_cell_start=nullptr,const Vector3i* subgrid_cell_end=nullptr)	////vector=0,tensor=1
	{
		if(!Update_Data_To_Render_Pre())return;
	
		if(subgrid_cell_start!=nullptr&&subgrid_cell_end!=nullptr){
			switch(data[0].store_type){
			case StoreType::Node:{
				Vector3i node_start=*subgrid_cell_start;Vector3i node_end=*subgrid_cell_end+Vector3i::Ones();
				iterate_range_d(iter,node_start,node_end,3){Update_Vertex_Lines_Helper((*field_val[0])(iter.Coord()),mac_grid.grid.Node(iter.Coord()));}}break;
			case StoreType::Cell:{
				iterate_range_d(iter,*subgrid_cell_start,*subgrid_cell_end,3){Update_Vertex_Lines_Helper((*field_val[0])(iter.Coord()),mac_grid.grid.Center(iter.Coord()));}}break;}		
		}
		else{
			switch(data[0].store_type){
			case StoreType::Node:{iterate_node_d(iter,mac_grid.grid,3){Update_Vertex_Lines_Helper((*field_val[0])(iter.Coord()),mac_grid.grid.Node(iter.Coord()));}}break;
			case StoreType::Cell:{iterate_cell_d(iter,mac_grid.grid,3){Update_Vertex_Lines_Helper((*field_val[0])(iter.Coord()),mac_grid.grid.Center(iter.Coord()));}}break;}}

		
		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color
		
		Update_Data_To_Render_Post();
	}

	template<class T_FIELD2,class T_MESH> void Update_Data_To_Render_Mesh_Lines(const OpenGLMesh<T_MESH>* opengl_mesh,const T_FIELD2& field_val)	////vector=0,tensor=1
	{
		if(!Update_Data_To_Render_Pre())return;

		if(opengl_mesh==nullptr){return;}
		auto& mesh=opengl_mesh->mesh;

		switch(data[0].store_type){
		case StoreType::Node:{
			for(size_type i=0;i<mesh.Vertices().size();i++){
				const auto& val=(*field_val[0])(Vec1i((int)i));
				const Vector3& pos=opengl_mesh->Vertex((int)i);
				Update_Vertex_Lines_Helper(val,pos);}
		}break;
		case StoreType::Cell:{
			for(size_type i=0;i<mesh.elements.size();i++){
				const auto& val=(*field_val[0])(Vec1i((int)i));
		        Vector3 pos=Vector3::Zero();
				for(auto j=0;j<mesh.Element_Dim();j++)
					pos+=opengl_mesh->Vertex(mesh.elements[i][j]);
				pos/=(real)mesh.Element_Dim();
                Update_Vertex_Lines_Helper(val,pos);}
		}break;}

		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,4,8,0);	////position
		Set_OpenGL_Vertex_Attribute(1,4,8,4);	////color

		Update_Data_To_Render_Post();
	}
};
#endif

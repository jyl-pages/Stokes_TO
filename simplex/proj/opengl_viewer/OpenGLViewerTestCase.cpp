//////////////////////////////////////////////////////////////////////////
// Opengl viewer test
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#ifdef USE_TINY_OBJ_LOADER
#include "TinyObjLoader.h"
#endif
#include "MeshFunc.h"
#include "AuxFunc.h"
#include "OpenGLViewer.h"
#include "OpenGLAabb.h"
#include "OpenGLGrid.h"
#include "OpenGLGridField.h"
#include "OpenGLGridVectors.h"
#include "OpenGLGridTensors.h"
#include "OpenGLGridHeightField.h"
#include "OpenGLMesh.h"
#include "OpenGLMeshField.h"
#include "OpenGLMeshVectors.h"
#include "OpenGLMeshTensors.h"
#include "OpenGLParticles.h"
#include "OpenGLScreenObjects.h"
#include "OpenGLViewerTestCase.h"

void OpenGLViewerTestCase::Initialize_Data()
{
	Generate_Test_Data();

	Initialize_Mesh_Data();
	//Initialize_Grid_Data();
	//Initialize_Volume_Data();
	//Initialize_2D_Data();

	Add_Object<OpenGLAxes>();
	//Add_Object<OpenGLBar>();
	//OpenGLColorBar* color_bar=new OpenGLColorBar(Vector2(256,256));Add_Object(color_bar);
}

void OpenGLViewerTestCase::Initialize_Mesh_Data()
{
	////Add p
	auto opengl_particles=Add_Object<OpenGLParticles<Particles<3> > >("p");

	////Add a segment mesh
	auto opengl_seg_mesh=Add_Object<OpenGLSegmentMesh>("segment_mesh");

	////Add a triangle mesh
	auto opengl_tri_mesh=Add_Object<OpenGLTriangleMesh>("triangle_mesh.txt");
	Set_Polygon_Mode(opengl_tri_mesh,PolygonMode::Wireframe);

	//auto opengl_bunny=Add_Object<OpenGLTriangleMesh>("bunny.txt");
	//Set_Polygon_Mode(opengl_bunny,PolygonMode::Wireframe);

	//////Add a tri field
	//{auto opengl_mesh_scalar_field=Add_Mesh_Scalar_Field(opengl_tri_mesh,{{"tri_field_scalar"}});
	//////The second parameter is the data initializer list: {{name},{key},{color_type},{store_type},{color}}, see OpenGLObject.h
	//if(opengl_mesh_scalar_field!=nullptr){
	//	opengl_mesh_scalar_field->visible=true;
	//	Bind_Callback_Key('1',&opengl_mesh_scalar_field->Toggle_Draw_Func,"draw tri scalar field");}}

	////Add a tet mesh
	auto opengl_tet_mesh=Add_Object<OpenGLTetrahedronMesh>("tetrahedron_mesh");	////init from file
	Set_Polygon_Mode(opengl_tet_mesh,PolygonMode::Wireframe);

	//////Add tet fields: scalar, vector, tensor; invisible by default
	//{auto opengl_mesh_scalar_field=Add_Mesh_Scalar_Field(opengl_tet_mesh,{{"tet_field_scalar"}});
	//if(opengl_mesh_scalar_field!=nullptr){
	//	opengl_mesh_scalar_field->visible=false;
	//	Bind_Callback_Key('2',&opengl_mesh_scalar_field->Toggle_Draw_Func,"draw tet scalar field");}}

	//{auto opengl_mesh_vector_field=Add_Mesh_Vector_Field(opengl_tet_mesh,{{"tet_field_vector"}});
	//if(opengl_mesh_vector_field!=nullptr){
	//	opengl_mesh_vector_field->visible=false;
	//	Bind_Callback_Key('3',&opengl_mesh_vector_field->Toggle_Draw_Func,"draw tet vector field");}}

	//{auto opengl_mesh_tensor_field=Add_Mesh_Tensor_Field(opengl_tet_mesh,{{"tet_field_tensor","",ColorType::Jet,StoreType::Cell}});
	//if(opengl_mesh_tensor_field!=nullptr){
	//	opengl_mesh_tensor_field->visible=false;
	//	Bind_Callback_Key('4',&opengl_mesh_tensor_field->Toggle_Draw_Func,"draw tet tensor field");}}

	//////Add a quad mesh
	//auto opengl_quad_mesh=Add_Object<OpenGLQuadMesh>("quad_mesh");
}

void OpenGLViewerTestCase::Initialize_Grid_Data()
{
	////Grid and fields
	{auto opengl_grid=Add_Object<OpenGLGrid>("grid");}

	{auto opengl_grid_scalar_field=Add_Grid_Scalar_Field("grid",{{"node_field_scalar","",ColorType::Jet,StoreType::Node}});
	Set_Visibility(opengl_grid_scalar_field,'1',false);
	if(opengl_grid_scalar_field!=nullptr)Bind_Callback_Key('1',&opengl_grid_scalar_field->Toggle_Draw_Func,"draw grid scalar field");}

	{auto opengl_grid_tensor_field=Add_Grid_Tensor_Field("grid",{{"cell_field_tensor","",ColorType::Jet,StoreType::Cell}});
	Set_Visibility(opengl_grid_tensor_field,'2',false);
	if(opengl_grid_tensor_field!=nullptr)Bind_Callback_Key('2',&opengl_grid_tensor_field->Toggle_Draw_Func,"draw grid tensor field");}
}

void OpenGLViewerTestCase::Initialize_Volume_Data()
{
	////Volumetric data
	auto opengl_vol=Add_Object<OpenGLVolume<real> >();
	{Vector3i counts=Vector3i::Ones()*64;counts[1]*=2;
	Field<float,3> field;field.Resize(counts);
	real dx=(real)1/(real)counts[0];Grid<3> grid;grid.Initialize(counts,dx);
	Vector3 center=dx*counts.cast<real>()*(real).5;real R=(real).25;
	iterate_cell_d(iter,grid,3){const Vector3i& cell=iter.Coord();
		Vector3 pos=grid.Center(cell);real r=(pos-center).norm();
		real phi=r-R;phi=AuxFunc::Clamp(phi,(real)-1,(real)1);phi=(real).5*phi+(real).5;field(cell)=(float)phi;}	////levelset

	OpenGLColorRamp color_ramp;
	color_ramp.Add_Color((real)0,OpenGLColor(0.f,0.f,0.f,0.f));
	color_ramp.Add_Color((real).49,OpenGLColor(0.f,0.f,0.f,0.f));
	color_ramp.Add_Color((real)0.5,OpenGLColor(.1f,.9f,.1f,.02f));
	color_ramp.Add_Color((real).51,OpenGLColor(0.f,0.f,0.f,0.f));
	color_ramp.Add_Color((real)1.,OpenGLColor(0.f,0.f,0.f,0.f));

	Array<ushort> transfer;int size=256;transfer.resize(size*4);AuxFunc::Fill(transfer,(ushort)0);
	for(auto i=0;i<size;i++){real c=(real)i/(real)(size-1);OpenGLColor color=color_ramp.Lookup(c);
		for(int j=0;j<4;j++){transfer[i*4+j]=(ushort)((real)0xffff*color.rgba[j]);}}
	opengl_vol->Set_Data_Pointers(&grid,&field,&transfer);
	opengl_vol->Set_Data_Refreshed();}
}

void OpenGLViewerTestCase::Initialize_2D_Data()
{
	////Grid and fields
	{auto opengl_grid=Add_Object<OpenGLGrid>("grid2");}

	{auto opengl_grid_scalar_field=Add_Grid_Scalar_Field("grid2",{{"node_field_scalar2","",ColorType::Jet,StoreType::Node}});
		Set_Visibility(opengl_grid_scalar_field,'1',false);
		if(opengl_grid_scalar_field!=nullptr)Bind_Callback_Key('1',&opengl_grid_scalar_field->Toggle_Draw_Func,"draw grid scalar field");}

	{auto opengl_grid_scalar_field2=Add_Grid_Scalar_Field("grid2",{{"cell_field_scalar2","",ColorType::Jet,StoreType::Cell}});
		Set_Visibility(opengl_grid_scalar_field2,'2',false);
		if(opengl_grid_scalar_field2!=nullptr)Bind_Callback_Key('2',&opengl_grid_scalar_field2->Toggle_Draw_Func,"draw grid tensor field");}

	{auto opengl_grid_height_field2=Add_Grid_Height_Field("grid2",{{"node_field_scalar2","",ColorType::Jet,StoreType::Node}});
		Set_Visibility(opengl_grid_height_field2,'3',true);
		if(opengl_grid_height_field2!=nullptr)Bind_Callback_Key('2',&opengl_grid_height_field2->Toggle_Draw_Func,"draw grid height field");}
}

void OpenGLViewerTestCase::Generate_Test_Data()
{
	////Create data folder: output_dir/0/
	int frame=0;
	if(frame==0){if(!File::Directory_Exists(output_dir.c_str()))File::Create_Directory(output_dir);}
	std::string frame_dir=output_dir+"/"+std::to_string(frame);
	if(!File::Directory_Exists(frame_dir.c_str()))File::Create_Directory(frame_dir);

	////Generate vertices
	Array<Vector3> vertices;
	Vector3 p=Vector3::Unit(0)*3;
	vertices.push_back(p+Vector3::Zero());
	vertices.push_back(p+Vector3::Unit(0));
	vertices.push_back(p+Vector3::Unit(1));
	vertices.push_back(p+Vector3::Unit(2));
	Vector3 p2=Vector3::Unit(0)*5;
	vertices.push_back(p2+Vector3::Zero());
	vertices.push_back(p2+Vector3::Unit(0));
	vertices.push_back(p2+Vector3::Unit(1));
	vertices.push_back(p2+Vector3::Unit(0)+Vector3::Unit(1));
	vertices.push_back(p2+Vector3::Unit(0)*2);
	vertices.push_back(p2+Vector3::Unit(0)*2+Vector3::Unit(1)*2);

	////Generate points, p, meshes
	Points<3> points;
	for(const auto& vtx:vertices){int i=points.Add_Element();points.X(i)=vtx;}
	{std::string file_name=frame_dir+"/points";
	File::Write_Binary_To_File(file_name,points);}

	Particles<3> particles;
	for(const auto& vtx:vertices){int i=particles.Add_Element();particles.X(i)=vtx;}
	{std::string file_name=frame_dir+"/particles";
	File::Write_Binary_To_File(file_name,particles);}

	TriangleMesh<3> tri_mesh;
	(*tri_mesh.vertices)=vertices;
	tri_mesh.elements.push_back(Vector3i(4,5,6));
	tri_mesh.elements.push_back(Vector3i(5,6,7));
	{std::string file_name=frame_dir+"/triangle_mesh";
	File::Write_Binary_To_File(file_name,tri_mesh);
	file_name=frame_dir+"/triangle_mesh.txt";
	File::Write_Text_To_File(file_name,tri_mesh);}

	//const std::string mesh_data_path=Path::Data()+"/meshes/";
	//Array<std::shared_ptr<TriangleMesh<3> > > meshes;
	//Obj::Read_From_Obj_File(mesh_data_path+"/bunny.obj",meshes);
	//TriangleMesh<3> bunny=*meshes[0];
	//MeshFunc::Rescale<3>(*bunny.vertices,(real)1);
	//MeshFunc::Translate_Center_To<3>(*bunny.vertices,Vector3::Zero());
	//File::Write_Text_To_File(frame_dir+"/bunny.txt",bunny);
		
	TetrahedronMesh<3> tet_mesh;
	(*tet_mesh.vertices)=vertices;tet_mesh.elements.push_back(Vector4i(0,1,2,3));
	{std::string file_name=frame_dir+"/tetrahedron_mesh";
	File::Write_Binary_To_File(file_name,tet_mesh);}

	QuadMesh<3> quad_mesh;
	(*quad_mesh.vertices)=vertices;
	quad_mesh.elements.push_back(Vector4i(4,5,7,6));
	quad_mesh.elements.push_back(Vector4i(5,8,9,7));
	{std::string file_name=frame_dir+"/quad_mesh";
	File::Write_Binary_To_File(file_name,quad_mesh);}

	////Generate mesh fields
	{Field<real,1> field;field.Resize(Vec1i((int)(*tri_mesh.vertices).size()));
	for(int i=0;i<field.counts[0];i++)field.array[i]=(real)i;
	std::string file_name=frame_dir+"/tri_field_scalar";
	File::Write_Binary_To_File(file_name,field);}

	{Field<real,1> field;field.Resize(Vec1i((int)(*tet_mesh.vertices).size()));
	for(int i=0;i<field.counts[0];i++)field.array[i]=(real)i;
	std::string file_name=frame_dir+"/tet_field_scalar";
	File::Write_Binary_To_File(file_name,field);}

	{Field<Vector3,1> field;field.Resize(Vec1i((int)(*tet_mesh.vertices).size()));
	for(int i=0;i<field.counts[0];i++)field.array[i]=Vector3::Ones();
	std::string file_name=frame_dir+"/tet_field_vector";
	File::Write_Binary_To_File(file_name,field);}

	{Field<Matrix3,1> field;field.Resize(Vec1i((int)(*tet_mesh.vertices).size()));
	Matrix3 m;m<<1,0,0,0,2,0,0,0,3;
	for(int i=0;i<field.counts[0];i++)field.array[i]=m;
	std::string file_name=frame_dir+"/tet_field_tensor";
	File::Write_Binary_To_File(file_name,field);}

	{Field<Matrix3,1> field;field.Resize(Vec1i((int)tet_mesh.elements.size()));
	Matrix3 m;m<<1,0,0,0,2,0,0,0,3;
	for(int i=0;i<field.counts[0];i++)field.array[i]=m;
	std::string file_name=frame_dir+"/tet_field_cell_tensor";
	File::Write_Binary_To_File(file_name,field);}

    ////Generate grid
	Grid<3> grid(Vector3i::Ones()*4,(real)1);
	{std::string file_name=frame_dir+"/grid";
	File::Write_Binary_To_File(file_name,grid);}

	////Generate grid fields
	{Field<Vector3,3> field;field.Resize(grid.node_counts);field.Fill(Vector3::Zero());
	Vector3i node=Vector3i::Ones();field(node)=Vector3::Ones()*(real)1;
	std::string file_name=frame_dir+"/displacement";
	File::Write_Binary_To_File(file_name,field);}

	{Field<real,3> field;field.Resize(grid.node_counts);
	iterate_node_d(iter,grid,3){const Vector3i& node=iter.Coord();
		field(node)=grid.Node(node).norm();}
	std::string file_name=frame_dir+"/node_field_scalar";
	File::Write_Binary_To_File(file_name,field);}

	{Field<Vector3,3> field;field.Resize(grid.node_counts);
	iterate_node_d(iter,grid,3){const Vector3i& node=iter.Coord();field(node)=node.cast<real>();}
	std::string file_name=frame_dir+"/node_field_vector";
	File::Write_Binary_To_File(file_name,field);}

	{Field<Matrix3,3> field;field.Resize(grid.node_counts);
	Matrix3 m;m<<1,0,0,0,2,0,0,0,3;
	iterate_node_d(iter,grid,3){const Vector3i& node=iter.Coord();field(node)=(real).1*m;}
	std::string file_name=frame_dir+"/node_field_tensor";
	File::Write_Binary_To_File(file_name,field);}

	{Field<real,3> field;field.Resize(grid.cell_counts);
	iterate_cell_d(iter,grid,3){const Vector3i& cell=iter.Coord();field(cell)=(real)1;}
	std::string file_name=frame_dir+"/cell_field_scalar";
    File::Write_Binary_To_File(file_name,field);}

	{Field<Vector3,3> field;field.Resize(grid.cell_counts);
	iterate_cell_d(iter,grid,3){const Vector3i& cell=iter.Coord();field(cell)=Vector3::Ones();}
	std::string file_name=frame_dir+"/cell_field_vector";
    File::Write_Binary_To_File(file_name,field);}

	{Field<Matrix3,3> field;field.Resize(grid.cell_counts);
	Matrix3 m;m<<1,0,0,0,2,0,0,0,3;
	iterate_cell_d(iter,grid,3){const Vector3i& cell=iter.Coord();field(cell)=m;}
	std::string file_name=frame_dir+"/cell_field_tensor";
    File::Write_Binary_To_File(file_name,field);}

	////2d data
	int n=4;Grid<2> grid2(Vector2i::Ones()*n,(real)1/(real)n);
	{std::string file_name=frame_dir+"/grid2";
	Grid<3> grid3;Dim_Conversion(grid2,grid3);
	File::Write_Binary_To_File(file_name,grid3);}

	{Field<real,2> field;field.Resize(grid2.node_counts);
	iterate_node_d(iter,grid2,2){const Vector2i& node=iter.Coord();const Vector2& pos=grid2.Node(node);
		field(node)=pos.norm();}
	std::string file_name=frame_dir+"/node_field_scalar2";
	Field<real,3> field3;Dim_Conversion(field,field3);
    File::Write_Binary_To_File(file_name,field3);}

	{Field<real,2> field;field.Resize(grid2.cell_counts);
	iterate_cell_d(iter,grid2,2){const Vector2i& cell=iter.Coord();const Vector2& pos=grid2.Center(cell);
		field(cell)=pos.norm();}
	std::string file_name=frame_dir+"/cell_field_scalar2";
	Field<real,3> field3;Dim_Conversion(field,field3);
    File::Write_Binary_To_File(file_name,field3);}
}
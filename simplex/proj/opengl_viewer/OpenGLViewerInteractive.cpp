//////////////////////////////////////////////////////////////////////////
// OpenGL viewer interactive
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <File.h>
#include "OpenGLViewerInteractive.h"
#include "OpenGLMesh.h"
#include "OpenGLUbos.h"
#include "OpenGLWindow.h"

void OpenGLViewerInteractive::Initialize_Data()
{
	opengl_tri_mesh=Add_Interactive_Object<OpenGLTriangleMesh>();
	auto tri_mesh=&opengl_tri_mesh->mesh;
	std::string file_name=Path::Data()+"/meshes/bunny.txt";
	File::Read_Text_From_File(file_name,*tri_mesh);
	//MeshFunc::Initialize_Sphere_Mesh((real).8,tri_mesh);
	//MeshFunc::Initialize_Cone_Mesh((real).4,(real).8,8,tri_mesh);
	//MeshFunc::Initialize_Cylinder_Mesh((real).4,(real).8,8,tri_mesh);
	MeshFunc::Update_Normals(*tri_mesh);
	opengl_tri_mesh->Set_Data_Refreshed();
	opengl_tri_mesh->Initialize();

	//Set_Polygon_Mode(opengl_tri_mesh,PolygonMode::Wireframe);
	//Set_Polygon_Mode(opengl_tri_mesh,PolygonMode::Fill);
	//Set_Shading_Mode(opengl_tri_mesh,ShadingMode::Lighting);

	Box<3> bounding_box=MeshFunc::Bounding_Box<3>(*tri_mesh->vertices);
	real dx=bounding_box.Edge_Lengths().maxCoeff()/32;
	spatial_hashing.reset(new SpatialHashing<3>());
	Grid<3> grid(Vector3i::Ones()*32,dx,bounding_box.min_corner);
	spatial_hashing->Initialize(grid);
	for(int i=0;i<tri_mesh->elements.size();i++){
		Vector3 center=MeshFunc::Element_Center(*tri_mesh->vertices,tri_mesh->elements,i);
		centers.push_back(center);}
	spatial_hashing->Update_Voxels(centers);

	opengl_point=Add_Interactive_Object<OpenGLPoint>();
	opengl_point->pos=Vector3::Unit(0);
	opengl_point->point_size=(real)4;
	opengl_point->Set_Data_Refreshed();
	opengl_point->Initialize();

	opengl_sphere=Add_Interactive_Object<OpenGLSphere>();
	opengl_sphere->pos=Vector3::Unit(0)*(real)2;
	opengl_sphere->radius=(real).02;
	opengl_sphere->Set_Data_Refreshed();
	opengl_sphere->Initialize();

	opengl_triangle=Add_Interactive_Object<OpenGLTriangle>();
	opengl_triangle->vtx={Vector3::Unit(0),Vector3::Unit(1),Vector3::Unit(2)};
	opengl_triangle->Set_Data_Refreshed();
	Set_Color(opengl_triangle,OpenGLColor::Green());
	Set_Line_Width(opengl_triangle,4.f);
	opengl_triangle->Initialize();

	opengl_arrow=Add_Interactive_Object<OpenGLArrow>();
	opengl_arrow->start=Vector3::Ones()*(real)1;
	opengl_arrow->end=Vector3::Ones()*(real)2;
	Set_Line_Width(opengl_arrow,4.f);
	opengl_arrow->Set_Data_Refreshed();
	opengl_arrow->Initialize();

	opengl_circle=Add_Interactive_Object<OpenGLCircle>();
	opengl_circle->pos=Vector3::Unit(0)*(real)2;
	opengl_circle->radius=(real).5;
	Set_Line_Width(opengl_circle,4.f);
	opengl_circle->Set_Data_Refreshed();
	opengl_circle->Initialize();

	auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
	OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
	OpenGLUbos::Update_Lights_Ubo();
	
	Reset_Interactive_Objects(0);
}

void OpenGLViewerInteractive::Move_Mesh(const Vector3& translate)
{
	//for(auto& v:(*tet_mesh->vertices))v+=translate;
	//opengl_tet_mesh->Update_Data_To_Render();
}

void OpenGLViewerInteractive::Move_Left()
{Move_Mesh(Vector3::Unit(0)*(real)-.1);}

void OpenGLViewerInteractive::Move_Right()
{Move_Mesh(Vector3::Unit(0)*(real)+.1);}

void OpenGLViewerInteractive::Initialize_Common_Callback_Keys()
{
	Base::Initialize_Common_Callback_Keys();
	Bind_Callback_Key('=',&Move_Right_Func,"move right");
	Bind_Callback_Key('-',&Move_Left_Func,"move left");
}

void OpenGLViewerInteractive::Reset_Interactive_Objects(const int idx)
{
	Vector3 center=centers[idx];

	ArrayF<Vector3,3> tri_vtx;
	for(int i=0;i<3;i++)tri_vtx[i]=(*opengl_tri_mesh->mesh.vertices)[opengl_tri_mesh->mesh.elements[idx][i]];
	opengl_triangle->vtx=tri_vtx;
	opengl_triangle->Set_Data_Refreshed();
	opengl_triangle->Update_Data_To_Render();

	Vector3 normal=(tri_vtx[1]-tri_vtx[0]).cross(tri_vtx[2]-tri_vtx[0]);normal.normalize();
	opengl_arrow->start=center+normal*(real).8;
	opengl_arrow->end=center+normal*(real).3;
	opengl_arrow->Set_Data_Refreshed();
	opengl_arrow->Update_Data_To_Render();

	opengl_point->pos=center+normal*(real).1;
	opengl_point->Set_Data_Refreshed();
	opengl_point->Update_Data_To_Render();

	opengl_sphere->pos=center+normal*(real).2;
	opengl_sphere->Set_Data_Refreshed();
	opengl_sphere->Update_Data_To_Render();
}

bool OpenGLViewerInteractive::Mouse_Click(int left,int right,int mid,int x,int y,int w,int h)
{
	if(left!=1){is_focus_captured=false;return false;}
	GLfloat depth=opengl_window->Win_Depth(x,y);
	if(depth<1.f){
		is_focus_captured=true;
		Vector3f pos=opengl_window->Win_Coord_To_World_Coord(x,y);
		int idx=spatial_hashing->Find_Nearest_Nb(pos.cast<real>(),centers);
		if(idx==-1){is_focus_captured=false;return false;}
		Reset_Interactive_Objects(idx);
		return true;}
	
	is_focus_captured=false;
	return false;
}

bool OpenGLViewerInteractive::Mouse_Drag(int x,int y,int w,int h)
{
	return is_focus_captured;
}

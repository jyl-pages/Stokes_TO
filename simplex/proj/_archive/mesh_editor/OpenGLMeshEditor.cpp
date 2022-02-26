//#####################################################################
// OpenGL Mesh Editor
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#include <File.h>
#include "OpenGLMeshEditor.h"
#include "OpenGLMesh.h"
#include "OpenGLUbos.h"
#include "OpenGLWindow.h"
#ifdef USE_TINY_OBJ_LOADER
#include "TinyObjLoader.h"
#endif

void OpenGLMeshEditor::Initialize_Data()
{
	opengl_tri_mesh=Add_Interactive_Object<OpenGLTriangleMesh>();
	auto tri_mesh=&opengl_tri_mesh->mesh;
	bool read=File::Read_Text_From_File(file_name,*tri_mesh);
	Array<std::shared_ptr<TriangleMesh<3> > > tri_meshes;
	Obj::Read_From_Obj_File<TriangleMesh<3> >(file_name,tri_meshes);
	bool read=tri_meshes.size()>0;
	*tri_mesh=*tri_meshes[0];

	if(!read){std::cerr<<"Error: [OpenGLMeshEditor] cannot read mesh file "<<file_name<<std::endl;return;}
	else{std::cout<<"Read mesh: "<<file_name<<std::endl;}

	////fix the 1-index problem
	//for(auto& e:(*tri_mesh).elements){
	//	for(int i=0;i<3;i++)e[i]-=1;}
	//File::Write_Text_To_File(file_name,*tri_mesh);

	MeshFunc::Update_Normals(*tri_mesh);
	opengl_tri_mesh->Set_Data_Refreshed();
	opengl_tri_mesh->Initialize();

	Box<3> bounding_box=MeshFunc::Bounding_Box<3>(*tri_mesh->vertices);
	real dx=bounding_box.Edge_Lengths().maxCoeff()/128;
	spatial_hashing.reset(new SpatialHashing64());
	Grid<3> grid(Vector3i::Ones()*32,dx,bounding_box.min_corner);
	spatial_hashing->Initialize(grid);
	for(int i=0;i<tri_mesh->elements.size();i++){
		Vector3 center=MeshFunc::Element_Center(*tri_mesh->vertices,tri_mesh->elements,i);
		centers.push_back(center);}
	spatial_hashing->Update_Voxels(centers);

	opengl_arrow=Add_Interactive_Object<OpenGLArrow>();
	opengl_arrow->start=Vector3::Ones()*(real)1;
	opengl_arrow->end=Vector3::Ones()*(real)2;
	Set_Line_Width(opengl_arrow,4.f);
	opengl_arrow->Set_Data_Refreshed();
	opengl_arrow->Initialize();

	auto dir_light=OpenGLUbos::Add_Directional_Light(glm::vec3(-1.f,-.1f,-.2f));
	OpenGLUbos::Set_Ambient(glm::vec4(.1f,.1f,.1f,1.f));
	OpenGLUbos::Update_Lights_Ubo();
	
	Reset_Interactive_Objects(0);
}

void OpenGLMeshEditor::Move_Mesh(const Vector3& translate)
{
	//for(auto& v:(*tet_mesh->vertices))v+=translate;
	//opengl_tet_mesh->Update_Data_To_Render();
}

void OpenGLMeshEditor::Move_Left()
{Move_Mesh(Vector3::Unit(0)*(real)-.1);}

void OpenGLMeshEditor::Move_Right()
{Move_Mesh(Vector3::Unit(0)*(real)+.1);}

void OpenGLMeshEditor::Initialize_Common_Callback_Keys()
{
	Base::Initialize_Common_Callback_Keys();
	Bind_Callback_Key('=',&Move_Right_Func,"move right");
	Bind_Callback_Key('-',&Move_Left_Func,"move left");
}

void OpenGLMeshEditor::Reset_Interactive_Objects(const int idx)
{
	if(idx>(centers.size()-1)){return;}

	Vector3 center=centers[idx];
	ArrayF<Vector3,3> tri_vtx;
	for(int i=0;i<3;i++)tri_vtx[i]=(*opengl_tri_mesh->mesh.vertices)[opengl_tri_mesh->mesh.elements[idx][i]];

	Vector3 normal=(tri_vtx[1]-tri_vtx[0]).cross(tri_vtx[2]-tri_vtx[0]);normal.normalize();
	opengl_arrow->start=center+normal*(real).8;
	opengl_arrow->end=center+normal*(real).3;
	opengl_arrow->Set_Data_Refreshed();
	opengl_arrow->Update_Data_To_Render();
}

bool OpenGLMeshEditor::Mouse_Click(int left,int right,int mid,int x,int y,int w,int h)
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

bool OpenGLMeshEditor::Mouse_Drag(int x,int y,int w,int h)
{
	return is_focus_captured;
}

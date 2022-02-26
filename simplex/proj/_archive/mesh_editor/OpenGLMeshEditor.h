//#####################################################################
// OpenGL Mesh Editor
// Author: 
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __OpenGLMeshEditor_h__
#define __OpenGLMeshEditor_h__
#include "Mesh.h"
#include "SpatialHashing.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"
#include "AuxFunc.h"

class OpenGLMeshEditor : public OpenGLViewer
{typedef OpenGLViewer Base;using SpatialHashing64=SpatialHashing<3,64>;
public:using Base::opengl_window;
	OpenGLTriangleMesh* opengl_tri_mesh=nullptr;
	OpenGLArrow* opengl_arrow=nullptr;

	std::shared_ptr<SpatialHashing64> spatial_hashing;
	Array<Vector3> centers;
	bool is_focus_captured=false;

	////mesh properties
	std::string file_name=Path::Data()+"/meshes/fish.txt";

	virtual void Initialize_Data();
	void Move_Mesh(const Vector3& translate);

	virtual void Move_Left();
	Define_Function_Object(OpenGLMeshEditor,Move_Left);

	virtual void Move_Right();
	Define_Function_Object(OpenGLMeshEditor,Move_Right);

	virtual void Initialize_Common_Callback_Keys();
	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h);
	virtual bool Mouse_Drag(int x,int y,int w,int h);

	void Reset_Interactive_Objects(const int idx);
};
#endif
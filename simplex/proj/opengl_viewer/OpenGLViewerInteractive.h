//////////////////////////////////////////////////////////////////////////
// OpenGL viewer interactive
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLViewerInteractive_h__
#define __OpenGLViewerInteractive_h__
#include "Mesh.h"
#include "SpatialHashing.h"
#include "OpenGLViewer.h"
#include "OpenGLMarkerObjects.h"

class OpenGLViewerInteractive : public OpenGLViewer
{typedef OpenGLViewer Base;
public:using Base::opengl_window;
	OpenGLTriangleMesh* opengl_tri_mesh=nullptr;
	OpenGLPoint* opengl_point=nullptr;
	OpenGLSphere* opengl_sphere=nullptr;
	OpenGLTriangle* opengl_triangle=nullptr;
	OpenGLArrow* opengl_arrow=nullptr;
	OpenGLCircle* opengl_circle=nullptr;

	std::shared_ptr<SpatialHashing<3> > spatial_hashing;
	Array<Vector3> centers;
	bool is_focus_captured=false;

	virtual void Initialize_Data();
	void Move_Mesh(const Vector3& translate);

	virtual void Move_Left();
	Define_Function_Object(OpenGLViewerInteractive,Move_Left);

	virtual void Move_Right();
	Define_Function_Object(OpenGLViewerInteractive,Move_Right);

	virtual void Initialize_Common_Callback_Keys();
	virtual bool Mouse_Click(int left,int right,int mid,int x,int y,int w,int h);
	virtual bool Mouse_Drag(int x,int y,int w,int h);

	void Reset_Interactive_Objects(const int idx);
};
#endif
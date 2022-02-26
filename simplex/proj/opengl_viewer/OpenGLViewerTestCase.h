//////////////////////////////////////////////////////////////////////////
// Opengl viewer test
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLViewerTestCase_h__
#define __OpenGLViewerTestCase_h__
#include "OpenGLViewer.h"

//////////////////////////////////////////////////////////////////////////
////OpenGLViewerTestCase contains test cases for all rendering primitives
//////////////////////////////////////////////////////////////////////////

class OpenGLViewerTestCase : public OpenGLViewer
{typedef OpenGLViewer Base;
public:
	OpenGLViewerTestCase():Base(){}

	virtual void Initialize_Data();
	void Initialize_Mesh_Data();
	void Initialize_Grid_Data();
	void Initialize_Volume_Data();
	void Initialize_2D_Data();
	void Generate_Test_Data();
};
#endif
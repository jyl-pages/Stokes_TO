//////////////////////////////////////////////////////////////////////////
// Opengl mesh vectors
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLMeshVectors_h__
#define __OpenGLMeshVectors_h__
#include "OpenGLMeshField.h"
template<class T,class T_MESH>
class OpenGLMeshVectors : public OpenGLMeshField<Vector<T,3>,T_MESH,1>
{typedef OpenGLMeshField<Vector<T,3>,T_MESH,1> Base;
public:
	using Base::name;using Base::field_draw_type;using Base::opengl_mesh;using Base::field_val;
	using Base::Update_Data_To_Render_Mesh_Lines;using Base::Display_Lines;
	OpenGLMeshVectors(const OpenGLMesh<T_MESH>* _opengl_mesh=nullptr):Base(_opengl_mesh){name="mesh_vectors";field_draw_type=FieldDrawType::Line;}
	virtual void Update_Data_To_Render(){Update_Data_To_Render_Mesh_Lines(opengl_mesh,field_val);}
	virtual void Display() const {Display_Lines();}
};
#endif
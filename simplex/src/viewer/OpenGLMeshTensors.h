//////////////////////////////////////////////////////////////////////////
// Opengl mesh tensors
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLMeshTensors_h__
#define __OpenGLMeshTensors_h__
#include "OpenGLMeshField.h"
template<class T,class T_MESH>
class OpenGLMeshTensors : public OpenGLMeshField<Matrix<T,3>,T_MESH,1>
{typedef OpenGLMeshField<Matrix<T,3>,T_MESH,1> Base;
public:
	using Base::field_draw_type;using Base::name;using Base::opengl_mesh;using Base::field_val;
	using Base::Display_Lines;using Base::Update_Data_To_Render_Mesh_Lines;
	OpenGLMeshTensors(const OpenGLMesh<T_MESH>* _opengl_mesh=nullptr):Base(_opengl_mesh){name="mesh_tensors";field_draw_type=FieldDrawType::Line;}
	virtual void Update_Data_To_Render(){Update_Data_To_Render_Mesh_Lines(opengl_mesh,field_val);}
	virtual void Display() const {Display_Lines();}
};
#endif
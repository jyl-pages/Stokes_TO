//////////////////////////////////////////////////////////////////////////
// Opengl vectors
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLVectors_h__
#define __OpenGLVectors_h__
#include "OpenGLObject.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLPrimitives.h"
#include "OpenGLUbos.h"

class OpenGLVectors : public OpenGLObject
{typedef OpenGLObject Base;
public:
    const Array<Vector3>* bases=nullptr;
    const Array<Vector3>* vectors=nullptr;
	bool draw_arrow=true;

    OpenGLVectors(const Array<Vector3>* _bases=nullptr,const Array<Vector3>* _vectors=nullptr)
		:bases(_bases),vectors(_vectors){name="vectors";color=OpenGLColor::Green();polygon_mode=PolygonMode::Fill;}

	void Set_Data_Pointers(const Array<Vector3>* _bases,const Array<Vector3>* _vectors){bases=_bases;vectors=_vectors;}

	virtual void Initialize()
	{
		Base::Initialize();
		Add_Shader_Program(OpenGLShaderLibrary::Get_Shader("psize_ucolor"));
	}

	virtual void Update_Data_To_Render()
	{
		if(!initialized)Initialize();

		Clear_OpenGL_Arrays();
		for(size_type i=0;i<(*bases).size();i++){
			const Vector3& p0=(*bases)[i];
			const Vector3 p1=p0+(*vectors)[i]*scale;
			OpenGL_Vertex(p0,opengl_vertices);	////position, 3 floats
			OpenGL_Vertex(p1,opengl_vertices);}	////position, 3 floats
		Set_OpenGL_Vertices();
		Set_OpenGL_Vertex_Attribute(0,3,0,0);	////position
		Clear_OpenGL_Arrays();
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
		glDrawArrays(GL_LINES,0,(GLsizei)(vtx_size/3));
		shader->End();}
    }

protected:
	////TODO: draw arrow
};
#endif
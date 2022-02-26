//////////////////////////////////////////////////////////////////////////
// Opengl tracker vectors
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OrientedParticle_h__
#define __OrientedParticle_h__
#include "Common.h"
#include "AuxFunc.h"
#include "OpenGLObject.h"
#include "OpenGLShaderLibrary.h"
#include "OpenGLUbos.h"
#include "OpenGLTexture.h"
#include "OpenGLPrimitives.h"
#include "ArrayIO.h"

class LagrangianQuantity {
public:
	std::string file_name;
	bool visible = true;
	OpenGLColor color;
	virtual void Initialize(std::string _filename, bool _visible = true, OpenGLColor _color = OpenGLColor::Yellow());
	virtual void Read_Data(const std::string& frame_path) = 0;
};

class LagrangianScalar :public LagrangianQuantity {
public:
	Array<real> data;
	virtual void Read_Data(const std::string& frame_path);
};

class LagrangianVector :public LagrangianQuantity {
public:
	Array<Vector3> data;
	virtual void Read_Data(const std::string& frame_path);
};

class OrientedParticle :public OpenGLObject
{
	typedef OpenGLObject Base;
public:
	LagrangianVector pos;
	LagrangianVector normal;
	Array<LagrangianScalar> scalar_qnts;
	Array<LagrangianVector> vector_qnts;

	GLfloat point_size = 6.f;

	virtual void Initialize(void);
	void Initialize_Data(std::string output_dir, const int frame);
	virtual void Refresh(const int frame);
	virtual void Update_Data_To_Render(void);
	void Update_Model_Matrix_Helper(const float* pos, const float* dir, const float radius, glm::mat4& model_matrix) const;
	void Display_Tracker_Points(void);
	virtual void Display()const;
	void Add_Scalar(std::string _name, bool _visible = true, OpenGLColor _color = OpenGLColor::Red());
	void Add_Vector(std::string _name, bool _visible = true, OpenGLColor _color = OpenGLColor::Blue());
};

#endif
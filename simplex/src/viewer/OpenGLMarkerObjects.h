//////////////////////////////////////////////////////////////////////////
// Opengl interactive object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLMarkerObjects_h__
#define __OpenGLMarkerObjects_h__
#include "glm.hpp"
#include "Mesh.h"
#include "OpenGLObject.h"

class OpenGLAxes : public OpenGLObject
{public:typedef OpenGLObject Base;
	real axis_length=(real)1;
	bool use_2d_display=false;

	OpenGLAxes(){name="axes";polygon_mode=PolygonMode::Wireframe;}

	virtual void Initialize();
	virtual void Display() const;
};

class OpenGLSquare: public OpenGLObject
{public:typedef OpenGLObject Base;
	real axis_length=(real)1;
	Vector3 pos=Vector3::Zero();
	glm::mat4 model=glm::mat4(1.f);

	OpenGLSquare(){name="box";color=OpenGLColor::Red();polygon_mode=PolygonMode::Wireframe;}

	virtual void Initialize();
	virtual void Display() const;
	virtual void Update_Model_Matrix();
};

class OpenGLSegments: public OpenGLObject
{public:typedef OpenGLObject Base;
	real axis_length=(real)1;
	Array<Vector3> vertices;
	glm::mat4 model=glm::mat4(1.f);

	OpenGLSegments(){name="segments";color=OpenGLColor::Black();}

	virtual void Initialize();
	virtual void Display() const;
};

class OpenGLPoint : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;
	Vector3 pos=Vector3::Zero();
	GLfloat point_size=16.f;

	OpenGLPoint(){name="point";color=OpenGLColor::Red();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
};

class OpenGLTriangle : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;using Base::line_width;
	ArrayF<Vector3,3> vtx;

	OpenGLTriangle(){name="triangle";color=OpenGLColor::Red();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
};

class OpenGLMarkerTriangleMesh : public OpenGLObject
{public:typedef OpenGLObject Base;using Base::color;
	TriangleMesh<3> mesh;
	glm::mat4 model=glm::mat4(1.f);

	OpenGLMarkerTriangleMesh(){name="interactive_triangle_mesh";color=OpenGLColor::Green();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
	virtual void Update_Model_Matrix(){}

protected:
	virtual void Update_Mesh_Data_To_Render();
};

class OpenGLCircle : public OpenGLObject
{public: typedef OpenGLObject Base;
	Vector3 pos=Vector3::Zero();
	real radius=(real).1;
	glm::mat4 model=glm::mat4(1.f);

	int n=32;
	Array<Vector3> vtx;

	OpenGLCircle(){name="circle";color=OpenGLColor::Green();polygon_mode=PolygonMode::Fill;}

	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Display() const;
	virtual void Update_Model_Matrix();
	virtual void Display_Multiple_Instances(const Array<Vector3>& centers) const;
	virtual void Display_Multiple_Instances(const Array<Vector3>& centers,const Array<real>& radii) const;

protected:
	void Update_Model_Matrix_Helper(const Vector3& pos,const real r,glm::mat4& model_matrix) const;
};

class OpenGLSphere : public OpenGLMarkerTriangleMesh
{public:typedef OpenGLMarkerTriangleMesh Base;
	using Base::color;using Base::mesh;using Base::model;
	Vector3 pos=Vector3::Zero();
	real radius=(real).1;

	OpenGLSphere(){name="sphere";color=OpenGLColor::Red();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Model_Matrix();
	virtual void Display_Multiple_Instances(const Array<Vector3>& centers) const;

protected:
	void Update_Model_Matrix_Helper(const Vector3& pos,glm::mat4& model_matrix) const;
};

class OpenGLArrow : public OpenGLMarkerTriangleMesh
{public:typedef OpenGLMarkerTriangleMesh Base;
	using Base::color;using Base::line_width;using Base::mesh;using Base::model;
	Vector3 start;
	Vector3 end;

	OpenGLArrow(){name="arrow";color=OpenGLColor::Blue();polygon_mode=PolygonMode::Fill;}
	
	virtual void Initialize();
	virtual void Update_Model_Matrix();
};
#endif
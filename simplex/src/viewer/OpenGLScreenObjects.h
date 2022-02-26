//////////////////////////////////////////////////////////////////////////
// Opengl screen object
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLScreenObjects_h__
#define __OpenGLScreenObjects_h__
#include "glm.hpp"
#include "OpenGLObject.h"

class OpenGLWindow;

class OpenGLScreenObject : public OpenGLObject
{public:typedef OpenGLObject Base;
	Box<2> box=Box<2>(Vector2::Zero(),Vector2::Ones());
	Vector2 pos=Vector2::Zero();	////starting from the bottom-left corner
	int res=8;
	glm::mat4 model;

	OpenGLScreenObject(const Vector2& _pos,const Vector2& size):
		box(Vector2::Zero(),size),pos(_pos){name="screen_object";Set_Pos(pos);}

	virtual void Initialize();
	void Set_Pos(const Vector2& _pos);
};

class OpenGLText : public OpenGLScreenObject
{public:typedef OpenGLScreenObject Base;
	Array<std::string> texts;
	Array<Vector2> offsets;
	const glm::mat4* ortho;
	
	OpenGLText(const Vector2& _pos=Vector2::Zero()):Base(_pos,Vector2(128,128)),ortho(nullptr){name="text";color=OpenGLColor::Red(.5f);}

	void Set_Data_Pointers(const glm::mat4* _ortho=nullptr){ortho=_ortho;}
	void Clear_Data(){texts.clear();offsets.clear();}
	virtual void Display() const;
};

class OpenGLText3D : public OpenGLText
{public:using Base=OpenGLText;
	OpenGLWindow* opengl_window=nullptr;
	Vector3 pos_3d=Vector3::Zero();

	OpenGLText3D(OpenGLWindow* _opengl_window=nullptr):Base(),opengl_window(_opengl_window){name="text_3d";}
	void Set_Pos(const Vector3& _pos_3d);
	virtual void Display() const;
protected:
	Vector2 Proj_To_Win_Pos(const Vector3& pos_3d) const;
};

class OpenGLTextArray3D : public OpenGLText3D
{public:using Base=OpenGLText3D;
	OpenGLTextArray3D(OpenGLWindow* _opengl_window=nullptr):Base(_opengl_window){name="text_array_3d";}
	Array<Vector3> pos_3d_array;

	virtual void Display() const;
	void Update_Vertex_Data(const Array<Vector3>& vertices);
	void Update_Vertex_Data(const Array<Vector3>& vertices,const Hashset<int>& selected_vertices);
	void Update_Element_Data(const Array<Vector3>& vertices,const Array<Vector3i>& elements);
};

class OpenGLBar : public OpenGLScreenObject
{public:typedef OpenGLScreenObject Base;
	real bar_width=(real)256;
	real bar_height=(real)16;

	OpenGLBar(const Vector2& _pos,const Vector2& size=Vector2(256,16))
		:Base(_pos,size),bar_width(size[0]),bar_height(size[1]){name="bar";color=OpenGLColor::Blue();}

	virtual void Initialize();
	virtual void Update_Data_To_Render();
	virtual void Set_Bar_Mesh();
	virtual void Display() const;
};

class OpenGLColorBar : public OpenGLBar
{public:typedef OpenGLBar Base;
protected:
	int bar_vtx_size;
	int scale_vtx_size;
	OpenGLColorMapper* color_mapper_ptr;
	int* data_idx_ptr;

	bool draw_scalars;
	OpenGLText scalar_texts;
	Array<Vector2> scalar_vertices;
public:
	OpenGLColorBar(const Vector2& _pos=Vector2(64,64),const Vector2& size=Vector2(16,256))
		:Base(_pos,size),bar_vtx_size(0),scale_vtx_size(0),color_mapper_ptr(nullptr),data_idx_ptr(nullptr),
		draw_scalars(true),scalar_texts(_pos)
	{name="color_bar";bar_width=size[0];bar_height=size[1];}

	void Set_Data_Pointer(OpenGLColorMapper* _color_mapper_ptr,int* _data_idx_ptr);
	virtual void Initialize();
	virtual void Set_Bar_Mesh();
	virtual void Display() const;
	void Update_Scalar_Texts(){}
};

class OpenGLBackground : public OpenGLObject
{public:typedef OpenGLObject Base;
	Box<2> box;	
	std::string tex_name;
	real depth;
	bool use_fbo_tex;
	Array<OpenGLColor> mix_colors;

	OpenGLBackground();

	void Set_Box(const Vector2& min_corner,const Vector2& max_corner){box=Box<2>(min_corner,max_corner);}
	void Set_Texture(const std::string& _tex_name){use_vtx_tex=true;tex_name=_tex_name;}
	void Set_Depth(const real _depth){depth=_depth;}
	void Set_Fbo(){}
	virtual void Initialize();
	virtual void Display() const;
};

#endif
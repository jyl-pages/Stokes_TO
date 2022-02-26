//////////////////////////////////////////////////////////////////////////
// Opengl common
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLCommon_h__
#define __OpenGLCommon_h__
#include "glm.hpp"
#include "gtc/type_ptr.hpp"
#include "Common.h"
#include "OpenGLColor.h"

////Callback function object
#define Define_Function_Object(class_name,function_name) \
std::function<void(void)> function_name##_Func=std::bind(&class_name::function_name,this);

////OpenGL Enums
enum class ColorType:int{Jet=0,Hot=1,Den=2,Mat=3,MC=4};
enum class StoreType:int{None=-1,Node=0,Cell=1};
enum class FieldDrawType:int{Color=0,Line=1};	
enum class TensorFieldDrawType:int{Eigenvector=0,Frame=1};
enum class TextureType:int{Tx1d,Tx2d,Tx3d,Tx3d2,Tx3d3,Tx3d4,TxCube};
enum class PolygonMode:int{Fill=0,Wireframe,SurfOnly,SurfWithEdge};
//Picking: encode id of object in color
enum class ShadingMode :int { None = 0, Multicolor, Lighting, TexOnly, TexLighting, Sprite, ColorSprite, SizeSprite, TexSprite, Shadow, EnvMapping, ColorEncoding };
enum class ParticleMode:int{Dot=0,Circle,Sphere};

glm::vec3 inline ToGlm(const Vector3& v){return glm::vec3((float)v[0],(float)v[1],(float)v[2]);}
glm::vec4 inline ToGlm(const Vector4& v){return glm::vec4((float)v[0],(float)v[1],(float)v[2],(float)v[3]);}
glm::mat3 inline ToGlm(const Matrix3& m){return glm::make_mat3(m.data());}
glm::mat4 inline ToGlm(const Matrix4& m){return glm::make_mat4(m.data());}

OPENGL_COLOR inline ID_To_Color(int idx) {
	int r = (idx & 0x000000FF) >> 0;
	int g = (idx & 0x0000FF00) >> 8;
	int b = (idx & 0x00FF0000) >> 16;
	return OPENGL_COLOR(r / 255.0f, g / 255.0f, b / 255.0f, 1.0f);
}

#endif
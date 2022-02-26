//////////////////////////////////////////////////////////////////////////
// Opengl primitives
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __OpenGLPrimitives_h__
#define __OpenGLPrimitives_h__
#include <GL/glew.h>
#include "Common.h"

inline void OpenGL_Vertex(const int& v,Array<GLuint>& vertices)
{vertices.push_back((GLuint)v);}

inline void OpenGL_Vertex(const Vector2i& v,Array<GLuint>& vertices)
{vertices.push_back((GLuint)v[0]);vertices.push_back((GLuint)v[1]);}

inline void OpenGL_Vertex(const Vector3i& v,Array<GLuint>& vertices)
{vertices.push_back((GLuint)v[0]);vertices.push_back((GLuint)v[1]);vertices.push_back((GLuint)v[2]);}

inline void OpenGL_Vertex(const Vector4i& v,Array<GLuint>& vertices)
{vertices.push_back((GLuint)v[0]);vertices.push_back((GLuint)v[1]);vertices.push_back((GLuint)v[2]);vertices.push_back((GLuint)v[3]);}

inline void OpenGL_Vertex(const int v0,const int v1,const int v2,Array<GLuint>& vertices)
{vertices.push_back((GLuint)v0);vertices.push_back((GLuint)v1);vertices.push_back((GLuint)v2);}

inline void OpenGL_Vertex(const GLfloat v,Array<GLfloat>& vertices)
{vertices.push_back(v);}

inline void OpenGL_Vertex(const Vector2& v,Array<GLfloat>& vertices)
{vertices.push_back((GLfloat)v[0]);vertices.push_back((GLfloat)v[1]);}

inline void OpenGL_Vertex(const Vector3& v,Array<GLfloat>& vertices)
{vertices.push_back((GLfloat)v[0]);vertices.push_back((GLfloat)v[1]);vertices.push_back((GLfloat)v[2]);}

inline void OpenGL_Vertex4(const Vector2& v,Array<GLfloat>& vertices,const GLfloat placeholder=(GLfloat)0)
{vertices.push_back((GLfloat)v[0]);vertices.push_back((GLfloat)v[1]);vertices.push_back(placeholder);vertices.push_back(placeholder);}

inline void OpenGL_Vertex4(const Vector3& v,Array<GLfloat>& vertices,const GLfloat placeholder=(GLfloat)0)
{vertices.push_back((GLfloat)v[0]);vertices.push_back((GLfloat)v[1]);vertices.push_back((GLfloat)v[2]);vertices.push_back(placeholder);}

inline void OpenGL_Color3(const GLfloat* color,Array<GLfloat>& colors)
{colors.push_back(color[0]);colors.push_back(color[1]);colors.push_back(color[2]);}

inline void OpenGL_Color4(const GLfloat* color,Array<GLfloat>& colors)
{colors.push_back(color[0]);colors.push_back(color[1]);colors.push_back(color[2]);colors.push_back(color[3]);}

inline void OpenGL_Vertex4_And_Color4(const Vector2& v,const GLfloat* c,Array<GLfloat>& vertices)
{OpenGL_Vertex4(v,vertices);OpenGL_Color4(c,vertices);}

inline void OpenGL_Vertex4_And_Color4(const Vector3& v,const GLfloat* c,Array<GLfloat>& vertices)
{OpenGL_Vertex4(v,vertices);OpenGL_Color4(c,vertices);}

inline void OpenGL_Color(const GLfloat* color,Array<GLfloat>& colors){OpenGL_Color4(color,colors);}

template<class T_ARRAY> inline void OpenGL_Quad(const T_ARRAY& nodes,Array<GLuint>& vertices)	////For triangle_strip
{int s[4]={0,2,1,3};for(int i=0;i<4;i++)vertices.push_back((GLuint)nodes[s[i]]);vertices.push_back(UINT_MAX);}

template<class T_ARRAY> inline void OpenGL_Cube(const T_ARRAY& nodes,Array<GLuint>& vertices)
{static const int s1[8]={5,7,1,3,0,2,4,6};static const int s2[8]={1,0,5,4,7,6,3,2};
for(int i=0;i<8;i++)vertices.push_back((GLuint)nodes[s1[i]]);vertices.push_back(UINT_MAX);
for(int i=0;i<8;i++)vertices.push_back((GLuint)nodes[s2[i]]);vertices.push_back(UINT_MAX);}

template<class T_ARRAY> inline void OpenGL_Cube_Wireframe(const T_ARRAY& nodes,Array<GLuint>& vertices)
{static const int s[24]={0,1,1,5,5,4,4,0,2,3,3,7,7,6,6,2,0,2,1,3,4,6,5,7};
for(int i=0;i<24;i++){vertices.push_back((GLuint)nodes[s[i]]);}}

template<class T_ARRAY> inline void OpenGL_Quad_Wireframe(const T_ARRAY& nodes,Array<GLuint>& vertices)
{static const int s[8]={0,1,1,3,3,2,2,0};
for(int i=0;i<8;i++){vertices.push_back((GLuint)nodes[s[i]]);}}

#endif
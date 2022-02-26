//#####################################################################
// Copyright 2002-2007, Robert Bridson, Geoffrey Irving, Igor Neverov.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
// Class OPENGL_COLOR
//##################################################################### 
#ifndef __OPENGL_COLOR__
#define __OPENGL_COLOR__

#include "Common.h"
#include "File.h"
#include <GL/glew.h>

class OPENGL_COLOR
{
public:
    GLfloat rgba[4];

    explicit OPENGL_COLOR(const float red=0,const float green=0,const float blue=0,const float alpha=1)
    {rgba[0]=red;rgba[1]=green;rgba[2]=blue;rgba[3]=alpha;}

	explicit OPENGL_COLOR(const ArrayF<float,4>& _rgba)
	{rgba[0]=_rgba[0];rgba[1]=_rgba[1];rgba[2]=_rgba[2];rgba[3]=_rgba[3];}

    explicit OPENGL_COLOR(const Vector3d& rgb,const float alpha=1)
    {rgba[0]=(float)rgb.x();rgba[1]=(float)rgb.y();rgba[2]=(float)rgb.z();rgba[3]=alpha;}

    explicit OPENGL_COLOR(const Vector3f& rgb,const float alpha=1)
    {rgba[0]=rgb.x();rgba[1]=rgb.y();rgba[2]=rgb.z();rgba[3]=alpha;}

    float Luminance() const
    {return (float)(.299*(float)rgba[0]+.587*(float)rgba[1]+.114*(float)rgba[2]);}

    OPENGL_COLOR Grayscale() const
    {float y=Luminance();return OPENGL_COLOR(y,y,y,rgba[3]);}

    void Send_To_GL_Pipeline() const
    {glColor4fv(rgba);}

    OPENGL_COLOR operator*(const double a)const
    { return OPENGL_COLOR((float)a*rgba[0],(float)a*rgba[1],(float)a*rgba[2],(float)a*rgba[3]); }

    OPENGL_COLOR& operator*=(const double a)
    { rgba[0]*=(GLfloat)a; rgba[1]*=(GLfloat)a; rgba[2]*=(GLfloat)a; rgba[3]*=(GLfloat)a; return *this;}

    OPENGL_COLOR operator+(const OPENGL_COLOR& c)const
    { return OPENGL_COLOR(rgba[0]+c.rgba[0],rgba[1]+c.rgba[1],rgba[2]+c.rgba[2],rgba[3]+c.rgba[3]); }

    bool operator==(const OPENGL_COLOR& c)const
    { return rgba[0]==c.rgba[0] && rgba[1]==c.rgba[1] && rgba[2]==c.rgba[2] && rgba[3]==c.rgba[3]; }

    static OPENGL_COLOR From_HSV(float h,const float s,const float v)
    {
        if(s==0) return OPENGL_COLOR(v,v,v);
        h/=360;h-=floor(h);h*=6; //h is now in [0,6)
        int i=(int)floor(h);
        float f=h-i,p=v*(1-s),q=v*(1-s*f),t=v*(1-s*(1-f));
        switch(i){
            case 0:return OPENGL_COLOR(v,t,p);
            case 1:return OPENGL_COLOR(q,v,p);
            case 2:return OPENGL_COLOR(p,v,t);
            case 3:return OPENGL_COLOR(p,q,v);
            case 4:return OPENGL_COLOR(t,p,v);
            case 5:return OPENGL_COLOR(v,p,q);}
    }

    template<class T>
    void To_HSV(T* h, T* s, T* v) const
    {   T r=rgba[0],g=rgba[1],b=rgba[2];
        T max_component=max(r,g,b), min_component=min(r,g,b), delta=max_component-min_component;
        *v=max_component;
        *s=(max_component!=0)?(delta/max_component):0;
        if(*s==0) h=0;//an arbitrary but sensible value
        else{
            if(r==max_component) *h=(g-b)/delta;
            else if(g==max_component) *h=2.0+(b-r)/delta;
            else *h=4.0+(r-g)/delta;
            (*h)*=60; if(*h<0)(*h)+=360.0;
        }
    }   

    static OPENGL_COLOR Random()
    {return OPENGL_COLOR(rand()/(float)RAND_MAX,rand()/(float)RAND_MAX,rand()/(float)RAND_MAX);}

    static OPENGL_COLOR Red(float r=1,float alpha=1) {return OPENGL_COLOR(r,0,0,alpha);}
    static OPENGL_COLOR Green(float g=1,float alpha=1) {return OPENGL_COLOR(0,g,0,alpha);}
    static OPENGL_COLOR Ground_Tan(float g=1,float alpha=1) {return OPENGL_COLOR(g,g*(float).775,g*(float).5431,alpha);}
    static OPENGL_COLOR Blue(float b=1,float alpha=1) {return OPENGL_COLOR(0,0,b,alpha);}
    static OPENGL_COLOR Yellow(float y=1,float alpha=1) {return OPENGL_COLOR(y,y,0,alpha);}
    static OPENGL_COLOR Cyan(float c=1,float alpha=1) {return OPENGL_COLOR(0,c,c,alpha);}
    static OPENGL_COLOR Magenta(float m=1,float alpha=1) {return OPENGL_COLOR(m,0,m,alpha);}
    static OPENGL_COLOR Gray(float y=.5,float alpha=1) {return OPENGL_COLOR(y,y,y,alpha);}
    static OPENGL_COLOR Orange(float o=1.,float alpha=1) {return OPENGL_COLOR(o,.64f,0.f,alpha);}
    static OPENGL_COLOR Bone_White(float r=.891,float g=.816,float b=.75,float alpha=1) {return OPENGL_COLOR(r,g,b,alpha);}
    static OPENGL_COLOR White(float alpha = 1) { return OPENGL_COLOR(1, 1, 1, alpha); }
    static OPENGL_COLOR Black() {return OPENGL_COLOR(0,0,0);}
    static OPENGL_COLOR Transparent() {return OPENGL_COLOR(0,0,0,0);}

    static OPENGL_COLOR Affine_Combination(float t,const OPENGL_COLOR& color1,const OPENGL_COLOR& color2)
    {return OPENGL_COLOR(color1.rgba[0]+t*(color2.rgba[0]-color1.rgba[0]),color1.rgba[1]+t*(color2.rgba[1]-color1.rgba[1]),
                            color1.rgba[2]+t*(color2.rgba[2]-color1.rgba[2]),color1.rgba[3]+t*(color2.rgba[3]-color1.rgba[3]));}

    template<class RW>
    void Read(std::istream& input_stream)
    {for(int i=0;i<4;i++){float tmp;File::Read_Binary<float>(input_stream,tmp);rgba[i]=tmp;}} // Read as float

    template<class RW>
    void Write(std::ostream& output_stream) const
    {for(int i=0;i<4;i++){File::Write_Binary<float>(output_stream,(float)rgba[i]);}}

//#####################################################################
    friend OPENGL_COLOR operator*(const double a,const OPENGL_COLOR& c)
    {return OPENGL_COLOR((float)a*c.rgba[0],(float)a*c.rgba[1],(float)a*c.rgba[2],(float)a*c.rgba[3]);}
};

using OpenGLColor=OPENGL_COLOR;
#endif
//#####################################################################
// Copyright 2008, Michael Lentine.
// This file is part of PhysBAM whose distribution is governed by the license contained in the accompanying file PHYSBAM_COPYRIGHT.txt.
//#####################################################################
#ifndef __H_OPENGL_ARCBALL__
#define __H_OPENGL_ARCBALL__
#include "Common.h"
#include "AuxFunc.h"

class OpenGLWindow;

class SPHERE {public:Vector3 center;float radius;SPHERE():center(Vector3::Zero()),radius(1.f){}};

class OpenGLArcball
{
    typedef Vector<float,3> TV;typedef Vector<float,2> TV2;
public:
    SPHERE sphere;
    Quaternionf qNow,qDown,qDrag;
    TV2 center,vDown;
    TV vFrom,vTo,vrFrom,vrTo;
    Matrix4f mNow,mDown,mDeltaNow; //note that these only work in world space
    bool dragging;
    OpenGLWindow* window;

    OpenGLArcball():dragging(false)
    {Reinitialize();}

    void Reinitialize()
    {sphere=SPHERE();
    qNow=qDown=qDrag=Quaternionf(1,0,0,0);
    center=vDown=TV2::Zero();
    vFrom=vTo=vrFrom=vrTo=TV::Zero();
    mNow=mDown=mDeltaNow=Matrix4f::Identity();}

    void Update(const TV2& vNow){Update_World(vNow);}

    Matrix4f From_Linear(const Matrix3f& M) // Create a homogeneous 4x4 matrix corresponding to a 3x3 transform
    {Matrix4f mat; mat<<M.data()[0],M.data()[3],M.data()[6],0,M.data()[1],M.data()[4],M.data()[7],0,M.data()[2],M.data()[5],M.data()[8],0,0,0,0,1;return mat;}

    void Update_World(const TV2 &vNow)
    {vFrom=MouseOnSphere(vDown,center,sphere.radius);vTo=MouseOnSphere(vNow,center,sphere.radius);
    if (dragging){qDrag=Qt_FromBallPoints(vFrom,vTo);qNow=qDrag*qDown;}
    Qt_ToBallPoints(qDown,vrFrom,vrTo);
    mNow=From_Linear(qNow.toRotationMatrix());
    mDeltaNow=From_Linear(qDrag.toRotationMatrix());}

    Matrix4f Value(){return mDeltaNow;}
    void Begin_Drag(const TV2 &vNow){dragging=true;vDown=vNow;Update(vNow);}
    void End_Drag(const TV2 &vNow){dragging=false;Update(vNow);qDown=qNow;mDown=mNow;}
    
private:
    TV MouseOnSphere(const TV2 &mouse,const TV2 &ballCenter,double ballRadius)
    {TV ballMouse=TV::Zero();float mag=(float)0;
    ballMouse.x()=float((mouse.x()-ballCenter.x())/ballRadius);
    ballMouse.y()=float((mouse.y()-ballCenter.y())/ballRadius);
    mag=ballMouse.squaredNorm();
    if (mag>1.0){float scale=float(1.0/sqrt(mag));ballMouse.x()*=scale;ballMouse.y()*=scale;ballMouse.z()=0.0;}
    else ballMouse.z()=float(sqrt(1-mag));
    return ballMouse;}

    Quaternionf Qt_FromBallPoints(const TV &from,const TV &to)
    {Quaternionf q;q.setFromTwoVectors(from,to);return q;}

    Vector3f Orthogonal_Vector(const Vector3f& v)
    {
        float abs_x=abs(v.x()),abs_y=abs(v.y()),abs_z=abs(v.z());
        if(abs_x<abs_y) return abs_x<abs_z?Vector3f((float)0,v.z(),-v.y()):Vector3f(v.y(),-v.x(),(float)0);
        else return abs_y<abs_z?Vector3f(-v.z(),(float)0,v.x()):Vector3f(v.y(),-v.x(),(float)0);
    }

    void Qt_ToBallPoints(const Quaternionf &q,TV &arcFrom,TV &arcTo)
    //{arcFrom=q.Get_Axis().Unit_Orthogonal_Vector();arcTo=q.Rotate(arcFrom);}
    {TV vec=(q.w()>=0?(float)1:(float)-1)*q.vec().normalized();arcFrom=Orthogonal_Vector(vec).normalized();arcTo=Transform3f(q)*arcFrom;}
};
#endif

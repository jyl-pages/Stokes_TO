//#####################################################################
// PhysBAM Adaptor
// Bo Zhu
//#####################################################################
#ifndef __PhysBAMAdaptor_h__
#define __PhysBAMAdaptor_h__
#include "Common.h"
////physbam class
template<class real> class PLANE
{typedef Vector3 TV;public:TV normal;TV x1;PLANE(const TV& normal_input=TV(0,1,0),const TV& x1_input=TV(0,0,0)):normal(normal_input),x1(x1_input){}};

template<class TV> class RANGE
{public:
    TV min_corner,max_corner;
    RANGE(const TV& _min=TV::Zero(),const TV& _max=TV::Zero()):min_corner(_min),max_corner(_max){}
    TV Center() const {return (real).5*(min_corner+max_corner);}
    TV Edge_Lengths() const {return max_corner-min_corner;}
	bool Inside(const TV& pos)
	{for(int i=0;i<pos.size();i++)if(pos[i]<min_corner[i]||pos[i]>max_corner[i])return false;return true;}
};

////physbam helper functions
template<class TV> TV Projected(const TV& vec_from,const TV& vec_to){return vec_from.dot(vec_to)/vec_to.squaredNorm()*vec_to;}

//inline std::string string_sprintf(const char *format_string,...) // Assumes a max string length of 2048, since Windows doesn't support safety
//{char tmp[2048];va_list marker;va_start(marker,format_string);vsprintf_s(tmp,format_string,marker);va_end(marker);return tmp;}

inline Matrix4 Rotation_X_4(const real radians)
{real c=cos(radians),s=sin(radians);Matrix4 m;m<<1,0,0,0, 0,c,-s,0, 0,s,c,0, 0,0,0,1;return m;}

inline Matrix4 Rotation_Y_4(const real radians)
{real c=cos(radians),s=sin(radians);Matrix4 m;m<<c,0,s,0, 0,1,0,0, -s,0,c,0, 0,0,0,1;return m;}

inline Matrix4 Rotation_Z_4(const real radians)
{real c=cos(radians),s=sin(radians);Matrix4 m;m<<c,-s,0,0, s,c,0,0, 0,0,1,0, 0,0,0,1;return m;}

inline Vector3 Homogeneous_Times(const Matrix4& m,const Vector3& v) // assumes w=1 is the 4th coordinate of v
{const real* x=m.data();real w=x[3]*v.x()+x[7]*v.y()+x[11]*v.z()+x[15];assert(w!=0);
if(w==1) return Vector3(x[0]*v.x()+x[4]*v.y()+x[8]*v.z()+x[12],x[1]*v.x()+x[5]*v.y()+x[9]*v.z()+x[13],x[2]*v.x()+x[6]*v.y()+x[10]*v.z()+x[14]);
real s=1/w;// rescale so w=1
return Vector3(s*(x[0]*v.x()+x[4]*v.y()+x[8]*v.z()+x[12]),s*(x[1]*v.x()+x[5]*v.y()+x[9]*v.z()+x[13]),s*(x[2]*v.x()+x[6]*v.y()+x[10]*v.z()+x[14]));}

inline int Dominant_Axis(const Vector3& v)
{real abs_x=abs(v.x()),abs_y=abs(v.y()),abs_z=abs(v.z());return (abs_x>abs_y && abs_x>abs_z)?1:((abs_y>abs_z)?2:3);}

////physbam type macros
template<class T1,class T2> using PAIR=std::pair<T1,T2>;
template<class T1,class T2> using Pair=std::pair<T1,T2>;
template<class T1> using ARRAY=Array<T1>;
using TV2=Vector2;
using TV3=Vector3;

#endif
//////////////////////////////////////////////////////////////////////////
// Geometry primitives
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __GeometryPrimitives_h__
#define __GeometryPrimitives_h__
#include <limits>
#include <iostream>
#include "Common.h"
#include "AuxFunc.h"
#include "Constants.h"
#include "MeshFunc.h"

//////////////////////////////////////////////////////////////////////////
////Base class
//////////////////////////////////////////////////////////////////////////
template<int d> class ImplicitGeometry
{Typedef_VectorDii(d);
public:
	ImplicitGeometry(){}
	virtual real Phi(const VectorD& pos) const = 0;
	virtual bool Inside(const VectorD& pos) const {return Phi(pos)<(real)0;}
	virtual VectorD Normal(const VectorD& pos) const = 0;
};

//////////////////////////////////////////////////////////////////////////
////Dimension-independent geometries
//////////////////////////////////////////////////////////////////////////

template<int d> class Bowl : public ImplicitGeometry<d>
{
Typedef_VectorDii(d);
public:
	VectorD center;
	real radius;

	Bowl(VectorD _center=VectorD::Zero(),real _radius=1.):center(_center),radius(_radius){}
	Bowl<d>& operator=(const Bowl<d>& copy){center=copy.center;radius=copy.radius;return *this;}
	Bowl(const Bowl<d>& copy){*this=copy;}

	virtual real Phi(const VectorD& pos) const {return radius-(pos-center).norm();}
	virtual VectorD Normal(const VectorD& pos) const {return (center-pos).normalized();}
};

template<int d> class HalfBowl : public ImplicitGeometry<d>
{
	Typedef_VectorDii(d);
public:
	VectorD center;
	real radius;
	VectorD bowl_normal;//pointing toward the "bowl side" half-space

	HalfBowl(VectorD _center = VectorD::Zero(), real _radius = 1., VectorD _bowl_normal = VectorD::Unit(1)) : center(_center), radius(_radius), bowl_normal(_bowl_normal) {
		bowl_normal.normalize();
	}
	HalfBowl<d>& operator=(const HalfBowl<d>& copy) { center = copy.center;radius = copy.radius;bowl_normal = copy.bowl_normal;return *this; }
	HalfBowl(const HalfBowl<d>& copy) { *this = copy; }

	VectorD Nearest_On_Edge(const VectorD& diff) const {
		VectorD offset = diff - diff.dot(bowl_normal) * bowl_normal;
		offset.normalize();
		if (offset.norm() == 0) offset = AuxFunc::Orthogonal_Vector(bowl_normal).normalized();
		return offset * radius + center;
	}
	virtual real Phi(const VectorD& pos) const { 
		VectorD diff = (pos - center);
		if (diff.dot(bowl_normal) >= 0) {
			return radius - diff.norm();
		}
		else {
			return (pos - Nearest_On_Edge(diff)).norm();
		}
	}
	virtual VectorD Normal(const VectorD& pos) const {
		VectorD diff = (pos - center);
		if (diff.dot(bowl_normal) >= 0) {
			return (-diff).normalized();
		}
		else {
			return (pos - Nearest_On_Edge(diff)).normalized();
		}
	}
};

template<int d> class Tube : public ImplicitGeometry<d>
{
	//It's a dedicated-use geometry. It may not be what you think it is. --Mengdi Wang
	Typedef_VectorDii(d);
public:
	VectorD center = VectorD::Zero();
	real radius = (real)1;
	VectorD normal = VectorD::Unit(0);
	real height = (real)1;
	Tube(const VectorD& _center, const real _radius, const VectorD& _normal, const real _height) :center(_center), radius(_radius), normal(_normal), height(_height) {
		normal.normalize();
	}
	virtual real Phi(const VectorD& pos) const {
		VectorD diff = pos - center;
		real diff_dot_normal = diff.dot(normal);
		if (fabs(diff_dot_normal) <= height / 2) {
			VectorD plane_diff = diff - diff_dot_normal * normal;
			return radius - plane_diff.norm();
		}
		else return radius;
	}
	virtual VectorD Normal(const VectorD& pos) const {
		VectorD diff = pos - center;
		VectorD plane_diff = diff - diff.dot(normal) * normal;
		return -plane_diff.normalized();
	}
};

template<int d> class Sphere : public ImplicitGeometry<d>
{Typedef_VectorDii(d);
public:
	VectorD center=VectorD::Zero();
	real radius=(real)1;
	Sphere(const VectorD& _center,const real _radius):center(_center),radius(_radius){}
	Sphere<d>& operator=(const Sphere<d>& copy){center=copy.center;radius=copy.radius;return *this;}
	Sphere(const Sphere<d>& copy){*this=copy;}

	virtual real Phi(const VectorD& pos) const {return (pos-center).norm()-radius;}
	virtual VectorD Normal(const VectorD& pos) const {return -(pos-center).normalized();}
};

template<int d> class Ellipsoid : public ImplicitGeometry<d>
{Typedef_VectorDii(d);
public:
	VectorD center;VectorD radius;
	Ellipsoid(const VectorD& _center,const VectorD& _radius):center(_center),radius(_radius){}
	Ellipsoid<d>& operator=(const Ellipsoid<d>& copy){center=copy.center;radius=copy.radius;return *this;}
	Ellipsoid(const Ellipsoid<d>& copy){*this=copy;}

	virtual real Phi(const VectorD& pos) const     ////approximate signed distance
	{VectorD p0=pos-center;VectorD scaled_p0=p0.cwiseQuotient(radius);real length=scaled_p0.norm();scaled_p0/=length;
	real d0=(p0-radius.cwiseProduct(scaled_p0)).norm();return length<(real)1?-d0:d0;}

	virtual VectorD Normal(const VectorD& pos) const {return (pos-center).normalized();}	////not accurate
};

template<int d> class Box : public ImplicitGeometry<d>
{Typedef_VectorDii(d);
public:
	VectorD min_corner,max_corner;

	Box(const VectorD& _min=VectorD::Zero(),const VectorD& _max=VectorD::Zero()):min_corner(_min),max_corner(_max){}
	Box<d>& operator=(const Box<d>& copy){min_corner=copy.min_corner;max_corner=copy.max_corner;return *this;}
	Box(const Box<d>& copy){*this=copy;}

	virtual bool Inside(const VectorD& pos) const {return AuxFunc::All_Greater_Equal(pos,min_corner)&&AuxFunc::All_Less_Equal(pos,max_corner);}
	virtual real Phi(const VectorD& pos) const 
	{VectorD phi=(pos-Center()).cwiseAbs()-(real).5*Edge_Lengths();VectorD zero=VectorD::Zero();
	if(!AuxFunc::All_Less_Equal(phi,zero)) return (phi.cwiseMax(zero)).norm();return phi.maxCoeff();}
	virtual VectorD Normal(const VectorD& pos)const { return VectorD::Zero(); }
	VectorD Edge_Lengths() const {return max_corner-min_corner;}
	VectorD Center() const {return (real).5*(min_corner+max_corner);}
	Box<d> Enlarged(const Box<d>& box2) const {return Box<d>(AuxFunc::Cwise_Min(min_corner,box2.min_corner),AuxFunc::Cwise_Max(max_corner,box2.max_corner));}
	Box<d> Enlarged(const VectorD& length) const {return Box<d>(min_corner-length,max_corner+length);}
	Box<d> Rescaled(const real factor) const {VectorD length=Edge_Lengths();return Box<d>(min_corner-length*factor*(real).5,max_corner+length*factor*(real).5);}
	VectorD Wall_Normal(const VectorD& pos) const
	{VectorD normal=VectorD::Zero();
	for(int i=0;i<d;i++){
		if(pos[i]<min_corner[i])normal[i]=min_corner[i]-pos[i];
		else if(pos[i]>max_corner[i])normal[i]=max_corner[i]-pos[i];}
	if(normal!=VectorD::Zero())return normal.normalized();return VectorD::Zero();}
	static Box<d> Infi_Min(){const real fmax=std::numeric_limits<real>::max();return Box<d>(Vector<real,d>::Ones()*fmax,Vector<real,d>::Ones()*(real)-fmax);}
};

template<int d> class NegativeBox : public Box<d>
{Typedef_VectorDii(d);using Base=Box<d>;
public:
	using Base::min_corner;using Base::max_corner;
	NegativeBox(const VectorD& _min=VectorD::Zero(),const VectorD& _max=VectorD::Zero()):Base(_min,_max){}
	NegativeBox<d>& operator=(const NegativeBox<d>& copy){min_corner=copy.min_corner;max_corner=copy.max_corner;return *this;}
	NegativeBox(const NegativeBox<d>& copy){*this=copy;}

	virtual bool Inside(const VectorD& pos) const {return !Base::Inside(pos);}
	virtual real Phi(const VectorD& pos) const {return -Base::Phi(pos);}
	VectorD Wall_Normal(const VectorD& pos) const {return -Base::Wall_Normal(pos);}
};

//////////////////////////////////////////////////////////////////////////
////Plane
//////////////////////////////////////////////////////////////////////////

template<int d> class Plane : public ImplicitGeometry<d>
{Typedef_VectorDii(d);
public:
	VectorD n;
	VectorD p;
	real b;

	Plane(const VectorD _n,const VectorD _p):n(_n),p(_p){n.normalize();b=n.dot(p);}
	Plane<d>& operator=(const Plane<d>& copy){n=copy.n;p=copy.p;b=copy.b;return *this;}
	Plane(const Plane<d>& copy){*this=copy;}

	virtual bool Inside(const VectorD& pos) const {return n.dot(pos)-b<(real)0;}
	virtual real Phi(const VectorD& pos) const {return (n.dot(pos)-b);}
	virtual VectorD Normal(const VectorD& pos) const { return n; }
};

template<> class Plane<3> : public ImplicitGeometry<3>
{Typedef_VectorDii(3);
public:
	VectorD n;
	VectorD p;
	real b;

	Plane(const VectorD _n, const VectorD _p) :n(_n), p(_p) { n.normalize(); b = n.dot(p); }
	Plane(const VectorD& p0,const VectorD& p1,const VectorD& p2){n=((p1-p0).cross(p2-p1)).normalized();p=p0;b=n.dot(p);}
	Plane<3>& operator=(const Plane<3>& copy){n=copy.n;p=copy.p;b=copy.b;return *this;}
	Plane(const Plane<3>& copy){*this=copy;}

	virtual bool Inside(const VectorD& pos) const {return n.dot(pos)-b<(real)0;}
	virtual real Phi(const VectorD& pos) const {return (n.dot(pos)-b);}
	virtual VectorD Normal(const VectorD& pos) const {return n;}

	static bool Inside(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{VectorD n=((p1-p0).cross(p2-p0));real bp=n.dot(p0);real bq=n.dot(pos);return bq<bp;}

	static VectorD Projection(const VectorD& p0,const VectorD& p1,const VectorD& p2,const VectorD& pos)
	{VectorD n=((p1-p0).cross(p2-p0)).normalized();VectorD p0q=pos-p0;VectorD proj_p0q=p0q.dot(n)*n;return pos-proj_p0q;}
};

//////////////////////////////////////////////////////////////////////////
////Line
//////////////////////////////////////////////////////////////////////////

template<int d> class Line{};

template<> class Line<2>: public ImplicitGeometry<2>
{Typedef_VectorDii(2);
public:
	VectorD p0;
	VectorD p1;

	Line(const VectorD& _p0,const VectorD& _p1):p0(_p0),p1(_p1){}

	virtual bool Inside(const VectorD& pos) const {return Inside(p0,p1,pos);}

	virtual real Phi(const VectorD& pos) const {return Phi(p0,p1,pos);}

	virtual VectorD Normal(const VectorD& pos) const {return Normal(p0,p1,pos);}

	static bool Inside(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{
		VectorD p0q=pos-p0;VectorD p01=p1-p0;
		return AuxFunc::Cross(p0q,p01)[0]<0;		
	}

	static real Phi(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{
		real coef=Inside(p0,p1,pos)?(real)-1:(real)1;
		return coef*Closest_Point_Vector(p0,p1,pos).norm();
	}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{
		real coef=Inside(p0,p1,pos)?(real)-1:(real)1;
		return coef*Closest_Point_Vector(p0,p1,pos).normalized();
	}

	//// All the Closest_Point_Vector functions (the same as below) return a vector 
	//// pointing *FROM" the closest point on the interface to the starting point
	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{
		VectorD n=(p1-p0).normalized();
		VectorD pq=pos-p0;
		return (pq-pq.dot(n)*n);
	}
};

template<> class Line<3>: public ImplicitGeometry<3>
{Typedef_VectorDii(3);
public:
	VectorD p0;
	VectorD p1;

	Line(const VectorD& _p0,const VectorD& _p1):p0(_p0),p1(_p1){}

	virtual bool Inside(const VectorD& pos) const {return false;}	////always outside	
	virtual real Phi(const VectorD& pos) const {return Phi(p0,p1,pos);}
	virtual VectorD Normal(const VectorD& pos) const {return Normal(p0,p1,pos);}

	static real Phi(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{VectorD v=Closest_Point_Vector(p0,p1,pos);return v.norm();}

	static VectorD Normal(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{VectorD v=Closest_Point_Vector(p0,p1,pos);return v.normalized();}

	static VectorD Closest_Point_Vector(const VectorD& p0,const VectorD& p1,const VectorD& pos)
	{
		VectorD n=(p1-p0).normalized();
		VectorD pq=pos-p0;
		return (pq-pq.dot(n)*n);
	}
};

//////////////////////////////////////////////////////////////////////////
////Torus
//////////////////////////////////////////////////////////////////////////

template<int d> class Torus{};

template<> class Torus<2> : public ImplicitGeometry<2>
{Typedef_VectorDii(2);
public:
	VectorD center=VectorD::Zero();
	Vector2d radius;

	Torus(const VectorD& _center, const Vector2d& _radius):center(_center),radius(_radius) {}
	Torus(const VectorD& _center, const Vector2d& _radius, const VectorD& _normal):center(_center),radius(_radius) {}
	Torus<2>& operator=(const Torus<2>& copy) {center=copy.center;radius=copy.radius;return *this;}
	Torus(const Torus<2>& copy) {*this=copy;}

	virtual bool Inside(const VectorD& pos) const {real dist=(pos-center).norm();return std::abs(dist-radius[0])<radius[1];}
	virtual real Phi(const VectorD& pos) const {real dist=(pos-center).norm();return std::abs(dist-radius[0])-radius[1];}
	virtual VectorD Normal(const VectorD& pos) const {return -(pos-center).normalized();}
};

template<> class Torus<3> : public ImplicitGeometry<3>
{Typedef_VectorDii(3);
public:
	VectorD center=VectorD::Zero();
	Vector2d radius;
	VectorD normal=VectorD::Unit(1);

	Torus(const VectorD& _center, const Vector2d& _radius, const VectorD& _normal):center(_center),radius(_radius),normal(_normal) {normal.normalized();}
	Torus<3>& operator=(const Torus<3>& copy) {center=copy.center;radius=copy.radius;normal=copy.normal;return *this;}
	Torus(const Torus<3>& copy) {*this=copy;}

	virtual bool Inside(const VectorD& pos) const {return Phi(pos)<=0;}
	virtual real Phi(const VectorD& pos) const 
	{
		VectorD v=pos-center;
		real dist=v.dot(normal);
		VectorD proj_pos=pos-dist*normal;
		Vector2d q((proj_pos-center).norm()-radius[0],dist);
		return q.norm()-radius[1];
	}
	virtual VectorD Normal(const VectorD& pos) const
	{
		VectorD v=pos-center;
		real dist=v.dot(normal);
		VectorD proj_pos=pos-dist*normal;
		Vector2d q((proj_pos-center).norm()-radius[0],dist);
		return (((proj_pos-center).normalized())*q(0)+dist*normal).normalized();
	}
};

//////////////////////////////////////////////////////////////////////////
////Ray
//////////////////////////////////////////////////////////////////////////

template<int d> class Ray
{Typedef_VectorDii(d);
public:
	VectorD point;VectorD direction;
	Ray(const VectorD& _p,const VectorD& _d):point(_p),direction(_d){}
};

template<int d> bool Intersection(const Ray<d>& ray,const ImplicitGeometry<d>& object,const real& t_min,const real& t_max,const real& step,/*rst*/real& t)
{Typedef_VectorDii(d);
	int n=(int)ceil((t_max-t_min)/step);
	for(int i=0;i<n;i++){
		real t0=t_min+(real)i*step;VectorD p0=ray.point+ray.direction*t0;real phi_0=object.Phi(p0);
		real t1=t_min+(real)(i+1)*step;VectorD p1=ray.point+ray.direction*t1;real phi_1=object.Phi(p1);
		if((phi_0<(real)0)!=(phi_1<(real)0)){real frac=abs(phi_0)/(abs(phi_0)+abs(phi_1));t=t0+frac*step;return true;}}
	return false;
}

#endif
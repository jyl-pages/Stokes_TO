#ifndef __AnalyticalFields_h__
#define __AnalyticalFields_h__
#include "Common.h"
#include "Constants.h"

//////////////////////////////////////////////////////////////////////////
////Analytical velocity fields

namespace AnaField
{
	Vector2 Taylor_Green_Velocity(const Vector2& pos)
	{
		real theta=(real)0;
		real u=(real)2/(real)sqrt((real)3)*sin(theta+two_thirds*pi)*sin(pos[0])*cos(pos[1]);
		real v=(real)2/(real)sqrt((real)3)*sin(theta-two_thirds*pi)*cos(pos[0])*sin(pos[1]);
		return Vector2(u,v);
	}

	Vector3 Taylor_Green_Velocity(const Vector3& pos)
	{
		real theta=(real)0;real coef=(real)2/(real)sqrt((real)3);
		real u=coef*sin(theta+two_thirds*pi)*sin(pos[0])*cos(pos[1])*cos(pos[2]);
		real v=coef*sin(theta-two_thirds*pi)*cos(pos[0])*sin(pos[1])*cos(pos[2]);
		real w=coef*sin(theta)*cos(pos[0])*cos(pos[1])*sin(pos[2]);
		return Vector3(u,v,w);
	}

	real Taylor_Green_Vorticity(const Vector2& pos)
	{
		real theta=(real)0;
		real wz=(real)2*cos(theta)*sin(pos[0])*sin(pos[1])*cos(pos[2]);
		return wz;
	}

	Vector3 Taylor_Green_Vorticity(const Vector3& pos)
	{
		real theta=(real)0;
		real wx=-(sqrt((real)3)*sin(theta)+cos(theta))*cos(pos[0])*sin(pos[1])*sin(pos[2]);
		real wy=(sqrt((real)3)*sin(theta)-cos(theta))*sin(pos[0])*cos(pos[1])*sin(pos[2]);
		real wz=(real)2*cos(theta)*sin(pos[0])*sin(pos[1])*cos(pos[2]);
		return Vector3(wx,wy,wz);	
	}
};

template<int d> class AnalyticalField
{Typedef_VectorDii(d);
public:
	virtual VectorD Velocity(const VectorD& pos,const real t=(real)0) const {return VectorD::Zero();}
};

template<int d> class VorticityField : public AnalyticalField<d>
{Typedef_VectorDii(d);
public:
	VectorD origin=VectorD::Zero();
	real strength=(real)1;
	int test=0;
	real period=(real)4;

	virtual VectorD Velocity(const VectorD& pos,const real time=(real)0) const 
	{
		switch(test){
		case 0:return VectorD::Zero();
		case 1:return Rigid_Rotation(pos);
		case 2:{	////steady-state vortex
			return Vortex_In_Box(pos,strength);}
		case 3:{	////temporal vortex
			real s=strength*sin(two_pi*time/period);
			return Vortex_In_Box(pos,s);}
		case 4:{	////translation
			return VectorD::Unit(0)*(real).2;}
		case 5:{
				
		}
		}
		return VectorD::Zero();
	}

	VectorD Rigid_Rotation(const VectorD& pos) const 
	{return VectorD::Unit(0)*strength*(origin[1]-pos[1])+VectorD::Unit(1)*strength*(pos[0]-origin[0]);}

	VectorD Vortex_In_Box(const VectorD& _pos,real strength) const 
	{
		VectorD pos=_pos-origin;real c=pi*(real).5;
		real ux=strength*c*pow(sin(c*pos[0]),2)*sin(c*pos[1])*cos(c*pos[1]);
		real uy=-strength*c*pow(sin(c*pos[1]),2)*sin(c*pos[0])*cos(c*pos[0]);
		return VectorD::Unit(0)*ux+VectorD::Unit(1)*uy;
	}
};

#endif
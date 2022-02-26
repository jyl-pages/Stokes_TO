//////////////////////////////////////////////////////////////////////////
// Random number
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Random_Number_h__
#define __Random_Number_h__
#include <random>
#include <chrono>
#include "Common.h"
#include "Constants.h"

namespace RandomFunc {
	real Random_Real(real a, real b);
	template<int d> Vector<real, d> Random_Vector_Cartesian(real a, real b);//just every component is in [a,b]
	template<int d> Vector<real, d> Random_Vector_Spherical(real max_r);//the norm of all components are <=max_r
	bool Sample_With_Probability(real prob);//return true with prob and false with 1-prob
	int Stochastic_Coarsed_Real(real val);//like, if val=3.7, return 3 with probability of 0.3 and 4 with probability 0.7
}

class RandomNumber
{
public:
	RandomNumber(real r_min=(real)0,real r_max=(real)1,bool time_dependent=false):dis(r_min,r_max){
		if (time_dependent) { unsigned seed = std::chrono::steady_clock::now().time_since_epoch().count(); gen.seed(seed); }
		else { gen.seed(); }
	}
	real Value(){return dis(gen);}
	template<int d> Vector<real,d> VectorValue(){Vector<real,d> v;for(int i=0;i<d;i++)v[i]=Value();return v;}

protected:
	 std::mt19937 gen;
	 real r_min,r_max;
	 std::uniform_real_distribution<real> dis;
};

class RandomInt {
public:
	std::mt19937 gen;
	std::uniform_int_distribution<int> uid;
	RandomInt(int _min, int _max) :uid(_min, _max) { gen.seed(); }
	int Value(void) { return uid(gen); }
};

template<int d> class RandomGeometry
{
	Typedef_VectorDii(d);
public:
	RandomNumber random;

	RandomGeometry(bool time_dependent = false) :random((real)0, (real)1, time_dependent) {}
	void Triangle(const real min_r,const real max_r,ArrayF<VectorD,3>& tri_vtx)
	{
		real start=random.Value()*two_thirds*pi;
		for(int i=0;i<3;i++){
			real min_theta=(real)i*two_thirds*pi+start;
			real theta=random.Value()*(real).2*pi+min_theta;
			real r=random.Value()*(max_r-min_r)+min_r;
			tri_vtx[i]=(VectorD::Unit(0)*cos(theta)+VectorD::Unit(1)*sin(theta))*r;}
	}

	VectorD Point(const VectorD& box_min,const VectorD& box_max)
	{VectorD point=box_min;for(int i=0;i<d;i++)point[i]+=random.Value()*(box_max[i]-box_min[i]);return point;}

	void Triangle(const real min_r,const real max_r,const VectorD& center,ArrayF<VectorD,3>& tri_vtx)
	{Triangle(min_r,max_r,tri_vtx);for(int i=0;i<3;i++)tri_vtx[i]+=center;}

};


#endif
//////////////////////////////////////////////////////////////////////////
// Interpolation
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Interpolation_h__
#define __Interpolation_h__
#include "Common.h"
#include "Constants.h"
#include "Grid.h"
#include "MacGrid.h"
#include "Field.h"
#include "FaceField.h"
#include "Constants.h"
#include "AuxFunc.h"

namespace IntpFunc
{
	//////////////////////////////////////////////////////////////////////////
	////Linear interpolation
	template<class T> T Nodes_To_Point_Linear(const Field<T,1>& field,const Vector1i& cell,const Vector1& frac)
	{return ((real)1-frac[0])*field(cell)+frac[0]*field(Vec1i(cell[0]+1));}

	template<class T> T Nodes_To_Point_Linear(const Field<T,2>& field,const Vector2i& cell,const Vector2& frac)
	{
		return ((real)1-frac[0])*((real)1-frac[1])*field(cell)
			+frac[0]*((real)1-frac[1])*field(cell[0]+1,cell[1])
			+((real)1-frac[0])*frac[1]*field(cell[0],cell[1]+1)
			+frac[0]*frac[1]*field(cell[0]+1,cell[1]+1);
	}

	template<class T> T Nodes_To_Point_Linear(const Field<T,3>& field,const Vector3i& cell,const Vector3& frac)
	{	
		return ((real)1-frac[0])*((real)1-frac[1])*((real)1-frac[2])*field(cell[0],cell[1],cell[2])
			+(frac[0])*((real)1-frac[1])*((real)1-frac[2])*field(cell[0]+1,cell[1],cell[2])
			+((real)1-frac[0])*(frac[1])*((real)1-frac[2])*field(cell[0],cell[1]+1,cell[2])
			+((real)1-frac[0])*((real)1-frac[1])*(frac[2])*field(cell[0],cell[1],cell[2]+1)
			+(frac[0])*(frac[1])*((real)1-frac[2])*field(cell[0]+1,cell[1]+1,cell[2])
			+((real)1-frac[0])*(frac[1])*(frac[2])*field(cell[0],cell[1]+1,cell[2]+1)
			+(frac[0])*((real)1-frac[1])*(frac[2])*field(cell[0]+1,cell[1],cell[2]+1)
			+(frac[0])*(frac[1])*(frac[2])*field(cell[0]+1,cell[1]+1,cell[2]+1);
	}

	template<class T> void Point_To_Nodes_Linear(const T& value,Field<T,2>& field,Field<real,2>& weight,const Vector2i& cell,const Vector2& frac)
	{
		real w00=((real)1-frac[0])*((real)1-frac[1]);
		real w10=frac[0]*((real)1-frac[1]);
		real w01=((real)1-frac[0])*frac[1];
		real w11=frac[0]*frac[1];

		field(cell[0],cell[1])+=w00*value;
		field(cell[0]+1,cell[1])+=w10*value;
		field(cell[0],cell[1]+1)+=w01*value;
		field(cell[0]+1,cell[1]+1)+=w11*value;

		weight(cell[0],cell[1])+=w00;
		weight(cell[0]+1,cell[1])+=w10;
		weight(cell[0],cell[1]+1)+=w01;
		weight(cell[0]+1,cell[1]+1)+=w11;
	}

	template<class T> void Point_To_Nodes_Linear(const T& value,Field<T,3>& field,Field<real,3>& weight,const Vector3i& cell,const Vector3& frac)
	{
		real w000=((real)1-frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		real w100=(frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		real w010=((real)1-frac[0])*(frac[1])*((real)1-frac[2]);
		real w001=((real)1-frac[0])*((real)1-frac[1])*(frac[2]);
		real w110=(frac[0])*(frac[1])*((real)1-frac[2]);
		real w011=((real)1-frac[0])*(frac[1])*(frac[2]);
		real w101=(frac[0])*((real)1-frac[1])*(frac[2]);
		real w111=(frac[0])*(frac[1])*(frac[2]);

		field(cell[0],cell[1],cell[2])+=w000*value;
		field(cell[0]+1,cell[1],cell[2])+=w100*value;
		field(cell[0],cell[1]+1,cell[2])+=w010*value;
		field(cell[0],cell[1],cell[2]+1)+=w001*value;
		field(cell[0]+1,cell[1]+1,cell[2])+=w110*value;
		field(cell[0],cell[1]+1,cell[2]+1)+=w011*value;
		field(cell[0]+1,cell[1],cell[2]+1)+=w101*value;
		field(cell[0]+1,cell[1]+1,cell[2]+1)+=w111*value;

		weight(cell[0],cell[1],cell[2])+=w000;
		weight(cell[0]+1,cell[1],cell[2])+=w100;
		weight(cell[0],cell[1]+1,cell[2])+=w010;
		weight(cell[0],cell[1],cell[2]+1)+=w001;
		weight(cell[0]+1,cell[1]+1,cell[2])+=w110;
		weight(cell[0],cell[1]+1,cell[2]+1)+=w011;
		weight(cell[0]+1,cell[1],cell[2]+1)+=w101;
		weight(cell[0]+1,cell[1]+1,cell[2]+1)+=w111;
	}

	//////////////////////////////////////////////////////////////////////////
	////Quadratic interpolation
	inline real Quadratic_Kernel(const real abs_x)
	{
		if(abs_x<(real).5) return (real).75-pow(abs_x,2);
		else if(abs_x<(real)1.5) return (real).5*pow((real)1.5-abs_x,2);
		else return (real)0;
	}

	template<class T> T Nodes_To_Point_Quadratic(const Field<T,2>& field,const Vector2i& cell,const Vector2& frac)
	{
		T val=Zero<T>();
		int i0,i1,j0,j1;
		if(frac[0]<(real).5){i0=-1;i1=1;}else{i0=0;i1=2;}
		if(frac[1]<(real).5){j0=-1;j1=1;}else{j0=0;j1=2;}
		for(int i=i0;i<=i1;i++){for(int j=j0;j<=j1;j++){
			if(!Grid<2>::Valid(cell+Vector2i(i,j),field.counts))continue;
			real w=Quadratic_Kernel(abs((real)i-frac[0]))*Quadratic_Kernel(abs((real)j-frac[1]));
			val+=w*field(cell+Vector2i(i,j));}}
		return val;
	}

	template<class T> T Nodes_To_Point_Quadratic(const Field<T,3>& field,const Vector3i& cell,const Vector3& frac)
	{
		T val=Zero<T>();
		int i0,i1,j0,j1,k0,k1;
		if(frac[0]<(real).5){i0=-1;i1=1;}else{i0=0;i1=2;}
		if(frac[1]<(real).5){j0=-1;j1=1;}else{j0=0;j1=2;}
		if(frac[2]<(real).5){k0=-1;k1=1;}else{k0=0;k1=2;}
		for(int i=i0;i<=i1;i++){for(int j=j0;j<=j1;j++){for(int k=k0;k<=k1;k++){
			if(!Grid<3>::Valid(cell+Vector3i(i,j,k),field.counts))continue;
			real w=Quadratic_Kernel(abs((real)i-frac[0]))*Quadratic_Kernel(abs((real)j-frac[1]))*Quadratic_Kernel(abs((real)k-frac[2]));
			val+=w*field(cell+Vector3i(i,j,k));}}}
		return val;
	}

	template<class T> void Point_To_Nodes_Quadratic(const T& value,Field<T,2>& field,Field<real,2>& weight,const Vector2i& cell,const Vector2& frac)
	{
		int i0,i1,j0,j1;
		if(frac[0]<(real).5){i0=-1;i1=1;}else{i0=0;i1=2;}
		if(frac[1]<(real).5){j0=-1;j1=1;}else{j0=0;j1=2;}
		for(int i=i0;i<=i1;i++){for(int j=j0;j<=j1;j++){
			if(!Grid<2>::Valid(cell+Vector2i(i,j),field.counts))continue;
			real w=Quadratic_Kernel(abs((real)i-frac[0]))*Quadratic_Kernel(abs((real)j-frac[1]));
			field(cell+Vector2i(i,j))+=w*value;
			weight(cell+Vector2i(i,j))+=w;}}
	}

	template<class T> void Point_To_Nodes_Quadratic(const T& value,Field<T,3>& field,Field<real,3>& weight,const Vector3i& cell,const Vector3& frac)
	{
		int i0,i1,j0,j1,k0,k1;
		if(frac[0]<(real).5){i0=-1;i1=1;}else{i0=0;i1=2;}
		if(frac[1]<(real).5){j0=-1;j1=1;}else{j0=0;j1=2;}
		if(frac[2]<(real).5){k0=-1;k1=1;}else{k0=0;k1=2;}
		for(int i=i0;i<=i1;i++){for(int j=j0;j<=j1;j++){for(int k=k0;k<=k1;k++){
			if(!Grid<3>::Valid(cell+Vector3i(i,j,k),field.counts))continue;
			real w=Quadratic_Kernel(abs((real)i-frac[0]))*Quadratic_Kernel(abs((real)j-frac[1]))*Quadratic_Kernel(abs((real)k-frac[2]));
			field(cell+Vector3i(i,j,k))+=w*value;
			weight(cell+Vector3i(i,j,k))+=w;}}}
	}

	//////////////////////////////////////////////////////////////////////////
	////Cubic interpolation
	inline real Cubic_Kernel(const real abs_x)
	{
		if(abs_x<1) return ((real).5*pow(abs_x,3)-pow(abs_x,2)+two_thirds);	
		else if(abs_x<2) return one_sixth*pow((real)2-abs_x,3);
		else return (real)0;
	}

	template<class T> T Nodes_To_Point_Cubic(const Field<T,2>& field,const Vector2i& cell,const Vector2& frac)
	{
		T val=Zero<T>();
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){
			if(!Grid<2>::Valid(cell+Vector2i(i,j),field.counts))continue;
			real w=Cubic_Kernel(abs((real)i-frac[0]))*Cubic_Kernel(abs((real)j-frac[1]));
			val+=w*field(cell+Vector2i(i,j));}}
		return val;
	}

	template<class T> T Nodes_To_Point_Cubic(const Field<T,3>& field,const Vector3i& cell,const Vector3& frac)
	{
		T val=Zero<T>();
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){for(int k=-1;k<=2;k++){
			if(!Grid<3>::Valid(cell+Vector3i(i,j,k),field.counts))continue;
			real w=Cubic_Kernel(abs((real)i-frac[0]))*Cubic_Kernel(abs((real)j-frac[1]))*Cubic_Kernel(abs((real)k-frac[2]));
			val+=w*field(cell+Vector3i(i,j,k));}}}
		return val;
	}

	template<class T> void Point_To_Nodes_Cubic(const T& value,Field<T,2>& field,Field<real,2>& weight,const Vector2i& cell,const Vector2& frac)
	{
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){
			if(!Grid<2>::Valid(cell+Vector2i(i,j),field.counts))continue;
			real w=Cubic_Kernel(abs((real)i-frac[0]))*Cubic_Kernel(abs((real)j-frac[1]));
			field(cell+Vector2i(i,j))+=w*value;
			weight(cell+Vector2i(i,j))+=w;}}
	}

	template<class T> void Point_To_Nodes_Cubic(const T& value,Field<T,3>& field,Field<real,3>& weight,const Vector3i& cell,const Vector3& frac)
	{
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){for(int k=-1;k<=2;k++){
			if(!Grid<3>::Valid(cell+Vector3i(i,j,k),field.counts))continue;
			real w=Cubic_Kernel(abs((real)i-frac[0]))*Cubic_Kernel(abs((real)j-frac[1]))*Cubic_Kernel(abs((real)k-frac[2]));
			field(cell+Vector3i(i,j,k))+=w*value;
			weight(cell+Vector3i(i,j,k))+=w;}}}
	}

	//////////////////////////////////////////////////////////////////////////
	////Catmull-Rom interpolation with clamp
	inline void Catmull_Rom_Coefficients_Helper(const real frac,real& a,real& b,real& c,real& d)
	{
		real f2=pow(frac,2);real f3=pow(frac,3);
		a=(real)-.5*frac+f2-(real).5*f3;
		b=(real)1-(real)2.5*f2+(real)1.5*f3;
		c=(real).5*frac+(real)2*f2-(real)1.5*f3;
		d=(real)-.5*f2+(real).5*f3;	
	}

	template<class T> T Nodes_To_Point_Catmull_Rom(const Field<T,1>& field,const Vector1i& cell,const Vector1& frac)
	{
		Vector4 cx;Catmull_Rom_Coefficients_Helper(frac[0],cx[0],cx[1],cx[2],cx[3]);int i=cell[0];
		T val=Zero<T>();
		T local_min=std::numeric_limits<T>::max();
		T local_max=-local_min;
		for(int i=-1;i<=2;i++){
			int nb=cell[0]+i;if(nb<0||nb>field.counts[0]-1)continue;
			val+=cx[i+1]*field(nb);
			if(field(nb)>local_max)local_max=field(nb);
			if(field(nb)<local_min)local_min=field(nb);}
		val=AuxFunc::Clamp(val,local_min,local_max);
		return val;
	}

	template<class T> T Nodes_To_Point_Catmull_Rom(const Field<T,2>& field,const Vector2i& cell,const Vector2& frac)
	{
		T val=Zero<T>();
		Vector4 c[2];for(int i=0;i<2;i++)Catmull_Rom_Coefficients_Helper(frac[i],c[i][0],c[i][1],c[i][2],c[i][3]);
		T local_min=std::numeric_limits<T>::max();
		T local_max=-local_min;
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){
			Vector2i nb=cell+Vector2i(i,j);if(!Grid<2>::Valid(nb,field.counts))continue;
			val+=c[0][i+1]*c[1][j+1]*field(nb);
			if(field(nb)>local_max)local_max=field(nb);
			if(field(nb)<local_min)local_min=field(nb);}}
		val=AuxFunc::Clamp(val,local_min,local_max);
		return val;
	}

	template<class T> T Nodes_To_Point_Catmull_Rom(const Field<T,3>& field,const Vector3i& cell,const Vector3& frac)
	{
		T val=Zero<T>();
		Vector4 c[3];for(int i=0;i<3;i++)Catmull_Rom_Coefficients_Helper(frac[i],c[i][0],c[i][1],c[i][2],c[i][3]);
		T local_min=std::numeric_limits<T>::max();
		T local_max=-local_min;
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){for(int k=-1;k<=2;k++){
			Vector3i nb=cell+Vector2i(i,j);if(!Grid<3>::Valid(nb,field.counts))continue;
			val+=c[0][i+1]*c[1][j+1]*c[2][k+1]*field(nb);
			if(field(nb)>local_max)local_max=field(nb);
			if(field(nb)<local_min)local_min=field(nb);}}}
		val=AuxFunc::Clamp(val,local_min,local_max);
		return val;
	}

	template<class T> void Point_To_Nodes_Catmull_Rom(const T& value,Field<T,2>& field,Field<real,2>& weight,const Vector2i& cell,const Vector2& frac)
	{
		Vector4 c[2];for(int i=0;i<2;i++)Catmull_Rom_Coefficients_Helper(frac[i],c[i][0],c[i][1],c[i][2],c[i][3]);
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){
			real w=c[0][i+1]*c[1][j+1];
			field(cell+Vector2i(i,j))+=w*value;
			weight(cell+Vector2i(i,j))+=w;}}
	}

	template<class T> void Point_To_Nodes_Catmull_Rom(const T& value,Field<T,3>& field,Field<real,3>& weight,const Vector3i& cell,const Vector3& frac)
	{
		Vector4 c[3];for(int i=0;i<3;i++)Catmull_Rom_Coefficients_Helper(frac[i],c[i][0],c[i][1],c[i][2],c[i][3]);
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){for(int k=-1;k<=2;k++){
			real w=c[0][i+1]*c[1][j+1]*c[2][k+1];
			field(cell+Vector3i(i,j,k))+=w*value;
			weight(cell+Vector3i(i,j,k))+=w;}}}
	}

	//////////////////////////////////////////////////////////////////////////
	////M4 interpolation
	inline real M4_Kernel(const real abs_x)
	{
		if(abs_x<(real)1) return (real)1-(real)2.5*pow(abs_x,2)+(real)1.5*pow(abs_x,3);
		else if(abs_x<(real)2) return (real).5*pow((real)2-abs_x,2)*((real)1-abs_x);
		else return (real)0;
	}

	template<class T> T Nodes_To_Point_M4(const Field<T,2>& field,const Vector2i& cell,const Vector2& frac)
	{
		T val=Zero<T>();
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){
			Vector2i c=cell+Vector2i(i,j);if(!Grid<2>::Valid(c,field.counts))continue;
			real w=M4_Kernel(abs((real)i-frac[0]))*M4_Kernel(abs((real)j-frac[1]));
			val+=w*field(c);}}
		return val;
	}

	template<class T> T Nodes_To_Point_M4(const Field<T,3>& field,const Vector3i& cell,const Vector3& frac)
	{
		T val=Zero<T>();
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){for(int k=-1;k<=2;k++){
			Vector3i c=cell+Vector3i(i,j,k);if(!Grid<3>::Valid(c,field.counts))continue;
			real w=M4_Kernel(abs((real)i-frac[0]))*M4_Kernel(abs((real)j-frac[1]))*M4_Kernel(abs((real)k-frac[2]));
			val+=w*field(c);}}}
		return val;
	}

	template<class T> void Point_To_Nodes_M4(const T& value,Field<T,2>& field,Field<real,2>& weight,const Vector2i& cell,const Vector2& frac)
	{
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){
			if(!Grid<2>::Valid(cell+Vector2i(i,j),field.counts))continue;
			real w=M4_Kernel(abs((real)i-frac[0]))*M4_Kernel(abs((real)j-frac[1]));
			field(cell+Vector2i(i,j))+=w*value;
			weight(cell+Vector2i(i,j))+=w;}}
	}

	template<class T> void Point_To_Nodes_M4(const T& value,Field<T,3>& field,Field<real,3>& weight,const Vector3i& cell,const Vector3& frac)
	{
		for(int i=-1;i<=2;i++){for(int j=-1;j<=2;j++){for(int k=-1;k<=2;k++){
			if(!Grid<3>::Valid(cell+Vector3i(i,j,k),field.counts))continue;
			real w=M4_Kernel(abs((real)i-frac[0]))*M4_Kernel(abs((real)j-frac[1]))*M4_Kernel(abs((real)k-frac[2]));
			field(cell+Vector3i(i,j,k))+=w*value;
			weight(cell+Vector3i(i,j,k))+=w;}}}
	}

	//////////////////////////////////////////////////////////////////////////
	////Quaternion + linear
	////assuming T is Quaternion
	////implementation refers to http://wiki.unity3d.com/index.php/Averaging_Quaternions_and_Vectors?_ga=2.9199221.1942645287.1577671110-718620556.1577671110
	inline Vector4 Weighted_Averaged_Quaternions(const Array<Vector4>& q,const Array<int>& idx,const Array<real>& w)
	{
		Vector4 cum_q=q[idx[0]]*w[0];
		for(int i=1;i<idx.size();i++){
			if(cum_q.dot(q[idx[i]])<(real)0){cum_q-=q[idx[i]]*w[i];}
			else cum_q+=q[idx[i]]*w[i];}
		return cum_q.normalized();
	}

	inline Vector4 Weighted_Averaged_Quaternions(const Array<Vector4>& q,const Array<real>& w)
	{
		Vector4 cum_q=q[0]*w[0];
		for(int i=1;i<w.size();i++){
			if(cum_q.dot(q[i])<(real)0){cum_q-=q[i]*w[i];}
			else cum_q+=q[i]*w[i];}
		return cum_q.normalized();
	}

	template<class T> T Nodes_To_Point_Quaternion(const Field<T,2>& field,const Vector2i& cell,const Vector2& frac)
	{
		const T& q00=field(cell);
		const T& q10=field(Vector2i(cell[0]+1,cell[1]));
		const T& q01=field(Vector2i(cell[0],cell[1]+1));
		const T& q11=field(Vector2i(cell[0]+1,cell[1]+1));
		real w00=((real)1-frac[0])*((real)1-frac[1]);
		real w10=frac[0]*((real)1-frac[1]);
		real w01=((real)1-frac[0])*frac[1];
		real w11=frac[0]*frac[1];

		Vector4 cum_q=Vector4::Zero();
		#define Accumulate_Q(q,w) \
		if(cum_q.dot(q)<(real)0){cum_q-=q*w;}else cum_q+=q*w;
		Accumulate_Q(q00,w00);
		Accumulate_Q(q10,w10);
		Accumulate_Q(q01,w01);
		Accumulate_Q(q11,w11);
		#undef Accumulate_Q

		return cum_q.normalized();
	}

	template<class T> T Nodes_To_Point_Quaternion(const Field<T,3>& field,const Vector3i& cell,const Vector3& frac)
	{	
		real w000=((real)1-frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		real w100=(frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		real w010=((real)1-frac[0])*(frac[1])*((real)1-frac[2]);
		real w001=((real)1-frac[0])*((real)1-frac[1])*(frac[2]);
		real w110=(frac[0])*(frac[1])*((real)1-frac[2]);
		real w011=((real)1-frac[0])*(frac[1])*(frac[2]);
		real w101=(frac[0])*((real)1-frac[1])*(frac[2]);
		real w111=(frac[0])*(frac[1])*(frac[2]);

		Vector4 cum_q=Vector4::Zero();
		#define Accumulate_Q(q,w) \
		cum_q+=q*w;
		//if(cum_q.dot(q)<(real)0){cum_q-=q*w;std::cout<<"n2p ";}else cum_q+=q*w;
		Accumulate_Q(field(cell[0],cell[1],cell[2]),w000);
		Accumulate_Q(field(cell[0]+1,cell[1],cell[2]),w100);
		Accumulate_Q(field(cell[0],cell[1]+1,cell[2]),w010);
		Accumulate_Q(field(cell[0],cell[1],cell[2]+1),w001);
		Accumulate_Q(field(cell[0]+1,cell[1]+1,cell[2]),w110);
		Accumulate_Q(field(cell[0],cell[1]+1,cell[2]+1),w011);
		Accumulate_Q(field(cell[0]+1,cell[1],cell[2]+1),w101);
		Accumulate_Q(field(cell[0]+1,cell[1]+1,cell[2]+1),w111);
		#undef Accumulate_Q

		return cum_q;
	}

	template<class T> void Point_To_Nodes_Quaternion(const T& value,Field<T,2>& field,Field<real,2>& weight,const Vector2i& cell,const Vector2& frac)
	{
		real w00=((real)1-frac[0])*((real)1-frac[1]);
		real w10=frac[0]*((real)1-frac[1]);
		real w01=((real)1-frac[0])*frac[1];
		real w11=frac[0]*frac[1];

		weight(cell[0],cell[1])+=w00;
		weight(cell[0]+1,cell[1])+=w10;
		weight(cell[0],cell[1]+1)+=w01;
		weight(cell[0]+1,cell[1]+1)+=w11;

		#define Accumulate_Q(cum_q,q,w) \
		if(cum_q.dot(q)<(real)0){cum_q-=q*w;}else cum_q+=q*w;
		Accumulate_Q(field(cell[0],cell[1]),value,w00);
		Accumulate_Q(field(cell[0]+1,cell[1]),value,w10);
		Accumulate_Q(field(cell[0],cell[1]+1),value,w01);
		Accumulate_Q(field(cell[0]+1,cell[1]+1),value,w11);
		#undef Accumulate_Q

	}

	template<class T> void Point_To_Nodes_Quaternion(const T& value,Field<T,3>& field,Field<real,3>& weight,const Vector3i& cell,const Vector3& frac)
	{
		real w000=((real)1-frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		real w100=(frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		real w010=((real)1-frac[0])*(frac[1])*((real)1-frac[2]);
		real w001=((real)1-frac[0])*((real)1-frac[1])*(frac[2]);
		real w110=(frac[0])*(frac[1])*((real)1-frac[2]);
		real w011=((real)1-frac[0])*(frac[1])*(frac[2]);
		real w101=(frac[0])*((real)1-frac[1])*(frac[2]);
		real w111=(frac[0])*(frac[1])*(frac[2]);

		weight(cell[0],cell[1],cell[2])+=w000;
		weight(cell[0]+1,cell[1],cell[2])+=w100;
		weight(cell[0],cell[1]+1,cell[2])+=w010;
		weight(cell[0],cell[1],cell[2]+1)+=w001;
		weight(cell[0]+1,cell[1]+1,cell[2])+=w110;
		weight(cell[0],cell[1]+1,cell[2]+1)+=w011;
		weight(cell[0]+1,cell[1],cell[2]+1)+=w101;
		weight(cell[0]+1,cell[1]+1,cell[2]+1)+=w111;

		#define Accumulate_Q(cum_q,q,w) \
		cum_q+=q*w;
		//if(cum_q.dot(q)<(real)0){cum_q-=q*w;std::cout<<"p2n ";}else cum_q+=q*w;
		Accumulate_Q(field(cell[0],cell[1],cell[2]),value,w000);
		Accumulate_Q(field(cell[0]+1,cell[1],cell[2]),value,w100);
		Accumulate_Q(field(cell[0],cell[1]+1,cell[2]),value,w010);
		Accumulate_Q(field(cell[0],cell[1],cell[2]+1),value,w001);
		Accumulate_Q(field(cell[0]+1,cell[1]+1,cell[2]),value,w110);
		Accumulate_Q(field(cell[0],cell[1]+1,cell[2]+1),value,w011);
		Accumulate_Q(field(cell[0]+1,cell[1],cell[2]+1),value,w101);
		Accumulate_Q(field(cell[0]+1,cell[1]+1,cell[2]+1),value,w111);
		#undef Accumulate_Q
	}

	//////////////////////////////////////////////////////////////////////////
	////Linear interpolation with phi function
	template<class T> T Nodes_To_Point_Linear_With_Phi(const Field<T,1>& field,const Vector1i& cell,const Grid<1>& grid,
		std::function<bool(const Vector1&)> Phi,const Vector1& frac)
	{
		real phi0=Phi(grid.Node(cell));
		real phi1=Phi(grid.Node(Vec1i(cell[0]+1)));
		real w0=phi0<(real)0?((real)1-frac[0]):(real)0;
		real w1=phi1<(real)0?frac[0]:(real)0;
		real total_w=w0+w1;
		if(total_w==(real)0)return (real)0;
		return (w0*field(cell)+w1*field(Vec1i(cell[0]+1)))/total_w;
	}

	template<class T> T Nodes_To_Point_Linear_With_Phi(const Field<T,2>& field,const Vector2i& cell,const Grid<2>& grid,
		std::function<real(const Vector2&)> Phi,const Vector2& frac)
	{
		real w[4],phi[4];
		w[0]=((real)1-frac[0])*((real)1-frac[1]);
		w[1]=frac[0]*((real)1-frac[1]);
		w[2]=((real)1-frac[0])*frac[1];
		w[3]=frac[0]*frac[1];

		phi[0]=Phi(grid.Node(cell));
		phi[1]=Phi(grid.Node(Vector2i(cell[0]+1,cell[1])));
		phi[2]=Phi(grid.Node(Vector2i(cell[0],cell[1]+1)));
		phi[3]=Phi(grid.Node(Vector2i(cell[0]+1,cell[1]+1)));

		for(int i=0;i<4;i++)w[i]=(phi[i]>=(real)0)?(real)0:w[i];
		real total_w=w[0]+w[1]+w[2]+w[3];

		////if all nodes are positive, return the minimal one
		if(total_w==(real)0){
			int min_idx=0;real min_phi=phi[0];
			for(int i=1;i<4;i++)if(phi[i]<min_phi){min_idx=i;min_phi=phi[i];}
			
			w[min_idx]=(real)1;total_w=(real)1;

			/////*dbg output*/
			//real result=(w[0]*field(cell)+w[1]*field(cell[0]+1,cell[1])+w[2]*field(cell[0],cell[1]+1)+w[3]*field(cell[0]+1,cell[1]+1))/total_w;
			//std::cout<<"total_w=0: cell: "<<cell[0]<<", "<<cell[1]<<", phi: "<<phi[0]<<", "<<phi[1]<<", "<<phi[2]<<", "<<phi[3]<<", frac: "<<frac[0]<<", "<<frac[1]<<std::endl;
			//std::cout<<"total_w==0: "<<min_idx<<", "<<min_phi<<", result: "<<result<<std::endl;
		}

		return (w[0]*field(cell)
			+w[1]*field(cell[0]+1,cell[1])
			+w[2]*field(cell[0],cell[1]+1)
			+w[3]*field(cell[0]+1,cell[1]+1))/total_w;
	}

	template<class T> T Nodes_To_Point_Linear_With_Phi(const Field<T,3>& field,const Vector3i& cell,const Grid<3>& grid,
		std::function<real(const Vector3&)> Phi,const Vector3& frac)
	{
		real w[8],phi[8];
		w[0]=((real)1-frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		w[1]=frac[0]*((real)1-frac[1])*((real)1-frac[2]);
		w[2]=((real)1-frac[0])*frac[1]*((real)1-frac[2]);
		w[3]=((real)1-frac[0])*((real)1-frac[1])*frac[2];
		w[4]=frac[0]*frac[1]*((real)1-frac[2]);
		w[5]=frac[0]*((real)1-frac[1])*frac[2];
		w[6]=((real)1-frac[0])*frac[1]*frac[2];
		w[7]=frac[0]*frac[1]*frac[2];

		phi[0]=Phi(grid.Node(cell));
		phi[1]=Phi(grid.Node(Vector3i(cell[0]+1,cell[1],cell[2])));
		phi[2]=Phi(grid.Node(Vector3i(cell[0],cell[1]+1,cell[2])));
		phi[3]=Phi(grid.Node(Vector3i(cell[0],cell[1],cell[2]+1)));
		phi[4]=Phi(grid.Node(Vector3i(cell[0]+1,cell[1]+1,cell[2])));
		phi[5]=Phi(grid.Node(Vector3i(cell[0]+1,cell[1],cell[2]+1)));
		phi[6]=Phi(grid.Node(Vector3i(cell[0],cell[1]+1,cell[2]+1)));
		phi[7]=Phi(grid.Node(Vector3i(cell[0]+1,cell[1]+1,cell[2]+1)));

		real total_w=(real)0;
		for(int i=0;i<8;i++){w[i]=(phi[i]>=(real)0)?(real)0:w[i]; total_w+=w[i];}
		
		////if all nodes are positive, return the minimal one
		if(total_w==(real)0){
			int min_idx=0;real min_phi=phi[0];
			for(int i=1;i<8;i++)if(phi[i]<min_phi){min_idx=i;min_phi=phi[i];}

			w[min_idx]=(real)1;total_w=(real)1;}

		return (w[0]*field(cell[0],cell[1],cell[2])
			+w[1]*field(cell[0]+1,cell[1],cell[2])
			+w[2]*field(cell[0],cell[1]+1,cell[2])
			+w[3]*field(cell[0],cell[1],cell[2]+1)
			+w[4]*field(cell[0]+1,cell[1]+1,cell[2])
			+w[5]*field(cell[0]+1,cell[1],cell[2]+1)
			+w[6]*field(cell[0],cell[1]+1,cell[2]+1)
			+w[7]*field(cell[0]+1,cell[1]+1,cell[2]+1))/total_w;
	}

	//////////////////////////////////////////////////////////////////////////
	////Linear interpolation with ghost value function
	template<class T> T Nodes_To_Point_Linear_With_Ghost(const Field<T,1>& field,const Vector1i& cell,const Grid<1>& grid,const Vector1& frac,
		std::function<T(const Vector1&,const int)> Ghost_Component,const int axis=0)
	{
		real w0=(real)1-frac[0];real w1=frac[0];
		Vector1i n0=cell, n1(cell[0]+1);
	
		T val0=grid.Valid_Node(n0)?field(n0):Ghost_Component(grid.Node(n0),axis);
		T val1=grid.Valid_Node(n1)?field(n1):Ghost_Component(grid.Node(n1),axis);
		return w0*val0+w1*val1;
	}

	template<class T> T Nodes_To_Point_Linear_With_Ghost(const Field<T,2>& field,const Vector2i& cell,const Grid<2>& grid,const Vector2& frac,
		std::function<T(const Vector2&,const int)> Ghost_Component,const int axis=0)
	{
		real w[4];
		w[0]=((real)1-frac[0])*((real)1-frac[1]);
		w[1]=frac[0]*((real)1-frac[1]);
		w[2]=((real)1-frac[0])*frac[1];
		w[3]=frac[0]*frac[1];

		Vector2i n0=cell;
		Vector2i n1(cell[0]+1,cell[1]);
		Vector2i n2(cell[0],cell[1]+1);
		Vector2i n3(cell[0]+1,cell[1]+1);

		T val0=grid.Valid_Node(n0)?field(n0):Ghost_Component(grid.Node(n0),axis);
		T val1=grid.Valid_Node(n1)?field(n1):Ghost_Component(grid.Node(n1),axis);
		T val2=grid.Valid_Node(n2)?field(n2):Ghost_Component(grid.Node(n2),axis);
		T val3=grid.Valid_Node(n3)?field(n3):Ghost_Component(grid.Node(n3),axis);

		//if(!grid.Valid_Node(n0)||!grid.Valid_Node(n1)||!grid.Valid_Node(n2)||!grid.Valid_Node(n3)){
		//	//std::cout<<grid.cell_counts.transpose()<<", "<<grid.node_counts.transpose()<<": "<<cell.transpose()<<std::endl;
		//	return (T)0;}

		return w[0]*val0+w[1]*val1+w[2]*val2+w[3]*val3;
	}

	template<class T> T Nodes_To_Point_Linear_With_Ghost(const Field<T,3>& field,const Vector3i& cell,const Grid<3>& grid,const Vector3& frac,
		std::function<T(const Vector3&,const int)> Ghost_Component,const int axis=0)
	{
		real w[8];
		w[0]=((real)1-frac[0])*((real)1-frac[1])*((real)1-frac[2]);
		w[1]=frac[0]*((real)1-frac[1])*((real)1-frac[2]);
		w[2]=((real)1-frac[0])*frac[1]*((real)1-frac[2]);
		w[3]=((real)1-frac[0])*((real)1-frac[1])*frac[2];
		w[4]=frac[0]*frac[1]*((real)1-frac[2]);
		w[5]=frac[0]*((real)1-frac[1])*frac[2];
		w[6]=((real)1-frac[0])*frac[1]*frac[2];
		w[7]=frac[0]*frac[1]*frac[2];

		Vector3i n0=cell;
		Vector3i n1(cell[0]+1,cell[1],cell[2]);
		Vector3i n2(cell[0],cell[1]+1,cell[2]);
		Vector3i n3(cell[0],cell[1],cell[2]+1);
		Vector3i n4(cell[0]+1,cell[1]+1,cell[2]);
		Vector3i n5(cell[0]+1,cell[1],cell[2]+1);
		Vector3i n6(cell[0],cell[1]+1,cell[2]+1);
		Vector3i n7(cell[0]+1,cell[1]+1,cell[2]+1);

		T val0=grid.Valid_Node(n0)?field(n0):Ghost_Component(grid.Node(n0),axis);
		T val1=grid.Valid_Node(n1)?field(n1):Ghost_Component(grid.Node(n1),axis);
		T val2=grid.Valid_Node(n2)?field(n2):Ghost_Component(grid.Node(n2),axis);
		T val3=grid.Valid_Node(n3)?field(n3):Ghost_Component(grid.Node(n3),axis);
		T val4=grid.Valid_Node(n4)?field(n4):Ghost_Component(grid.Node(n4),axis);
		T val5=grid.Valid_Node(n5)?field(n5):Ghost_Component(grid.Node(n5),axis);
		T val6=grid.Valid_Node(n6)?field(n6):Ghost_Component(grid.Node(n6),axis);
		T val7=grid.Valid_Node(n7)?field(n7):Ghost_Component(grid.Node(n7),axis);

		return w[0]*val0+w[1]*val1+w[2]*val2+w[3]*val3+w[4]*val4+w[5]*val5+w[6]*val6+w[7]*val7;
	}
};

#define Declare_Intp_Func_Grid_To_Point(FuncName,IntpName) \
template<class T> T Interpolate_Nodes##FuncName(const Field<T,d>& field,const VectorD& pos) const	\
{VectorD frac;VectorDi cell=mac_grid.grid.Clamped_Cell_Coord_With_Fraction(pos,frac);return IntpFunc::Nodes_To_Point##IntpName(field,cell,frac);} \
template<class T> T Interpolate_Centers##FuncName(const Field<T,d>& field,const VectorD& pos) const \
 {VectorD frac;VectorDi cell=center_grid.Clamped_Cell_Coord_With_Fraction(pos,frac);return IntpFunc::Nodes_To_Point##IntpName(field,cell,frac);} \
template<class T> T Interpolate_Faces##FuncName(const FaceField<T,d>& field,const VectorD& pos,const int axis) const \
{VectorD frac;VectorDi cell=mac_grid.face_grids[axis].Clamped_Cell_Coord_With_Fraction(pos,frac);return IntpFunc::Nodes_To_Point##IntpName(field.face_fields[axis],cell,frac);} \
VectorD Interpolate_Face_Vectors##FuncName(const FaceField<real,d>& field,const VectorD& pos) const \
{VectorD v;for(int i=0;i<d;i++)v[i]=Interpolate_Faces##FuncName(field,pos,i);return v;}

#define Declare_Intp_Func_Point_To_Grid(FuncName,IntpName) \
template<class T> void Interpolate_Point_To_Cells##FuncName(const VectorD& pos,const T& value,Field<T,d>& field,Field<real,d>& weights) const \
{VectorD frac;VectorDi cell=center_grid.Clamped_Cell_Coord_With_Fraction(pos,frac);IntpFunc::Point_To_Nodes##IntpName(value,field,weights,cell,frac);} \
template<class T> void Interpolate_Point_To_Nodes##FuncName(const VectorD& pos,const T& value,Field<T,d>& field,Field<real,d>& weights) const \
{VectorD frac;VectorDi cell=mac_grid.grid.Clamped_Cell_Coord_With_Fraction(pos,frac);IntpFunc::Point_To_Nodes##IntpName(value,field,weights,cell,frac);} \
template<class T> void Interpolate_Point_To_Faces##FuncName(const VectorD& pos,const Vector<T,d>& value,FaceField<T,d>& face_field,FaceField<real,d>& face_weights) const \
{for(int axis=0;axis<d;axis++){VectorD frac;VectorDi cell=mac_grid.face_grids[axis].Clamped_Cell_Coord_With_Fraction(pos,frac);IntpFunc::Point_To_Nodes##IntpName(value[axis],face_field.face_fields[axis],face_weights.face_fields[axis],cell,frac);}}

template<int d> class Interpolation
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	MacGrid<d> mac_grid;
	Grid<d> center_grid;

	Interpolation(const MacGrid<d>& _mac_grid):mac_grid(_mac_grid){center_grid=mac_grid.grid.Cell_Grid_To_Node_Grid();}
	Interpolation(const Grid<d>& _grid):mac_grid(_grid){center_grid=mac_grid.grid.Cell_Grid_To_Node_Grid();}

	//////////////////////////////////////////////////////////////////////////
	////interpolate from a grid to a point
	Declare_Intp_Func_Grid_To_Point(,_Linear);
	Declare_Intp_Func_Grid_To_Point(_Linear,_Linear);
	Declare_Intp_Func_Grid_To_Point(_Quadratic,_Quadratic);
	Declare_Intp_Func_Grid_To_Point(_Cubic,_Cubic);
	Declare_Intp_Func_Grid_To_Point(_Catmull_Rom,_Catmull_Rom);
	Declare_Intp_Func_Grid_To_Point(_M4,_M4);
	Declare_Intp_Func_Grid_To_Point(_Quaternion,_Quaternion);

	//////////////////////////////////////////////////////////////////////////
	////interpolate from a point to a grid
	Declare_Intp_Func_Point_To_Grid(,_Linear);
	Declare_Intp_Func_Point_To_Grid(_Linear,_Linear);
	Declare_Intp_Func_Point_To_Grid(_Quadratic,_Quadratic);
	Declare_Intp_Func_Point_To_Grid(_Cubic,_Cubic);
	Declare_Intp_Func_Point_To_Grid(_Catmull_Rom,_Catmull_Rom);
	Declare_Intp_Func_Point_To_Grid(_M4,_M4);
	Declare_Intp_Func_Point_To_Grid(_Quaternion,_Quaternion);

	//////////////////////////////////////////////////////////////////////////
	////interpolate from a grid from to another grid form
	template<class T> T Interpolate_Faces_To_Cell(const FaceField<T,d>& face_field,const VectorDi& cell,const int axis) const
	{return (real).5*(face_field.face_fields[axis](cell)+face_field.face_fields[axis](cell+VectorDi::Unit(axis)));}

	//////////////////////////////////////////////////////////////////////////
	////interpolate between two grids
	////The following functions can all be accelerated by performing better grid indexing and memory access and grid boundary treatment

	template<class T> void Interpolate_Faces_To_Nodes(const FaceField<T,d>& face_field,Field<Vector<T,d>,d>& node_field) const
	{
		int node_num=mac_grid.grid.Number_Of_Nodes();
		#pragma omp parallel for
		for(int i=0;i<node_num;i++){const VectorDi& node=mac_grid.grid.Node_Coord(i);		
			Vector<T,d> v=Vector<T,d>::Zero();
				for(int axis=0;axis<d;axis++){int n=0;for(int i=0;i<MacGrid<d>::Number_Of_Node_Incident_Faces_Per_Axis();i++){
					const VectorDi face=MacGrid<d>::Node_Incident_Face_Per_Axis(node,i,axis);
					if(mac_grid.Valid_Face(axis,face)){v[axis]+=face_field(axis,face);n++;}}if(n>0)v[axis]/=(T)n;}node_field(node)=v;}
	}

	template<class T> void Interpolate_Faces_To_Cells(const FaceField<T,d>& face_field,Field<Vector<T,d>,d>& cell_field) const
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){const VectorDi& cell=mac_grid.grid.Cell_Coord(i);Vector<T,d>v;
			for(int axis=0;axis<d;axis++){v[axis]=(T)((real).5*(face_field(axis,cell)+face_field(axis,cell+Vector<int,d>::Unit(axis))));}cell_field(cell)=v;}
	}

	template<class T> void Interpolate_Cells_To_Faces(const Field<Vector<T,d>,d>& cell_field,FaceField<T,d>& face_field) const
	{
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.face_grids[axis].Number_Of_Nodes();
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				T value=Zero<T>();int n=0;
				for(int i=0;i<2;i++){
					VectorDi cell=mac_grid.Face_Incident_Cell(axis,face,i);
					if(mac_grid.grid.Valid_Cell(cell)){value+=cell_field(cell)[axis];n++;}}
				if(n>0)value=(T)(value/(real)n);face_field(axis,face)=value;}}
	}

	template<class T> void Interpolate_Cells_To_Faces_Fast(const Field<Vector<T,d>,d>& cell_field,FaceField<T,d>& face_field) const
	{
		const VectorDi& cell_counts=mac_grid.grid.cell_counts;

		for(int axis=0;axis<d;axis++){
			const int a1=(axis+1)%d;const int a2=(axis+2)%d;
			const VectorDi& face_counts=mac_grid.face_counts[axis];
			Vector3i s,e;
			s[axis]=1;
			s[a1]=s[a2]=0;
			e[axis]=face_counts[axis]-1;
			e[a1]=face_counts[a1];
			e[a2]=face_counts[a2];

			#pragma omp parallel for
			for(int i=s[0];i<e[0];i++){
				for(int j=s[1];j<e[1];j++){
					for(int k=s[2];k<e[2];k++){
						Vector3i v(i,j,k);v[axis]--;
						face_field(axis,i,j,k)=(real)(.5)*(cell_field(i,j,k)[axis]+cell_field(v[0],v[1],v[2])[axis]);}}}}

	}

	template<class T> void Interpolate_Cells_To_Nodes(const Field<T,d>& cell_field,Field<T,d>& node_field) const
	{
		int node_num=mac_grid.grid.Number_Of_Nodes();
		#pragma omp parallel for
		for(int i=0;i<node_num;i++){const VectorDi& node=mac_grid.grid.Node_Coord(i);
			T value=Zero<T>();int n=0;
			for(int i=0;i<Grid<d>::Number_Of_Node_Incident_Cells();i++){
				VectorDi cell=mac_grid.grid.Node_Incident_Cell(node,i);
				if(mac_grid.grid.Valid_Cell(cell)){value+=cell_field(cell);n++;}}
			if(n>0)value=(T)(value/(real)n);node_field(node)=value;}
	}

	template<class T> void Interpolate_Nodes_To_Cells(const Field<T,d>& node_field,Field<T,d>& cell_field) const
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			T value=Zero<T>();int n=0;
			for(int i=0;i<Grid<d>::Number_Of_Cell_Incident_Nodes();i++){
				VectorDi node=mac_grid.grid.Cell_Incident_Node(cell,i);if(mac_grid.grid.Valid_Node(node)){value+=node_field(node);n++;}}
			if(n>0)value/=(real)n;cell_field(cell)=value;}
	}

	//////////////////////////////////////////////////////////////////////////
	////interpolate gradient from a grid to a point
	VectorD Gradient_Faces(const FaceField<real,d>& field,const VectorD& pos,const int axis) const
	{
		VectorD g;for(int i=0;i<d;i++){
			VectorD pos_neg=pos-VectorD::Unit(i)*mac_grid.grid.dx;real value_neg=Interpolate_Faces(field,pos_neg,axis);
			VectorD pos_pos=pos+VectorD::Unit(i)*mac_grid.grid.dx;real value_pos=Interpolate_Faces(field,pos_pos,axis);
			g[i]=(value_pos-value_neg)/((real)2*mac_grid.grid.dx);}return g;	
	}

	MatrixD Gradient_Face_Vectors(const FaceField<real,d>& field,const VectorD& pos) const
	{MatrixD g;for(int i=0;i<d;i++){VectorD g_axis=Gradient_Faces(field,pos,i);g.row(i)=g_axis.transpose();}return g;}

	//////////////////////////////////////////////////////////////////////////
	////interpolate with phi function
	VectorD Interpolate_Face_Vectors_With_Phi(const FaceField<real,d>& field,std::function<real(const VectorD&)> Phi,const VectorD& pos) const
	{
		VectorD v;for(int i=0;i<d;i++)v[i]=Interpolate_Faces_With_Phi(field,Phi,pos,i);return v;
	}

	template<class T> T Interpolate_Faces_With_Phi(const FaceField<T,d>& field,std::function<real(const VectorD&)> Phi,const VectorD& pos,const int axis) const \
	{
		VectorD frac;VectorDi cell=mac_grid.face_grids[axis].Clamped_Cell_Coord_With_Fraction(pos,frac);
		return IntpFunc::Nodes_To_Point_Linear_With_Phi(field.face_fields[axis],cell,mac_grid.face_grids[axis],Phi,frac);
	}

	template<class T> T Interpolate_Faces_With_Phi_Linear(const FaceField<T,d>& field,std::function<real(const VectorD&)> Phi,const VectorD& pos,const int axis) const \
	{
		VectorD frac;VectorDi cell=mac_grid.face_grids[axis].Clamped_Cell_Coord_With_Fraction(pos,frac);
		return IntpFunc::Nodes_To_Point_Linear_With_Phi(field.face_fields[axis],cell,mac_grid.face_grids[axis],Phi,frac);
	}

	template<class T> T Interpolate_Centers_With_Phi(const Field<T,d>& field,std::function<real(const VectorD&)> Phi,const VectorD& pos) const \
	{
		VectorD frac;VectorDi cell=center_grid.Clamped_Cell_Coord_With_Fraction(pos,frac);
		return IntpFunc::Nodes_To_Point_Linear_With_Phi(field,cell,mac_grid.grid,Phi,frac);
	}

	//////////////////////////////////////////////////////////////////////////
	////interpolate with ghost values
	VectorD Interpolate_Face_Vectors_With_Ghost(const FaceField<real,d>& field,const VectorD& pos,std::function<real(const VectorD&,const int)> Ghost_Component) const
	{
		VectorD v;for(int axis=0;axis<d;axis++){
			VectorD frac;VectorDi cell=mac_grid.face_grids[axis].Cell_Coord_With_Fraction_No_Clamp(pos,frac);	////no clamp
			v[axis]=IntpFunc::Nodes_To_Point_Linear_With_Ghost<real>(field.face_fields[axis],cell,mac_grid.face_grids[axis],frac,Ghost_Component,axis);}return v;
	}

	template<class T> T Interpolate_Faces_With_Ghost_Linear(const FaceField<T,d>& field,const VectorD& pos,std::function<T(const VectorD&,const int)> Ghost_Component,const int axis) const \
	{
		VectorD frac;VectorDi cell=mac_grid.face_grids[axis].Cell_Coord_With_Fraction_No_Clamp(pos,frac);	////no clamp
		return IntpFunc::Nodes_To_Point_Linear_With_Ghost<T>(field.face_fields[axis],cell,mac_grid.face_grids[axis],frac,Ghost_Component,axis);
	}

	template<class T> T Interpolate_Centers_With_Ghost_Linear(const Field<T,d>& field,const VectorD& pos,std::function<T(const VectorD&,const int)> Ghost_Component) const \
	{
		VectorD frac;VectorDi cell=center_grid.Cell_Coord_With_Fraction_No_Clamp(pos,frac);	////no clamp
		return IntpFunc::Nodes_To_Point_Linear_With_Ghost<T>(field,cell,center_grid,frac,Ghost_Component);
	}
protected:
};
#endif

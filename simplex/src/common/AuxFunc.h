//////////////////////////////////////////////////////////////////////////
// Common auxiliary functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __AuxFunc_h__
#define __AuxFunc_h__
#include <algorithm>
#include "Common.h"

namespace AuxFunc{
    ////Aux functions

	////Scalar Math
	real Quick_Pow(real a, int n);//must ensure n>=0

    ////Comparison
    template<class ArrayT> bool All_Less(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]>=a1[i])return false;}return true;}
    template<class ArrayT> bool All_Less_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]>a1[i])return false;}return true;}
    template<class ArrayT> bool All_Greater(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]<=a1[i])return false;}return true;}
    template<class ArrayT> bool All_Greater_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++){if(a0[i]<a1[i])return false;}return true;}
    template<class ArrayT> bool Has_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++)if(a0[i]==a1[i])return true;return false;}
    template<class ArrayT> bool Has_Less_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++)if(a0[i]<=a1[i])return true;return false;}
    template<class ArrayT> bool Has_Greater_Equal(const ArrayT& a0,const ArrayT& a1){for(auto i=0;i<a0.size();i++)if(a0[i]>=a1[i])return true;return false;}
	template<int d> bool Outside(const Vector<int,d>& domain,const Vector<int,d>& p);
	
    ////Clamp
    template<class T> T inline Clamp(const T x,const T x_min,const T x_max)
    {return std::max(x_min,std::min(x,x_max));}
    template<class T,int d> inline Vector<T,d> Clamp_Vector(const Vector<T,d>& x,const Vector<T,d>& x_min,const Vector<T,d>& x_max)
    {Vector<T,d> clamped;for(int i=0;i<d;i++)clamped[i]=Clamp(x[i],x_min[i],x_max[i]);return clamped;}

    ////Vector math
    Vector1 Cross(const Vector2& v1,const Vector2& v2);
	Vector2 Cross(const Vector1& v1,const Vector2& v2);
	Vector2 Cross(const real& v1,const Vector2& v2);
	Matrix3 Skew(const Vector3& v);
	Vector1 Orthogonal_Vector(const Vector1& v);
    Vector2 Orthogonal_Vector(const Vector2& v);	////this is the orthogonal vector on the *left* hand
	Vector3 Orthogonal_Vector(const Vector3& v);
	Vector4 Orthogonal_Vector(const Vector4& v);
	template<class T, int d> Vector<T, d> Component_Along(const Vector<T, d>& vec, const Vector<T, d>& dir) { return dir * vec.dot(dir) / dir.squaredNorm(); }
	template<class T,int d> Vector<T,d> Eliminate_Unit_Component(const Vector<T,d>& vec, const Vector<T,d>& unit_norm) { return vec - vec.dot(unit_norm) * unit_norm; }//length of norm MUST be 1

    template<int d> inline real Norm(const Vector<real,d>& v);
    template<class T,int d> inline Vector<T,d> Unit_Orthogonal_Vector(const Vector<T,d>& v);

	////Vector angles, assuming v1 and v2 are unit vectors
	real Angle_Between(const Vector2& v1,const Vector2& v2);	////[0,pi]
	real Angle_Between(const Vector3& v1,const Vector3& v2);	////[0,pi]
	real Angle_From_To(const Vector2& v1,const Vector2& v2);	////[-pi,pi], default is +z
	real Angle_From_To(const Vector3& v1,const Vector3& v2);	////[-pi,pi], no axis specified, [0,pi]
	real Angle_From_To_With_Specified_Axis(const Vector3& v1,const Vector3& v2,const Vector3& axis);
	void Angle_And_Axis_Between(const Vector3d& v1,const Vector3d& v2,double& theta,Vector3d& axis);		////theta [0,pi]
	void Angle_And_Axis_Between(const Vector3f& v1,const Vector3f& v2,float& theta,Vector3f& axis);	////float version
	real Dihedral_Angle(const Vector2& a,const Vector2& c,const Vector2& b,const Vector2& d);	////{0,pi}
	real Dihedral_Angle(const Vector3& a,const Vector3& c,const Vector3& b,const Vector3& d);	////dihedral angle between abc and acd, [0,pi]
	
	////Distances and closest points
	template<int d> Vector<real,d> Closest_Point_On_Segment(const Vector<real,d>& p,const Vector<real,d>& v1,const Vector<real,d>& v2,real& t);
	template<int d> real Distance_From_Point_To_Segment(const Vector<real,d>& p,const Vector<real,d>& v1,const Vector<real,d>& v2,real& t);

    ////Array operations
	template<class T_ARRAY, class T> void Fill(T_ARRAY& arr, const T& value) { std::fill(arr.begin(), arr.end(), value); }
    template<class T> size_type Find(Array<T>& array,const T& value){auto iter=std::find(array.begin(),array.end(),value);if(iter==array.end())return -1;else return (size_type)(iter-array.begin());}
    template<class T,int d> inline Vector<T,d> Cwise_Prod(const Vector<T,d>& v0,const Vector<T,d>& v1){return v0.cwiseProduct(v1);}
    template<class T,int d> inline Vector<T,d> Cwise_Dvid(const Vector<T,d>& v0,const Vector<T,d>& v1){return v0.cwiseQuotient(v1);}
    template<class T,int dim> inline Vector<T,dim> Cwise_Min(const Vector<T,dim>& v1,const Vector<T,dim>& v2){return v1.cwiseMin(v2);}
    template<class T,int dim> inline Vector<T,dim> Cwise_Max(const Vector<T,dim>& v1,const Vector<T,dim>& v2){return v1.cwiseMax(v2);}
    template<class T> T Min(const Array<T>& v){T v_min=v[0];for(auto i=1;i<v.size();i++)if(v[i]<v_min)v_min=v[i];return v_min;}
    template<class T> T Max(const Array<T>& v){T v_max=v[0];for(auto i=1;i<v.size();i++)if(v[i]>v_max)v_max=v[i];return v_max;}
    template<class T> void Min_And_Max(const Array<T>& v,T& v_min,T& v_max){v_min=v[0];v_max=v[0];for(auto i=1;i<v.size();i++){if(v[i]<v_min)v_min=v[i];if(v[i]>v_max)v_max=v[i];}}
    template<class T> T Mean(const Array<T>& v){T c=Zero<T>();for(auto i=0;i<v.size();i++)c+=v[i];c/=(real)v.size();return c;}
	template<class T> T Sum(const Array<T>& v) { T c = Zero<T>(); for (auto i = 0; i < v.size(); i++)c += v[i]; return c; }
	template<class T,int dim> inline T Min(const Vector<T,dim>& v){return v.minCoeff();}
    template<class T,int dim> inline T Max(const Vector<T,dim>& v){return v.maxCoeff();}
    template<class T,int dim> inline void Min_And_Max(const Vector<T,dim>& v,T& v_min,T& v_max){v_min=v[0];v_max=v[0];for(int i=1;i<dim;i++){if(v[i]<v_min)v_min=v[i];if(v[i]>v_max)v_max=v[i];}}
    template<class T,int dim> inline void Min_And_Max(const Array<Vector<T,dim> >& a,Vector<T,dim>& v_min,Vector<T,dim>& v_max)
    {v_min=a[0];v_max=a[0];for(auto i=1;i<a.size();i++){v_min=Cwise_Min(a[i],v_min);v_max=Cwise_Max(a[i],v_max);}}
    template<class T,int dim> inline T Abs_Min(const Vector<T,dim>& v){T v_min=abs(v[0]);for(int i=1;i<dim;i++)if(abs(v[i])<v_min)v_min=v[i];return v_min;}
    template<class T,int dim> inline T Abs_Max(const Vector<T,dim>& v){T v_max=abs(v[0]);for(int i=1;i<dim;i++)if(abs(v[i])>v_max)v_max=v[i];return v_max;}
    template<class T,int dim> inline int Min_Index(const Vector<T,dim>& v){int i_min=0;for(int i=1;i<dim;i++)if(v[i]<v[i_min])i_min=i;return i_min;}
    template<class T,int dim> inline int Max_Index(const Vector<T,dim>& v){int i_max=0;for(int i=1;i<dim;i++)if(v[i]>v[i_max])i_max=i;return i_max;}
    template<class T> inline T Min(const T& v0,const T& v1){if(v0<v1)return v0;else return v1;}
    template<class T> inline T Max(const T& v0,const T& v1){if(v0>v1)return v0;else return v1;}

	template<class T> real Norm(const T& v);

	template<class ARRAY_T> real Linf(const ARRAY_T& x)
	{real Linf=(real)0;for(auto i=0;i<x.size();i++){real L=(real)abs(x[i]);if(L>Linf)Linf=L;}return Linf;}
	
	template<class T> real Linf_Norm(const Array<T>& x)	
	{real Linf=(real)0;for(const auto& xi:x){real L=Norm<T>(xi);if(L>Linf)Linf=L;}return Linf;}

	template<class T> real Difference_Linf(const Array<T>& x,const Array<T>& x_old)
	{real Linf=(real)0;for(auto i=0;i<x.size();i++){T dif=x[i]-x_old[i];real dif_abs=Norm<T>(dif);if(dif_abs>Linf)Linf=dif_abs;}return Linf;}

    template<class T,int dim> int Min_Index(const ArrayF<T,dim>& v){int i_min=0;for(int i=1;i<dim;i++)if(v[i]<v[i_min])i_min=i;return i_min;}
    template<class T,int dim> int Max_Index(const ArrayF<T,dim>& v){int i_max=0;for(int i=1;i<dim;i++)if(v[i]>v[i_max])i_max=i;return i_max;}
    template<class T,int dim> void Min_And_Max_Index(const ArrayF<T,dim>& v,int& i_min,int& i_max)
    {i_min=0;i_max=0;for(int i=1;i<dim;i++){if(v[i]<v[i_min])i_min=i;if(v[i]>v[i_max])i_max=i;}}

    template<class T,int dim> int Abs_Min_Index(const Vector<T,dim>& v){int i_min=0;for(int i=1;i<dim;i++)if(abs(v[i])<abs(v[i_min]))i_min=i;return i_min;}
    template<class T,int dim> int Abs_Max_Index(const Vector<T,dim>& v){int i_max=0;for(int i=1;i<dim;i++)if(abs(v[i])>abs(v[i_max]))i_max=i;return i_max;}

    template<class T,int dim> bool Has(const Vector<T,dim>& v,const T value){for(int i=0;i<dim;i++)if(v[i]==value)return true;return false;}
    template<class T> bool Has(const Array<T>& v,const T value){for(int i=0;i<(int)v.size();i++)if(v[i]==value)return true;return false;}
    template<class T,int dim> bool Has(const ArrayF<T,dim>& v,const T value){for(int i=0;i<dim;i++)if(v[i]==value)return true;return false;}

	template<class ArrayT,class TV> void Copy_Scalar_Array_To_Vector_Array(const ArrayT& from,const int d,Array<TV>& to) ////Array<real> to Array<TV>
	{
		const int n=(int)from.size()/d;to.resize(n);
		#pragma omp parallel for
		for(int i=0;i<n;i++)for(int j=0;j<d;j++)to[i][j]=from[i*d+j];
	}

	template<class ArrayT,class TV> void Copy_Scalar_Array_To_Vector_Array(const ArrayT& from,const int size,const int d,Array<TV>& to)	////Array<real> to Array<TV>
	{
		const int n=size/d;
		#pragma omp parallel for
		for(auto i=0;i<n;i++)for(int j=0;j<d;j++)to[i][j]=from[i*d+j];
	}

    template<class ArrayT,class TV> void Copy_Indirect(const ArrayT& from,const int d,const Array<int>& map,/*rst*/Array<TV>& to)
    {
		const int n=(int)from.size()/d;
		#pragma omp parallel for
		for(auto i=0;i<n;i++)for(int j=0;j<d;j++)to[map[i]][j]=from[i*d+j];
	}

    template<class T> void Copy_Indirect(const Array<T>& from,const Array<int>& map,/*rst*/Array<T>& to)
    {int c=0;for(auto i=0;i<map.size();i++)if(map[i]>=0)c++;to.resize(c);for(auto i=0;i<from.size();i++)if(map[i]>=0)to[map[i]]=from[i];}
    template<class VecT,class ArrayT> VecT Get_Indirect(const ArrayT& from,const int d,const int idx)
    {VecT rst;for(auto i=0;i<d;i++)rst[i]=from[idx*d+i];return rst;}
    template<class VecT,class ArrayT> void Set_Indirect(ArrayT& to,const VecT& value,const int d,const int idx)
    {for(auto i=0;i<d;i++)to[idx*d+i]=value[i];}

	////Dim conversion for vectors and vector arrays
	template<class T,int d1,int d2> void Dim_Conversion(const Vector<T,d1>& input,Vector<T,d2>& output,const T filled_value=(T)0)
	{
		constexpr int n=d1<d2?d1:d2;
		for(int i=0;i<n;i++)output[i]=input[i];
		if /*constexpr*/ (n<d2){
			for(int i=n;i<d2;i++)output[i]=filled_value;}
	}

	template<class T,int d1,int d2> void Dim_Conversion_Array(const Array<Vector<T,d1> >& input,Array<Vector<T,d2> >& output,const T filled_value=(T)0)
	{
		const int n=(int)input.size();
		#pragma omp parallel for
		for(auto i=0;i<n;i++){Dim_Conversion<T,d1,d2>(input[i],output[i],filled_value);}
	}
	
	////Dim conversion for matrices and matrix arrays
	template<class T,int d1,int d2> void Dim_Conversion(const Matrix<T,d1>& input,Matrix<T,d2>& output,const T filled_value=(T)0)
	{
		constexpr int n=d1<d2?d1:d2;
		output=Matrix<T,d2>::Constant(filled_value);
		for(int i=0;i<n;i++)for(int j=0;j<n;j++)output(i,j)=input(i,j);
	}

	template<class T,int d1,int d2> void Dim_Conversion_Array(const Array<Matrix<T,d1> >& input,Array<Matrix<T,d2> >& output)
	{
		const int n=(int)input.size();
		#pragma omp parallel for
		for(auto i=0;i<n;i++)Dim_Conversion<T,d1,d2>(input[i],output[i]);
	}

	////create vectors with compatible dimensions
	template<int d> Vector<real,d> V(const real x=(real)0,const real y=(real)0,const real z=(real)0);
	template<int d> Vector<int, d> Vi(const int x=0,const int y=0,const int z=0,const int w=0);
	template<int d> Vector<real,d> V(const Vector2& v2);
	template<int d> Vector<real,d> V(const Vector3& v2);

	////eigenvectors, svd, polar decomposition
	////eigenvector/value corresponding to the eigenvalue with the max magnitude
	Vector2 Principal_Eigenvector(const Matrix2& v); 
	Vector3 Principal_Eigenvector(const Matrix3& v);
	real Principal_Eigenvalue(const Matrix2& v); 
	real Principal_Eigenvalue(const Matrix3& v); 

	////eigenvector/value corresponding to the eigenvalue with the min magnitude
	Vector2 Min_Eigenvector(const Matrix2& v);
	Vector3 Min_Eigenvector(const Matrix3& v);
	real Min_Eigenvalue(const Matrix2& v);
	real Min_Eigenvalue(const Matrix3& v);

	template<int dim> void Svd(const Matrix<real,dim>& A,Vector<real,dim>& D,Matrix<real,dim>& U,Matrix<real,dim>& V);
	template<int dim> void Svd_Rot_UV(const Matrix<real,dim>& A,Vector<real,dim>& D,Matrix<real,dim>& U,Matrix<real,dim>& V);
	template<int dim> void Polar_Decomposition(const Matrix<real,dim>& A,Matrix<real,dim>& R,Matrix<real,dim>& S);

	////Verbose
	void Seperation(int n=2);
	void Seperation_Line(int n=2);

	////String Operation
	Array<std::string> Split_String(const std::string& s, const std::string& delimiters = " \t\n\v\f\r");

	////Debug
	void Crash_With_Info(const std::string& s, int ret = 1);
	void Assert(bool flag, const std::string& s = "", int ret = 1);
};

////System paths with values passed from CMake
namespace Path{
	std::string Script();
	std::string Data();
};

#endif
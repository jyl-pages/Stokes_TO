//////////////////////////////////////////////////////////////////////////
// Aux functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "AuxFunc.h"
#include "Constants.h"
#include <iostream>
#include <cmath>

namespace AuxFunc{
	template<> bool Outside<1>(const Vector1i& domain,const Vector1i& p){return p[0]<0||p[0]>=domain[0];}
	template<> bool Outside<2>(const Vector2i& domain,const Vector2i& p){return p[0]<0||p[0]>=domain[0]||p[1]<0||p[1]>=domain[1];}
	template<> bool Outside<3>(const Vector3i& domain,const Vector3i& p){return p[0]<0||p[0]>=domain[0]||p[1]<0||p[1]>=domain[1]||p[2]<0||p[2]>=domain[2];}

	real Quick_Pow(real a, int n)
	{
		real ans = 1;
		while(n){
			if (n & 1) ans *= a;
			n >>= 1;
			a *= a;
		}
		return ans;
	}

	Vector1 Cross(const Vector2& v1,const Vector2& v2){Vector1 v;v[0]=(v1[0]*v2[1]-v1[1]*v2[0]);return v;}
	Vector2 Cross(const Vector1& v1,const Vector2& v2){return Orthogonal_Vector(v2)*v1[0];}
	Vector2 Cross(const real& v1,const Vector2& v2){return Orthogonal_Vector(v2)*v1;}

	Matrix3 Skew(const Vector3& v)
	{
		Matrix3 m;m<<0.,-v[2],v[1], v[2],0.,-v[0], -v[1],v[0],0.;return m;
	}

	template<int d> real Norm(const Vector<real, d>& v) { return v.norm(); }
    template<> real Norm<1>(const Vector1& v){return abs(v[0]);}
    template real Norm<2>(const Vector2& v);
    template real Norm<3>(const Vector3& v);
	
    Vector1 Orthogonal_Vector(const Vector1& v){return v;}
    Vector2 Orthogonal_Vector(const Vector2& v){return {-v.y(),v.x()};}
    Vector3 Orthogonal_Vector(const Vector3& v)
    {
        real abs_x=abs(v.x()),abs_y=abs(v.y()),abs_z=abs(v.z());
        if(abs_x<abs_y) return abs_x<abs_z?Vector3((real)0,v.z(),-v.y()):Vector3(v.y(),-v.x(),(real)0);
        else return abs_y<abs_z?Vector3(-v.z(),(real)0,v.x()):Vector3(v.y(),-v.x(),(real)0);
    }
    Vector4 Orthogonal_Vector(const Vector4& v)
    {
        Vector4 n=Vector4::Zero();int min_axis=0;real min_abs=abs(v[0]);for(int i=1;i<4;i++)if(abs(v[i])<min_abs){min_abs=abs(v[i]);min_axis=i;}
        Vector3 v3;{int c=0;for(int i=0;i<4;i++){if(i==min_axis)continue;v3[c++]=v[i];}}
        v3=Orthogonal_Vector(v3);{int c=0;for(int i=0;i<4;i++){if(i==min_axis)continue;n[i]=v3[c++];}}return n;
    }

    template<int d> Vector<real,d> Unit_Orthogonal_Vector(const Vector<real,d>& v)
    {return Orthogonal_Vector(v).normalized();}
	template Vector2 Unit_Orthogonal_Vector<2>(const Vector2& v);
    template Vector3 Unit_Orthogonal_Vector<3>(const Vector3& v);
    template Vector4 Unit_Orthogonal_Vector<4>(const Vector4& v);

	real Angle_Between(const Vector2& v1,const Vector2& v2)
	{
		real c=v1.dot(v2);return acos(c);
	}

	real Angle_Between(const Vector3& v1,const Vector3& v2)
	{
		real c=v1.dot(v2);return acos(c);
	}

	real Angle_From_To(const Vector2& v1,const Vector2& v2)
	{
		real s=Cross(v1,v2)[0];
		real c=v1.dot(v2);
		return atan2(s,c);
	}

	real Angle_From_To(const Vector3& v1,const Vector3& v2)
	{
		real s=(v1.cross(v2)).norm();
		real c=v1.dot(v2);
		return atan2(s,c);
	}

	real Angle_From_To_With_Specified_Axis(const Vector3& v1,const Vector3& v2,const Vector3& axis)
	{
		Vector3 a;real theta;Angle_And_Axis_Between(v1,v2,theta,a);
		if(a.dot(axis)>0)return theta;else return two_pi-theta;
	}

	void Angle_And_Axis_Between(const Vector3d& v1,const Vector3d& v2,double& theta,Vector3d& axis)
	{
		Vector3d a=v1.cross(v2);double s=a.norm();
		double c=v1.dot(v2);theta=atan2(s,c);axis=a/s;
	}

	void Angle_And_Axis_Between(const Vector3f& v1,const Vector3f& v2,float& theta,Vector3f& axis)
	{
		Vector3f a=v1.cross(v2);float s=a.norm();
		float c=v1.dot(v2);theta=atan2f(s,c);
		if (s == 0) axis = v1;
		else axis = a / s;
	}

	real Dihedral_Angle(const Vector2& a, const Vector2& c, const Vector2& b, const Vector2& d)
	{
		Vector1 n_abc=Cross(Vector2(a-c),Vector2(b-c));
		Vector1 n_acd=Cross(Vector2(d-c),Vector2(a-c));
		return n_abc[0]*n_acd[0]>=0?pi:(real)0;
	}

	real Dihedral_Angle(const Vector3& a,const Vector3& c,const Vector3& b,const Vector3& d)
	{
		Vector3 n_abc=(a-c).cross(b-c).normalized();
		Vector3 n_acd=(d-c).cross(a-c).normalized();
		Vector3 ac=c-a;
		Vector3 n_axis;real theta;Angle_And_Axis_Between(n_abc,n_acd,theta,n_axis);
		if(ac.dot(n_axis)>(real)0)return pi-theta;else return pi+theta;
	}

	template<int d> Vector<real, d> Closest_Point_On_Segment(const Vector<real, d>& p, const Vector<real, d>& v1, const Vector<real, d>& v2, real& t)
	{                  
		Vector<real,d> v=v2-v1;
		real denominator=v.squaredNorm();
		if(denominator==(real)0){t=(real)0;return v1;}
		else{
			t=(p-v1).dot(v)/denominator;
			if(t<=(real)0) return v1;
			else if(t>=1) return v2;
			else{v=v1+(v2-v1)*t;return v;}}
	}
	template Vector2 Closest_Point_On_Segment<2>(const Vector2&,const Vector2&,const Vector2&,real&);
	template Vector3 Closest_Point_On_Segment<3>(const Vector3&,const Vector3&,const Vector3&,real&);

	template<int d> real Distance_From_Point_To_Segment(const Vector<real,d>& p,const Vector<real,d>& v1,const Vector<real,d>& v2,real& t)
	{
		Vector<real,d> v=Closest_Point_On_Segment<d>(p,v1,v2,t);
		return (p-v).norm();
	}
	template real Distance_From_Point_To_Segment<2>(const Vector2&,const Vector2&,const Vector2&,real&);
	template real Distance_From_Point_To_Segment<3>(const Vector3&,const Vector3&,const Vector3&,real&);
	
	Vector2 Principal_Eigenvector(const Matrix2& v)
	{Eigen::SelfAdjointEigenSolver<Matrix2> eig(v);return eig.eigenvectors().col(Abs_Max_Index(eig.eigenvalues()));}
	Vector3 Principal_Eigenvector(const Matrix3& v)
	{Eigen::SelfAdjointEigenSolver<Matrix3> eig(v);return eig.eigenvectors().col(Abs_Max_Index(eig.eigenvalues()));}
	real Principal_Eigenvalue(const Matrix2& v)
	{Eigen::SelfAdjointEigenSolver<Matrix2> eig(v);return Abs_Max(eig.eigenvalues());}
	real Principal_Eigenvalue(const Matrix3& v)
	{Eigen::SelfAdjointEigenSolver<Matrix3> eig(v);return Abs_Max(eig.eigenvalues());}

	Vector2 Min_Eigenvector(const Matrix2& v)
	{Eigen::SelfAdjointEigenSolver<Matrix2> eig(v);return eig.eigenvectors().col(Abs_Min_Index(eig.eigenvalues()));}
	Vector3 Min_Eigenvector(const Matrix3& v)
	{Eigen::SelfAdjointEigenSolver<Matrix3> eig(v);return eig.eigenvectors().col(Abs_Min_Index(eig.eigenvalues()));}
	real Min_Eigenvalue(const Matrix2& v)
	{Eigen::SelfAdjointEigenSolver<Matrix2> eig(v);return Abs_Min(eig.eigenvalues());}
	real Min_Eigenvalue(const Matrix3& v)
	{Eigen::SelfAdjointEigenSolver<Matrix3> eig(v);return Abs_Min(eig.eigenvalues());}

	template<int dim> void Svd(const Matrix<real,dim>& A,Vector<real,dim>& D,Matrix<real,dim>& U,Matrix<real,dim>& V)
	{Eigen::JacobiSVD<Matrix<real,dim> > svd(A,Eigen::ComputeFullU|Eigen::ComputeFullV);D=svd.singularValues();U=svd.matrixU();V=svd.matrixV();}
	template void Svd<2>(const Matrix2&,Vector2&,Matrix2&,Matrix2&);
	template void Svd<3>(const Matrix3&,Vector3&,Matrix3&,Matrix3&);

	template<int dim> void Svd_Rot_UV(const Matrix<real,dim>& A,Vector<real,dim>& D,Matrix<real,dim>& U,Matrix<real,dim>& V)
	{Eigen::JacobiSVD<Matrix<real,dim> > svd(A,Eigen::ComputeFullU|Eigen::ComputeFullV);D=svd.singularValues();U=svd.matrixU();V=svd.matrixV();
	real det_U=U.determinant();real det_V=V.determinant();int flip_n=0;if(det_U<(real)0){U.col(dim-1)*=(real)-1;flip_n++;}if(det_V<(real)0){V.col(dim-1)*=(real)-1;flip_n++;}if(flip_n==1)D[dim-1]*=(real)-1;}
	template void Svd_Rot_UV<2>(const Matrix2&,Vector2&,Matrix2&,Matrix2&);
	template void Svd_Rot_UV<3>(const Matrix3&,Vector3&,Matrix3&,Matrix3&);

	template<int dim> void Polar_Decomposition(const Matrix<real,dim>& A,Matrix<real,dim>& R,Matrix<real,dim>& S)
	{Matrix<real,dim> U,V;Vector<real,dim> D;Svd<dim>(A,D,U,V);R=U*V.transpose();S=V*D.asDiagonal()*V.transpose();}
	template void Polar_Decomposition<2>(const Matrix2&,Matrix2&,Matrix2&);
	template void Polar_Decomposition<3>(const Matrix3&,Matrix3&,Matrix3&);

	void Seperation(int n) { for (int i = 0; i < n; i++)std::cout << "----------------"; }
	void Seperation_Line(int n){std::cout<<"\n";Seperation(n);std::cout<<"\n";}

	Array<std::string> Split_String(const std::string& s, const std::string& delimiters)
	{
		//https://stackoverflow.com/questions/26328793/how-to-split-string-with-delimiter-using-c
		Array<std::string> tokens; tokens.clear();
		std::string::size_type lastPos = s.find_first_not_of(delimiters, 0);
		std::string::size_type pos = s.find_first_of(delimiters, lastPos);
		while (std::string::npos != pos || std::string::npos != lastPos) {
			tokens.push_back(s.substr(lastPos, pos - lastPos));
			lastPos = s.find_first_not_of(delimiters, pos);
			pos = s.find_first_of(delimiters, lastPos);
		}
		return tokens;
	}

	void Crash_With_Info(const std::string& s, int ret)
	{
		std::cerr << s << "\n";
		exit(ret);
	}

	void Assert(bool flag, const std::string& s, int ret)
	{
		if (!flag) Crash_With_Info(s, ret);
	}

	template<> Vector1 V<1>(const real x,const real y,const real z){return Vector1(x);}
	template<> Vector2 V<2>(const real x,const real y,const real z){return Vector2(x,y);}
	template<> Vector3 V<3>(const real x,const real y,const real z){return Vector3(x,y,z);}	
	
	template<> Vector<real,1> V<1>(const Vector2& v2){return Vector1(v2[0]);}
	template<> Vector<real,2> V<2>(const Vector2& v2){return v2;}
	template<> Vector<real,3> V<3>(const Vector2& v2){return Vector3(v2[0],v2[1],(real)0);}
	template<> Vector<real,1> V<1>(const Vector3& v3){return Vector1(v3[0]);}
	template<> Vector<real,2> V<2>(const Vector3& v3){return Vector2(v3[0],v3[1]);}
	template<> Vector<real,3> V<3>(const Vector3& v3){return v3;}

	template<> Vector1i Vi<1>(const int x,const int y,const int z,const int w){return Vector1i(x);}
	template<> Vector2i Vi<2>(const int x,const int y,const int z,const int w){return Vector2i(x,y);}
	template<> Vector3i Vi<3>(const int x,const int y,const int z,const int w){return Vector3i(x,y,z);}
	template<> Vector4i Vi<4>(const int x,const int y,const int z,const int w){return Vector4i(x,y,z,w);}

	template<class T> real Norm(const T& v){return v.norm();}
	template<> real Norm<real>(const real& v){return abs(v);}
	template real Norm<Vector2>(const Vector2&);
	template real Norm<Vector3>(const Vector3&);

}

namespace Path
{
#define STRINGIZE(x) #x
	std::string Script()
	{
#ifdef SCRIPT_PATH
		return STRINGIZE(SCRIPT_PATH);
#else
		return "../../../../script/";	////assuming the current path is in the exe folder
#endif	
	}

	std::string Data()
	{
#ifdef DATA_PATH
		return STRINGIZE(DATA_PATH);
#else
		return "../../../../data/";	////assuming the current path is in the exe folder
#endif
	}
}
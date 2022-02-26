#ifndef __NonlinearFem_h__
#define __NonlinearFem_h__
#include "Common.h"
#include <iostream>

template<int d> class NonlinearFemFunc
{Typedef_VectorDii(d);Typedef_MatrixD(d);Typedef_VectorEi(d+1);
public:
	static void First_PK_Stress_To_Cauchy_Stress(const MatrixD& P,const MatrixD& F,/*rst*/MatrixD& C,const real* J=0) 
	{real det_F=J==0?abs(F.determinant()):*J;C=(real)1/det_F*P*F.transpose();}

	static void Cauchy_Stress_To_First_PK_Stress(const MatrixD& C,const MatrixD& F,/*rst*/MatrixD& P,const real* J=0) 
	{real det_F=J==0?abs(F.determinant()):*J;P=det_F*P*F.inverse().transpose();}

	static void First_PK_Stress_Linear(const real mu,const real lambda,const MatrixD& F,/*rst*/MatrixD& P)
	{P=mu*(F+F.transpose()-(real)2*MatrixD::Identity())+lambda*(F-MatrixD::Identity()).trace()*MatrixD::Identity();}
	static void Cauchy_Stress_Linear(const real mu,const real lambda,const MatrixD& F,/*rst*/MatrixD& C)
	{MatrixD P;First_PK_Stress_Linear(mu,lambda,F,P);First_PK_Stress_To_Cauchy_Stress(P,F,C);}
	
	static void First_PK_Stress_Neohookean(const real mu,const real lambda,const MatrixD& F,/*rst*/MatrixD& P) 
	{real J=F.determinant();MatrixD F_invT=F.inverse().transpose();P=mu*(F-F_invT)+lambda*log(J)*F_invT;}

	static real Elastic_Energy_Neohookean(const real mu,const real lambda,const MatrixD&F)
	{real I1=(F.transpose()*F).trace();real J=F.determinant();real log_J=log(J);return (real).5*mu*(I1-(real)3)-mu*log_J+(real).5*lambda*log_J*log_J;}

	//////////////////////////////////////////////////////////////////////////
	////Hex element
	static void Deformation_Gradient(const real dx,const VectorX& cell_u,const VectorD& natural_coord,/*rst*/MatrixD& F);
	////for implicit solver: 
	////assembling K
	static void dPdx_Neohookean(const real mu,const real lambda,const VectorD& natural_coord,const MatrixD& F,const VectorD& dNidX_T,/*rst*/Array<MatrixD>& dPdx);
	////updating K
	static void Cell_Stiffness_Matrix_And_f_Nonlinear(const real mu,const real lambda,const real dx,const VectorX& cell_u,/*rst*/MatrixX& K_e,/*rst*/VectorX* f_elas=0,/*rst*/real* e_elas=0);

	//////////////////////////////////////////////////////////////////////////
	////Tet element
	static void D(const Vector2& x1,const Vector2& x2,const Vector2& x3,/*rst, ds or dm*/Matrix2& ds)
	{ds.col(0)=x1-x3;ds.col(1)=x2-x3;}
	static void D(const Vector3& x1,const Vector3& x2,const Vector3& x3,const Vector3& x4,/*rst*/Matrix3& ds) 
	{ds.col(0)=x1-x4;ds.col(1)=x2-x4;ds.col(2)=x3-x4;}

	static void D(const ArrayF<Vector<real,2>,3>& x,Matrix<real,2>& ds)
	{return D(x[0],x[1],x[2],ds);}
	static void D(const ArrayF<Vector<real,3>,4>& x,Matrix<real,3>& ds)
	{return D(x[0],x[1],x[2],x[3],ds);}

	static void D_Inv_And_Vol(const Vector2& X1,const Vector2& X2,const Vector2& X3,/*rst*/Matrix2& dm_inv,/*rst*/real& vol)
	{D(X1,X2,X3,dm_inv);vol=(real).5*dm_inv.determinant();dm_inv=dm_inv.inverse().eval();}
	static void D_Inv_And_Vol(const Vector3& X1,const Vector3& X2,const Vector3& X3,const Vector3& X4,/*rst*/Matrix3& dm_inv,/*rst*/real& vol)
	{D(X1,X2,X3,X4,dm_inv);vol=(real)1/(real)6*abs(dm_inv.determinant());dm_inv=dm_inv.inverse().eval();}

	static void D_Inv_And_Vol(const ArrayF<Vector<real,2>,3>& x,/*rst*/Matrix2& dm_inv,/*rst*/real& vol)
	{D_Inv_And_Vol(x[0],x[1],x[2],dm_inv,vol);}
	static void D_Inv_And_Vol(const ArrayF<Vector<real,3>,4>& x,/*rst*/Matrix3& dm_inv,/*rst*/real& vol)
	{D_Inv_And_Vol(x[0],x[1],x[2],x[3],dm_inv,vol);}

	static void Elastic_Force_Neohookean(const Vector3& x1,const Vector3& x2,const Vector3& x3,const Vector3& x4,const Matrix3& Dm_inv,const real vol,const real mu,const real lambda,/*rst*/Matrix3& H)
	{Matrix3 _Ds,F,P;D(x1,x2,x3,x4,_Ds);F=_Ds*Dm_inv;NonlinearFemFunc<3>::First_PK_Stress_Neohookean(mu,lambda,F,P);H=-vol*P*Dm_inv.transpose();}
};

#endif
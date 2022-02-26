//////////////////////////////////////////////////////////////////////////
// Krylov solver
// Copyright(c)(2018-),Bo Zhu,Jinyuan Liu,Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//NOTE: DO NOT instantiate this file in any .cu!
//lu()in GMRES()will cause confusion for NVCC.
//This is why I moved that function to KrylovSolver.cpp
// -- Mengdi Wang


#ifndef __KrylovSolver_h__
#define __KrylovSolver_h__
#include <iostream>
#include "Common.h"
#include "SparseFunc.h"
#include <Eigen/Dense>

////A sample for Krylov vector base:
////Member functions to implement:
//void Axpy(real a,const KrylovX<KX>& x,const KrylovX<KX>& y);
//void Copy(const KrylovX<KX>& x);
//real Dot(const KrylovX<KX>& x);
//void Set_Zero();

template<class KX> class KrylovX
{using KXW=KrylovX<KX>;
public:
	KX* kx=nullptr;
	bool own_data=true;
	KrylovX(KX* _kx):kx(_kx){own_data=false;}
	KrylovX(const KX& _kx){if(kx==nullptr)kx=new KX();*kx=_kx;}
	KXW& operator=(const KXW& copy){if(kx==nullptr)kx=new KX();*kx=*copy.kx;return *this;}
	KrylovX(const KXW& copy){*this=copy;}
	KrylovX(int n){if(kx!=nullptr&&own_data)delete kx;kx=new KX(n);Set_Zero();}
	KrylovX(){if(kx==nullptr)kx=new KX();}
	~KrylovX(){if(own_data){delete kx;}}
	
	////*this=ax+y
	void Axpy(real a,const KXW& x,const KXW& y)
	{*kx=a*(*x.kx)+(*y.kx);}		
	////*this=x
	void Copy(const KXW& x)
	{*kx=*x.kx;}
	////(*this).dot(x)
	real Dot(const KXW& x)const
	{return(*kx).dot(*x.kx);}	
	////fill all the elements to be zero
	void Set_Zero()
	{kx->fill((real)0);}

	void Print(){}
};

////A sample for Krylov matrix base:
////Member functions:
////void Compute_Precond();
////template<class KX> void Multiply(const KX& x,KX& rst);
////template<class KX> void Res(const KX& b,const KX& x,KX& rst);

template<class KA> class KrylovA
{public:
	const KA* ka=nullptr;
	Eigen::IncompleteCholesky<real,Eigen::Upper,Eigen::NaturalOrdering<int> > ic;

	KrylovA(const KA* _ka):ka(_ka){Compute_Precond();}

	void Compute_Precond(){ic.compute(*ka);}

	////rst=Ax
	template<class KXW> void Multiply(const KXW& x,KXW& rst)const 
	{*rst.kx=(*ka)*(*x.kx);}							
	////rst=b-Ax
	template<class KXW> void Res(const KXW& b,const KXW& x,KXW& rst)const 
	{*rst.kx=*b.kx-(*ka)*(*x.kx);}			

	void Print(){}
};

////A sample for Krylov preconditioner
////Member functions:
////template<class KXW> void Precond(const KXW& r,KXW& rst)

template<class KA> class KrylovPrecIdentity
{public:
	const KA* ka=nullptr;
	KrylovPrecIdentity(const KA* _ka):ka(_ka){Initialize();}
	void Initialize(){}
	template<class KXW> void Precond(const KXW& r,KXW& rst){rst=r;}
};

template<class KA> class KrylovPrecIC
{public:
	const KA* ka=nullptr;
	Eigen::IncompleteCholesky<real,Eigen::Upper,Eigen::NaturalOrdering<int> > ic;

	KrylovPrecIC(const KA* _ka):ka(_ka){Initialize();}

	void Initialize()
	{
		ic.compute(*ka);
	}

	template<class KXW> void Precond(const KXW& r,KXW& rst)
	{
		*rst.kx=ic.solve(*r.kx);
		if(ic.info()!=Eigen::Success){std::cerr<<"Precond failed";}
	}
};

namespace KrylovSolver
{
	class Params
	{public:
		real tolerance=(real)1e-5;
		int max_iter_num=1000;
		bool verbose=true;
		int restart_iter_num=10;
	};

	template<class KA,class KX,class KP> bool Preconditioned_Conjugate_Gradient(const KA& A,KX& x,const KX& b,KX& r,KX& d,KX& q,KX& s,KP& precond,const Params params=Params())
	{
		A.Res(b,x,r);													////r=b-A*x
		precond.Precond(r,d);											////d=precond(r)
		real delta_0=r.Dot(d);
		real delta=delta_0;
		real epsilon_0=pow(params.tolerance,2)*delta_0;
		int iter=1;

		while(iter<params.max_iter_num && delta>epsilon_0){
			A.Multiply(d,q);											////q=A*d
			real alpha=delta/(d.Dot(q));
			x.Axpy(alpha,d,x);											////x=x+alpha*d

			if(iter%params.restart_iter_num==0)A.Res(b,x,r);			////r=b-A*x
			else r.Axpy(-alpha,q,r);									////r=r-alpha*q
			s.Set_Zero();												////set s to zero before precond!
			precond.Precond(r,s);										////s=precond(r)
			real delta_old=delta;
			delta=r.Dot(s);
			real beta=delta/delta_old;
			d.Axpy(beta,d,s);											////d=s+beta*d;
			iter++;}

		bool converged=delta<=epsilon_0;
		
		if(params.verbose){
			real res=delta_0==0?delta:(delta/delta_0);
			std::cout<<"KrylovSolve::PCG: "<<(converged ? "converge" : "fail")
				<<" in "<<(iter-1)<<" iterations,with residual "<<sqrt(res)<<"."<<std::endl;}

		return converged;
	}

	inline bool CG(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params=Params())
	{
		KrylovA<SparseMatrix<real> > kA(&A);
		KrylovX<VectorN<real> > kx(&x),kb(b);
		KrylovPrecIdentity<SparseMatrix<real> > kP(&A);
		int nc=(int)x.rows(); KrylovX<VectorN<real> > kr(nc),kd(nc),kq(nc),ks(nc);
		return Preconditioned_Conjugate_Gradient(kA,kx,kb,kr,kd,kq,ks,kP);
	}

	inline bool ICPCG(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params=Params())
	{
		KrylovA<SparseMatrix<real> > kA(&A);
		KrylovX<VectorN<real> > kx(&x),kb(b);
		KrylovPrecIC<SparseMatrix<real> > kP(&A);
		int nc=(int)x.rows();KrylovX<VectorN<real> > kr(nc),kd(nc),kq(nc),ks(nc);
		return Preconditioned_Conjugate_Gradient(kA,kx,kb,kr,kd,kq,ks,kP,params);
	}

	//////////////////////////////////////////////////////////////////////////
	////GMRES
	bool GMRES(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params=Params());

};

#endif
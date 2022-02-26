//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid
// Copyright (F) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "GmgPcgSolverCPU.h"

using namespace GeometricMultiGrid;

//////////////////////////////////////////////////////////////////////////
////solver object
//////////////////////////////////////////////////////////////////////////

template<int d> void GMGPCG_Solver_CPU<d>::Initialize(const SparseMatrix<real>& A, 
	const Vector<int,d>& counts,const MultiGrid::Params& params,const Field<short,d>* mat_id/*=nullptr*/)
{
	if(initialized&&!update_A_levels){return;}
	kP=std::make_shared<T_GMG>(&A,counts,params,mat_id);
	initialized=true;
}

template<int d> bool GMGPCG_Solver_CPU<d>::Solve(VectorN<real>& x,const VectorN<real>& b)
{
	if(update_A_levels)kP->multigrid.Update_A_Levels();

	KrylovA<SparseMatrix<real> > kA(kP->multigrid.A[0]);
	KrylovX<VectorN<real> > kx(&x),kb(b);
	int nc=(int)x.rows();KrylovX<VectorN<real> > kr(nc),kd(nc),kq(nc),ks(nc);
	KrylovSolver::Params krylov_params;
	krylov_params.tolerance=kP->params.tolerance;
	krylov_params.max_iter_num=kP->params.max_cg_iter_num;
	return KrylovSolver::Preconditioned_Conjugate_Gradient(kA,kx,kb,kr,kd,kq,ks,(*kP),krylov_params);
}

template class GMGPCG_Solver_CPU<2>;
template class GMGPCG_Solver_CPU<3>;

//////////////////////////////////////////////////////////////////////////
////single function API
//////////////////////////////////////////////////////////////////////////

template<int d> bool GMGPCG_CPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int, d>& counts,const Params params/*=Params()*/)
{
	KrylovA<SparseMatrix<real> > kA(&A);
	KrylovPreGMG<d,real,SparseMatrix<real>,VectorN<real> > kP(&A,counts,params);

	KrylovX<VectorN<real> > kx(&x),kb(b);
	int nc=(int)x.rows();KrylovX<VectorN<real> > kr(nc),kd(nc),kq(nc),ks(nc);
	return KrylovSolver::Preconditioned_Conjugate_Gradient(kA,kx,kb,kr,kd,kq,ks,kP);
}

template<int d> bool GMGPCG_CPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const Field<short,d>& mat_id,const Params params/*=Params()*/)
{
	KrylovA<SparseMatrix<real> > kA(&A);
	KrylovPreGMG<d,real,SparseMatrix<real>,VectorN<real> > kP(&A,counts,params,&mat_id);

	KrylovX<VectorN<real> > kx(&x),kb(b);
	int nc=(int)x.rows();KrylovX<VectorN<real> > kr(nc),kd(nc),kq(nc),ks(nc);
	return KrylovSolver::Preconditioned_Conjugate_Gradient(kA,kx,kb,kr,kd,kq,ks,kP);
}

#define Inst_Helper(d) \
template bool GMGPCG_CPU<d>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,d>&,const Params);\
template bool GMGPCG_CPU<d>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,d>&,const Field<short,d>&,const Params);
Inst_Helper(2);
Inst_Helper(3);
#undef Inst_Helper
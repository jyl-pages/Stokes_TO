//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
//// This is the CPU implementation of the geometric multigrid class
//////////////////////////////////////////////////////////////////////////

#ifndef __GmgPcgSolverCPU_h__
#define __GmgPcgSolverCPU_h__
#include "GeometricMultiGrid.h"

//////////////////////////////////////////////////////////////////////////
////CPU geometric multigrid solver object
template<int d> class GMGPCG_Solver_CPU
{
	using T_GMG=GeometricMultiGrid::KrylovPreGMG<d,real,SparseMatrix<real>,VectorN<real> >;
	std::shared_ptr<T_GMG> kP=nullptr;
public:
	bool update_A_levels=false;
	bool initialized=false;

	virtual void Initialize(const SparseMatrix<real>& A,const Vector<int,d>& counts,const MultiGrid::Params& params,const Field<short,d>* mat_id=nullptr);
	virtual bool Solve(VectorN<real>& x,const VectorN<real>& b);
};

//////////////////////////////////////////////////////////////////////////
////CPU geometric multigrid single function API
template<int d> bool GMGPCG_CPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const MultiGrid::Params params=MultiGrid::Params());
template<int d> bool GMGPCG_CPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const Field<short,d>& mat_id,const MultiGrid::Params params=MultiGrid::Params());

#endif

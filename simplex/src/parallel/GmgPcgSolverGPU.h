//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid PCG GPU
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
//// This is the GPU implementation of the geometric multigrid class
//////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA
#ifndef __GmgPcgSolverGPU_h__
#define __GmgPcgSolverGPU_h__

#include "GeometricMultiGrid.h"
#include "MultiGridCuda.h"

//////////////////////////////////////////////////////////////////////////
//// GPU geometric multigrid solver object
//// ATTENTION: this function does not support DYNAMIC, irregular domain (the DoFs on the grid should not change during the solve).
template<class T,int d> class GMGPCG_Solver_GPU
{
	using T_MG_HOST=GeometricMultiGrid::GeometricMultiGridSystem<d,T,SparseMatrix<T>,VectorN<T> >;
	T_MG_HOST* mg_host=nullptr;
	MultiGridCuda::MultiGridSystemCuda<T>* mg_dev=nullptr;
public:
	bool update_A_levels=false;				////The flag specifies whether to recalculate A levels based on A[0]. ATTENTION: it assumes R and P do not change and does not support dynamic domain!
	bool initialized=false;
	GeometricMultiGrid::Params params;

	GMGPCG_Solver_GPU(){}
	~GMGPCG_Solver_GPU();

	////The way to initialize a GPU mg is to initialize a CPU mg data hierarchy based on the input A and mat_id first and then copy the hierarchy to GPU.
	virtual void Initialize(const SparseMatrix<T>& A,const Vector<int,d>& counts,const GeometricMultiGrid::Params& params,const Field<short,d>* mat_id=nullptr);
	virtual bool Solve(VectorN<T>& x,const VectorN<T>& b);

protected:
	T* b_dev=nullptr;
	T* x_dev=nullptr;
	T* p_dev=nullptr;
	T* q_dev=nullptr;
};

//////////////////////////////////////////////////////////////////////////
//// GPU geometric multigrid single function API
template<int d> bool GMGPCG_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const MultiGrid::Params params=MultiGrid::Params(),bool verbose=false);
template<int d> bool GMGPCG_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const Field<short,d>& mat_id,const MultiGrid::Params params=MultiGrid::Params(),bool verbose=false);

#endif
#endif
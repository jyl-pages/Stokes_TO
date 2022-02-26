//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "GmgPcgSolverGPU.h"
#include "ContextCuda.h"
#include "AuxFuncCuda.h"
using namespace AuxFuncCuda;
using namespace MultiGridCuda;

//////////////////////////////////////////////////////////////////////////
//// solver object
template<class T,int d> GMGPCG_Solver_GPU<T,d>::~GMGPCG_Solver_GPU()
{
	delete mg_host;
	delete mg_dev;
	if(b_dev!=nullptr)Global_Free(b_dev, DEVICE);
	if(x_dev!=nullptr)Global_Free(x_dev, DEVICE);
	if(p_dev!=nullptr)Global_Free(p_dev, DEVICE);
	if(q_dev!=nullptr)Global_Free(q_dev, DEVICE);
}

template<class T,int d> void GMGPCG_Solver_GPU<T,d>::Initialize(
	const SparseMatrix<T>& A,const Vector<int,d>& counts,const GeometricMultiGrid::Params& params,const Field<short,d>* mat_id/*=nullptr*/)
{
	if(initialized)return;
	mg_host=new T_MG_HOST();
	mg_host->init_A_levels=!params.init_hier_on_gpu;				////ATTENTION: for GPU multigrid, we calculate RAP on device by default
	mg_host->Initialize(A,counts,params,mat_id);					////
	if(!params.use_color)mg_host->Initialize_Color();				////initialize color anyhow on CPU for GPU. TODO: this can be improved by initializing color on GPU directly
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();
	MultiGridCuda::Malloc_And_Copy_Host_To_Device<T,T_MG_HOST>(mg_dev,mg_host);	////copy CPU data to GPU
	int n=mg_dev->A[0]->m;
	b_dev=Global_Malloc<T>(n, DEVICE);
	x_dev=Global_Malloc<T>(n, DEVICE);
	p_dev=Global_Malloc<T>(n, DEVICE);
	q_dev=Global_Malloc<T>(n, DEVICE);
	initialized=true;	
}

template<class T,int d> bool GMGPCG_Solver_GPU<T,d>::Solve(VectorN<T>& x,const VectorN<T>& b)
{
	if(!initialized||update_A_levels){
		T* dev_A_val=mg_dev->A[0]->valuePtr();
		T* host_A_val=mg_host->A[0]->valuePtr();
		Global_Copy_Array<T>(dev_A_val, host_A_val, mg_host->A[0]->nonZeros(), DEVICE, HOST);
		Update_A_Levels_On_Device(mg_dev,mg_host->one_over_intp);}
	int n=mg_dev->A[0]->m;
	auto* bptr=const_cast<VectorN<real>*>(&b);
	T* b_host=bptr->data();
	T* x_host=x.data();
	Global_Copy_Array<T>(b_dev, b_host, n, DEVICE, HOST);
	Global_Copy_Array<T>(x_dev, x_host, n, DEVICE, HOST);
	GMGPCG_Device(mg_dev,b_dev,x_dev,p_dev,q_dev,params.max_cg_iter_num,params.tolerance);
	Global_Copy_Array<T>(x_host, x_dev, n, HOST, DEVICE);
	return true;
}

template class GMGPCG_Solver_GPU<real,2>;
template class GMGPCG_Solver_GPU<real,3>;

//////////////////////////////////////////////////////////////////////////
//// single function API
////gpu mg on regular domain
template<int d> bool GMGPCG_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const MultiGrid::Params params/*=Params()*/,bool verbose/*=false*/)
{
	using namespace GeometricMultiGrid;
#ifndef USE_CUDA
	std::cerr<<"Error: [GMGPCG_GPU] USE_CUDA disabled"<<std::endl;
	return false;
#else
	Timer<real> timer;timer.Reset();
	using T_MG=GeometricMultiGridSystem<d,real,SparseMatrix<real>,VectorN<real> >;
	T_MG mg_host;
	mg_host.init_A_levels=!params.init_hier_on_gpu;		////ATTENTION: for GPU multigrid, we calculate RAP on device by default
	auto* bptr=const_cast<VectorN<real>*>(&b);
	mg_host.Initialize(A,counts,params);
	if(!params.use_color)mg_host.Initialize_Color();	////initialize color anyhow for GPU
	if(verbose) timer.Elapse_And_Output_And_Reset("MG: build system on host");
	MultiGridCuda::GMGPCG<real,T_MG>(&mg_host,bptr->data(),x.data(),params.max_cg_iter_num,params.tolerance,verbose);
	if(verbose) timer.Elapse_And_Output_And_Reset("MG: solve system on dev");
	return true;
#endif
}

template bool GMGPCG_GPU<2>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,2>&,const MultiGrid::Params,bool verbose);
template bool GMGPCG_GPU<3>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,3>&,const MultiGrid::Params,bool verbose);

////gpu mg on irregular domain
template<int d> bool GMGPCG_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const Field<short,d>& mat_id,const MultiGrid::Params params/*=MultiGrid::Params()*/,bool verbose/*=false*/)
{
	using namespace GeometricMultiGrid;
#ifndef USE_CUDA
	std::cerr<<"Error: [GMGPCG_GPU] USE_CUDA disabled"<<std::endl;
	return false;
#else
	Timer<real> timer;timer.Reset();
	using T_MG=GeometricMultiGridSystem<d,real,SparseMatrix<real>,VectorN<real> >;
	T_MG mg_host;
	mg_host.init_A_levels=false;		////ATTENTION: for GPU multigrid, we calculate RAP on device by default
	auto* bptr=const_cast<VectorN<real>*>(&b);
	mg_host.Initialize(A,counts,params,&mat_id);
	if(!params.use_color)mg_host.Initialize_Color();	////initialize color anyhow for GPU
	if(verbose) timer.Elapse_And_Output_And_Reset("MG: build system on host");
	MultiGridCuda::GMGPCG<real,T_MG>(&mg_host,bptr->data(),x.data(),params.max_cg_iter_num,params.tolerance,verbose);
	if(verbose) timer.Elapse_And_Output_And_Reset("MG: solve system on dev");
	return true;
#endif
}

template bool GMGPCG_GPU<2>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,2>&,const Field<short,2>&,const MultiGrid::Params,bool verbose);
template bool GMGPCG_GPU<3>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,3>&,const Field<short,3>&,const MultiGrid::Params,bool verbose);




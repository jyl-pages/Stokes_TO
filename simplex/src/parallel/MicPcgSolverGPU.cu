//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "MicPcgSolverGPU.h"
#include "ContextCuda.h"
#include "AuxFuncCuda.h"
using namespace AuxFuncCuda;
using namespace MultiGridCuda;

////////////////////////////////////////////////////////////////////////////
////// solver object
//template<class T,int d> MICPCG_Solver_GPU<T,d>::~MICPCG_Solver_GPU()
//{
//	delete mg_host;
//	delete mg_dev;
//	if(b_dev!=nullptr)Free_On_Device(b_dev); 
//	if(x_dev!=nullptr)Free_On_Device(x_dev);
//	if(p_dev!=nullptr)Free_On_Device(p_dev);
//	if(q_dev!=nullptr)Free_On_Device(q_dev);	
//}
//
//template<class T,int d> void MICPCG_Solver_GPU<T,d>::Initialize(
//	const SparseMatrix<T>& A,const Vector<int,d>& counts,const GeometricMultiGrid::Params& params,const Field<short,d>* mat_id/*=nullptr*/)
//{
//	if(initialized)return;
//	mg_host=new T_MG_HOST();
//	mg_host->init_A_levels=false;									////ATTENTION: for GPU multigrid, we calculate RAP on device by default
//	mg_host->Initialize(A,counts,params,mat_id);					////
//	if(!params.use_color)mg_host->Initialize_Color();				////initialize color anyhow on CPU for GPU. TODO: this can be improved by initializing color on GPU directly
//	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();
//	MultiGridCuda::Malloc_And_Copy_Host_To_Device<T,T_MG_HOST>(mg_dev,mg_host);	////copy CPU data to GPU
//	int n=mg_dev->A[0]->m;
//	Malloc_On_Device(b_dev,n);
//	Malloc_On_Device(x_dev,n);
//	Malloc_On_Device(p_dev,n);
//	Malloc_On_Device(q_dev,n);
//	initialized=true;	
//}
//
//template<class T,int d> bool MICPCG_Solver_GPU<T,d>::Solve(VectorN<T>& x,const VectorN<T>& b)
//{
//	if(!initialized||update_A_levels){
//		T* dev_A_val=mg_dev->A[0]->valuePtr();
//		T* host_A_val=mg_host->A[0]->valuePtr();
//		Copy_Array_Host_To_Device(dev_A_val,host_A_val,mg_host->A[0]->nonZeros());
//		Update_A_Levels_On_Device(mg_dev,mg_host->one_over_intp);}
//	int n=mg_dev->A[0]->m;
//	auto* bptr=const_cast<VectorN<real>*>(&b);
//	T* b_host=bptr->data();
//	T* x_host=x.data();
//	Copy_Array_Host_To_Device(b_dev,b_host,n);
//	Copy_Array_Host_To_Device(x_dev,x_host,n);
//	GMGPCG_Device(mg_dev,b_dev,x_dev,p_dev,q_dev,params.max_cg_iter_num,params.tolerance);
//	Copy_Array_Device_To_Host(x_host,x_dev,n);	
//	return true;
//}
//
//template class MICPCG_Solver_GPU<real,2>;
//template class MICPCG_Solver_GPU<real,3>;

//////////////////////////////////////////////////////////////////////////
//// single function API

////gpu micpcg on regular domain
template<int d> bool MICPCG_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const MultiGrid::Params params/*=Params()*/,bool verbose/*=false*/)
{
	using namespace GeometricMultiGrid;
#ifndef USE_CUDA
	std::cerr<<"Error: [MICPCG_GPU] USE_CUDA disabled"<<std::endl;
	return false;
#else
	Timer<real> timer;timer.Reset();
	using T_MG=GeometricMultiGridSystem<d,real,SparseMatrix<real>,VectorN<real> >; 
	T_MG mg_host;

	mg_host.init_A_levels=false;		////ATTENTION: for GPU multigrid, we calculate RAP on device by default
	auto* bptr=const_cast<VectorN<real>*>(&b);
	MultiGrid::Params prm=params;
	prm.use_auto_calculated_levels=false;
	prm.levels=1;
	prm.use_color=false;
	prm.max_cg_iter_num=10000;
	mg_host.Initialize(A,counts,prm);
	MultiGridCuda::MICPCG<real,T_MG>(&mg_host,bptr->data(),x.data(),prm.max_cg_iter_num,prm.tolerance,verbose);
	if(verbose) timer.Elapse_And_Output_And_Reset("MICPCG: solve system on dev");
	return true;
#endif
}

template bool MICPCG_GPU<2>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,2>&,const MultiGrid::Params,bool verbose);
template bool MICPCG_GPU<3>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,3>&,const MultiGrid::Params,bool verbose);

////gpu micpcg on irregular domain
template<int d> bool MICPCG_GPU(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Vector<int,d>& counts,const Field<short,d>& mat_id,const MultiGrid::Params params/*=MultiGrid::Params()*/,bool verbose/*=false*/)
{
	using namespace GeometricMultiGrid;
#ifndef USE_CUDA
	std::cerr<<"Error: [MICPCG_GPU] USE_CUDA disabled"<<std::endl;
	return false;
#else
	Timer<real> timer;timer.Reset();
	using T_MG=GeometricMultiGridSystem<d,real,SparseMatrix<real>,VectorN<real> >;
	T_MG mg_host;
	mg_host.init_A_levels=false;		////ATTENTION: for GPU multigrid, we calculate RAP on device by default
	auto* bptr=const_cast<VectorN<real>*>(&b);
	MultiGrid::Params prm=params;
	prm.use_auto_calculated_levels=false;
	prm.levels=1;
	prm.use_color=false;
	prm.max_cg_iter_num=2000;
	mg_host.Initialize(A,counts,prm,&mat_id);
	MultiGridCuda::MICPCG<real,T_MG>(&mg_host,bptr->data(),x.data(),prm.max_cg_iter_num,prm.tolerance,verbose);
	if(verbose) timer.Elapse_And_Output_And_Reset("MICPCG irregular: solve system on dev");
	return true;
#endif
}

template bool MICPCG_GPU<2>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,2>&,const Field<short,2>&,const MultiGrid::Params,bool verbose);
template bool MICPCG_GPU<3>(const SparseMatrix<real>&,VectorN<real>&,const VectorN<real>&,const Vector<int,3>&,const Field<short,3>&,const MultiGrid::Params,bool verbose);




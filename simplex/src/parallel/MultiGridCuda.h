//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid on GPU
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA
#ifndef __MultiGridCuda_h__
#define __MultiGridCuda_h__
#include "AuxFuncCuda.h"

namespace MultiGridCuda
{
//////////////////////////////////////////////////////////////////////////
//// GPU multigrid data structure: stores a hierarchy of matrices
//// The pointers support both host and device
template<class T> class MultiGridSystemCuda
{public:
	int levels;
	SparseMatrixCuda<T>** A=nullptr;
	SparseMatrixCuda<T>** P=nullptr;
	SparseMatrixCuda<T>** R=nullptr;
	SparseMatrixCuda<T>** AP=nullptr;
	T** b=nullptr;
	T** x=nullptr;
	T** r=nullptr;
	bool on_device;
	//DataHolder data_holder;

	////solve the coarsest level on host
	SparseMatrixCuda<T>* A_coarsest_host=nullptr;	////on host
	T* b_coarsest_host=nullptr;						////on host
	T* x_coarsest_host=nullptr;						////on host
	bool coarsest_on_host=false;					////by default not put coarsest solve on host

	MultiGridSystemCuda(){}
	~MultiGridSystemCuda(){Clear();}

	////A CUDA MG class can be initialized from a host class, e.g., GeometricMultiGridSystem
	template<typename T_MG> void Initialize(T_MG& multi_grid_system,bool _on_device)
	{
		on_device=_on_device;
		if(on_device)Malloc_And_Copy_Host_To_Device(this,&multi_grid_system);
		else Malloc_And_Copy_Host_To_Host(this,&multi_grid_system);
	}

	void Sync_Pointers_And_Flags(int _levels,SparseMatrixCuda<T>**& _A,SparseMatrixCuda<T>**& _P,SparseMatrixCuda<T>**& _R,T**& _b,T**& _x,T**& _r,bool _on_device)
	{levels=_levels;A=_A;P=_P;R=_R;b=_b;x=_x;r=_r;on_device=_on_device;}

	void Clear();
	size_t Memory_Size();
};

////Aux functions
template<typename T> void Malloc_And_Copy_Host_To_Device(MultiGridSystemCuda<T>*& ptr_dev,MultiGridSystemCuda<T>*& ptr_host);
template<typename T, typename T_MG> void Malloc_And_Copy_Host_To_Device(MultiGridSystemCuda<T>*& ptr_dev,T_MG*& ptr_host);
template<typename T, typename T_MG> void Malloc_And_Copy_Host_To_Host(MultiGridSystemCuda<T>*& ptr_dev,T_MG*& ptr_host);

//////////////////////////////////////////////////////////////////////////
////Solver APIs
////Jacobi
template<class T> void Jacobi(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance);
template<class T> void Jacobi_Device(SparseMatrixCuda<T>*& A_dev,T*& b,T*& x_cur,T*& x_next,T*& r,int iter_num,T tolerance);

////G-S
template<class T,int thread_n,int nnz_max> void Multi_Color_Gauss_Seidel(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,int order=0);
template<class T,int thread_n,int nnz_max> void Multi_Color_Gauss_Seidel_Device(SparseMatrixCuda<T>* A_dev,T* b,T* x,T* r,int iter_num,T tolerance,int order=0);
template<class T> void Serialized_Gauss_Seidel(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,int order=0);
template<class T> void Serialized_Gauss_Seidel_Device(SparseMatrixCuda<T>* A_dev,T* b,T* x,T* r,int iter_num,T tolerance,int order=0);

////CG
template<class MtxT,class T> void Preconditioned_Conjugate_Gradient(MtxT* A_host,T* _b,T* _x,int iter_num,T tolerance,bool use_precond=true);

////MG
template<class T> void Multi_Grid(MultiGridSystemCuda<T>* mg_host,int iter_num,T tolerance);
template<class T> void V_Cycle_Device(MultiGridSystemCuda<T>* mg_dev,int iter_num,T tolerance);
template<class T> void Update_A_Levels_On_Device(MultiGridSystemCuda<T>* mg_dev,T scale,bool resize=false);	////update A_{n+1}=R_n*A_n*P_n

////Direct
int Direct_Solve_Device(SparseMatrixCuda<double>* A_dev,double* b,double* x,double tolerance);
int Direct_Solve_Device(SparseMatrixCuda<float>* A_dev,float* b,float* x,float tolerance);
int Direct_Solve_Host(SparseMatrixCuda<double>* A_dev,double* b,double* x,double tolerance);
int Direct_Solve_Host(SparseMatrixCuda<float>* A_dev,float* b,float* x,float tolerance);

////Diagonal preconditioner

//////////////////////////////////////////////////////////////////////////
////Geometric CUDA multigrid API
template<class T,class T_MG> void GMGPCG(T_MG* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose=false);
template<class T> void GMGPCG_Device(MultiGridSystemCuda<T>* mg_dev,T* b,T* x,T* p,T* q,int iter_num,T tolerance,bool verbose=true);

template<class T,class T_MG> void MICPCG(T_MG* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose=false);
template<class T> void MICPCG_Device(MultiGridSystemCuda<T>* mg_dev,T* b,T* x,T* p,T* q,int iter_num,T tolerance,bool verbose=true);
};

#endif
#endif

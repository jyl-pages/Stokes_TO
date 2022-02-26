//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include <iostream>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <cusolverSp.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#include <cuda_profiler_api.h>
#include "ContextCuda.h"
#include "GeometricMultiGrid.h"
#include "MultiGridCuda.h"
#include "Timer.h"
#include "AuxFunc.h"

namespace MultiGridCuda
{
using namespace AuxFuncCuda;
using namespace GeometricMultiGrid;
size_t Mb(){return 1048576;}

//////////////////////////////////////////////////////////////////////////
////Multigrid data structures
template<class T> void MultiGridSystemCuda<T>::Clear()
{
	if(A){for(int i=0;i<levels;i++)delete A[i];delete[] A;}
	if(P){for(int i=0;i<levels-1;i++)delete P[i];delete[] P;}
	if(R){for(int i=0;i<levels-1;i++)delete R[i];delete[] R;}
	if(AP){for(int i=0;i<levels-1;i++)delete AP[i];delete[] AP;}

	if(!on_device){
		if(b){for(int i=0;i<levels;i++)if(b[i])delete[] b[i];delete[] b;}
		if(x){for(int i=0;i<levels;i++)if(x[i])delete[] x[i];delete[] x;}
		if(r){for(int i=0;i<levels;i++)if(r[i])delete[] r[i];delete[] r;}}
	else{
		if(b){for(int i=0;i<levels;i++)Global_Free(b[i],DEVICE);delete[] b;}
		if(x){for(int i=0;i<levels;i++)Global_Free(x[i],DEVICE);delete[] x;}
		if(r){for(int i=0;i<levels;i++)Global_Free(r[i],DEVICE);delete[] r;}}
}

template<class T> size_t MultiGridSystemCuda<T>::Memory_Size()
{
	size_t size=0;
	for(int i=0;i<levels;i++){
		if(A[i])size+=A[i]->Memory_Size();
		if(b[i])size+=A[i]->m*sizeof(T);	////b,x,r
		if(x[i])size+=A[i]->m*sizeof(T);
		if(r[i])size+=A[i]->m*sizeof(T);}
	for(int i=0;i<levels-1;i++){
		if(P[i])size+=P[i]->Memory_Size();
		if(R[i])size+=R[i]->Memory_Size();}
	return size;
}

//////////////////////////////////////////////////////////////////////////
////Aux functions
////MG Cuda (on host) -> MG Cuda (on dev)
template<typename T> void Malloc_And_Copy_Host_To_Device(MultiGridSystemCuda<T>*& ptr_dev,MultiGridSystemCuda<T>*& ptr_host)
{
	ptr_dev=new MultiGridSystemCuda<T>();
	int levels=ptr_host->levels;
	////allocate and copy A, P, R
	SparseMatrixCuda<T>** A_dev=new SparseMatrixCuda<T>*[levels];
	for (int i = 0; i < levels; i++) { A_dev[i] = new SparseMatrixCuda<T>(*ptr_host->A[i], DEVICE); }
	SparseMatrixCuda<T>** P_dev=new SparseMatrixCuda<T>*[levels-1];
	for (int i = 0; i < levels - 1; i++) { P_dev[i] = new SparseMatrixCuda<T>(*ptr_host->P[i], DEVICE);}
	SparseMatrixCuda<T>** R_dev=new SparseMatrixCuda<T>*[levels-1];
	for (int i = 0; i < levels - 1; i++) { R_dev[i] = new SparseMatrixCuda<T>(*ptr_host->R[i], DEVICE);}
	////allocate b, x, r
	T** b_dev=new T*[levels];
	for(int i=0;i<levels;i++){ b_dev[i]=Global_Malloc<T>(ptr_host->A[i]->m,DEVICE);}
	T** x_dev=new T*[levels];
	for(int i=0;i<levels;i++){ x_dev[i]=Global_Malloc<T>(ptr_host->A[i]->m,DEVICE);}
	T** r_dev=new T*[levels];
	for(int i=0;i<levels;i++){ r_dev[i]=Global_Malloc<T>(ptr_host->A[i]->m,DEVICE);}
	////copy b0 and x0
	Global_Copy_Array<T>(b_dev[0], ptr_host->b[0], ptr_host->A[0]->m, DEVICE, HOST);
	Global_Copy_Array<T>(x_dev[0], ptr_host->x[0], ptr_host->A[0]->m, DEVICE, HOST);
	////sync pointers
	ptr_dev->Sync_Pointers_And_Flags(levels,A_dev,P_dev,R_dev,b_dev,x_dev,r_dev,true);
	////coarsest on host
	ptr_dev->coarsest_on_host=ptr_host->coarsest_on_host;
	if (ptr_dev->coarsest_on_host) {
		ptr_dev->A_coarsest_host = new SparseMatrixCuda<T>(*ptr_host->A[levels - 1], HOST);
		int m = ptr_dev->A_coarsest_host->m;
		ptr_dev->b_coarsest_host = Global_Malloc<T>(m, HOST);
		ptr_dev->x_coarsest_host = Global_Malloc<T>(m, HOST);
	}
	////no color update
}

////*****************************************************************************
//// MG CPU (on host) -> MG Cuda (on dev), *** the one that is being used now! ***
////*****************************************************************************
template<typename T, typename T_MG> void Malloc_And_Copy_Host_To_Device(MultiGridSystemCuda<T>*& ptr_dev,T_MG*& ptr_host)
{
	ptr_dev=new MultiGridSystemCuda<T>();
	int levels=ptr_host->levels;
	////allocate and copy A, P, R from host to device
	SparseMatrixCuda<T>** P_dev=new SparseMatrixCuda<T>*[levels-1];
	for (int i = 0; i < levels - 1; i++) { P_dev[i] = new SparseMatrixCuda<T>(*ptr_host->P[i], DEVICE);}
	SparseMatrixCuda<T>** R_dev=new SparseMatrixCuda<T>*[levels-1];
	for (int i = 0; i < levels - 1; i++) { R_dev[i] = new SparseMatrixCuda<T>(*ptr_host->R[i], DEVICE);}

	SparseMatrixCuda<T>** A_dev=new SparseMatrixCuda<T>*[levels];
	if (ptr_host->init_A_levels) {
		for (int i = 0; i < levels; i++) { A_dev[i] = new SparseMatrixCuda<T>(*ptr_host->A[i], DEVICE); }
	}
	else { A_dev[0] = new SparseMatrixCuda<T>(*ptr_host->A[0], DEVICE);}

	////allocate b, x, r on device
	T** b_dev=new T*[levels];
	for(int i=0;i<levels;i++){ b_dev[i]= Global_Malloc<T>(ptr_host->A_size(i),DEVICE);}
	T** x_dev=new T*[levels];
	for(int i=0;i<levels;i++){ x_dev[i]=Global_Malloc<T>(ptr_host->A_size(i), DEVICE);}
	T** r_dev=new T*[levels];
	for(int i=0;i<levels;i++){ r_dev[i]=Global_Malloc<T>(ptr_host->A_size(i), DEVICE);}
	////copy b0 and x0
	T* b0_host=ptr_host->b[0]->data();T* x0_host=ptr_host->x[0]->data();
	Global_Copy_Array<T>(b_dev[0], b0_host, ptr_host->A[0]->rows(), DEVICE, HOST);
	Global_Copy_Array<T>(x_dev[0], x0_host, ptr_host->A[0]->cols(), DEVICE, HOST);
	////sync pointers
	ptr_dev->Sync_Pointers_And_Flags(levels,A_dev,P_dev,R_dev,b_dev,x_dev,r_dev,true);

	////update A[1] to A[levels-1] on device
	if(!ptr_host->init_A_levels){
		Update_A_Levels_On_Device(ptr_dev,ptr_host->one_over_intp);}

	////coarsest on host if ptr_host has A levels
	if(ptr_host->init_A_levels){
		ptr_dev->coarsest_on_host=true;
		ptr_dev->A_coarsest_host = new SparseMatrixCuda<T>(*ptr_host->A[levels - 1], HOST);
		int m=ptr_dev->A_coarsest_host->m;
		ptr_dev->b_coarsest_host = Global_Malloc<T>(m, DEVICE);
		ptr_dev->b_coarsest_host=Global_Malloc<T>(m,HOST);
		ptr_dev->x_coarsest_host=Global_Malloc<T>(m,HOST);}
	else ptr_dev->coarsest_on_host=false;

	////copy color from host to device
	if(ptr_host->color_n.size()>0){
		for(int i=0;i<levels;i++)A_dev[i]->Initialize_Color(ptr_host->color_n[i],ptr_host->block_size,ptr_host->color_ptr[i],ptr_host->color[i]);}
}

////MG CPU (on host) -> MG Cuda (on host)
template<typename T, typename T_MG> void Malloc_And_Copy_Host_To_Host(MultiGridSystemCuda<T>*& ptr_dev,T_MG*& ptr_host)
{
	ptr_dev=new MultiGridSystemCuda<T>();
	int levels=ptr_host->levels;
	////allocate and copy A, P, R
	SparseMatrixCuda<T>** A_dev=new SparseMatrixCuda<T>*[levels];
	for (int i = 0; i < levels; i++) { A_dev[i] = new SparseMatrixCuda<T>(*ptr_host->A[i], HOST); }
	SparseMatrixCuda<T>** P_dev=new SparseMatrixCuda<T>*[levels-1];
	for (int i = 0; i < levels - 1; i++) { P_dev[i] = new SparseMatrixCuda<T>(*ptr_host->P[i], HOST); }
	SparseMatrixCuda<T>** R_dev=new SparseMatrixCuda<T>*[levels-1];
	for (int i = 0; i < levels - 1; i++) { R_dev[i] = new SparseMatrixCuda<T>(*ptr_host->R[i], HOST);}
	////allocate b, x, r
	T** b_dev=new T*[levels];
	for (int i = 0; i < levels; i++) { b_dev[i] = Global_Malloc<T>(ptr_host->A[i]->rows(), HOST);}
	T** x_dev=new T*[levels];
	for (int i = 0; i < levels; i++) { x_dev[i] = Global_Malloc<T>(ptr_host->A[i]->cols(), HOST); }
	T** r_dev=new T*[levels];
	for (int i = 0; i < levels; i++) { r_dev[i] = Global_Malloc<T>(ptr_host->A[i]->rows(), HOST);}
	////copy b0 and x0
	T* b0_host=ptr_host->b[0]->data();T* x0_host=ptr_host->x[0]->data();
	Global_Copy_Array<T>(b_dev[0], b0_host, ptr_host->A[0]->rows(), HOST, HOST);
	Global_Copy_Array<T>(x_dev[0], x0_host, ptr_host->A[0]->cols(), HOST, HOST);
	////sync pointers
	ptr_dev->Sync_Pointers_And_Flags(levels,A_dev,P_dev,R_dev,b_dev,x_dev,r_dev,/*on host*/false);
	////coarsest on host
	ptr_dev->coarsest_on_host=true;
	ptr_dev->A_coarsest_host = new SparseMatrixCuda<T>(*ptr_host->A[levels - 1], HOST);
	int m=ptr_dev->A_coarsest_host->m;
	ptr_dev->b_coarsest_host = Global_Malloc<T>(m, DEVICE);
	ptr_dev->b_coarsest_host = Global_Malloc<T>(m, HOST);
	ptr_dev->x_coarsest_host = Global_Malloc<T>(m, HOST);
	////color
	if(ptr_host->color_n.size()>0){
		for(int i=0;i<levels;i++)A_dev[i]->Initialize_Color(ptr_host->color_n[i],ptr_host->block_size,ptr_host->color_ptr[i],ptr_host->color[i]);}
}

#define Inst_Helper(T) \
template void Malloc_And_Copy_Host_To_Device<T>(MultiGridSystemCuda<T>*&,MultiGridSystemCuda<T>*&); \
template void Malloc_And_Copy_Host_To_Device<T,GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> > >(MultiGridSystemCuda<T>*&,GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> >*&); \
template void Malloc_And_Copy_Host_To_Device<T,GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> > >(MultiGridSystemCuda<T>*&,GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> >*&); \
template void Malloc_And_Copy_Host_To_Host<T,GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> > >(MultiGridSystemCuda<T>*&,GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> >*&); \
template void Malloc_And_Copy_Host_To_Host<T,GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> > >(MultiGridSystemCuda<T>*&,GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> >*&);
Inst_Helper(float);
Inst_Helper(double);
#undef Inst_Helper

//////////////////////////////////////////////////////////////////////////
//// Jacobi

template<class T> __device__ void Jacobi_One_Iteration_Helper(int i,int* ptr,int* col,T* val,T* rhs,T* x,T* x_new)
{
	int first=ptr[i];
	int row_nnz=ptr[i+1]-first;
	T sum=rhs[i];
	T dia=(T)0;
	for(int k=0;k<row_nnz;k++){
		int idx=first+k;
		int j=col[idx];
		T A_ij=val[idx];
		T x_j=x[j];
		if(i==j)dia=A_ij;
		else sum-=A_ij*x_j;}
	sum/=dia;
	x_new[i]=sum;
}

template<class T> __global__ void Jacobi_One_Iteration(int n,int* ptr,int* col,T* val,T* rhs,T* x,T* x_new)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	if(i>=n)return;
	
	//using shared memory
	__shared__ T sval[64][7];

	int first=ptr[i];
	int row_nnz=ptr[i+1]-first;
	
	int si=threadIdx.x;
	for(int k=0;k<row_nnz;k++){sval[si][k]=val[first+k];}	////row_nnz read
	__syncthreads();

	T sum=rhs[i];
	T dia=(T)0;
	for(int k=0;k<row_nnz;k++){
		int j=col[first+k];		////1 read
		if(i==j)dia=sval[si][k];
		else sum-=sval[si][k]*x[j];}
	sum/=dia;
	x_new[i]=sum;	////1 write
}

template<class T> void Jacobi_Device(SparseMatrixCuda<T>*& A_dev,T*& b,T*& x_cur,T*& x_next,T*& r,int iter_num,T tolerance)
{
	int n=A_dev->m;
	int thread_num=64;
	int block_num=(n-1)/thread_num+1;

	for(int i=0;i<iter_num;i++){
		Jacobi_One_Iteration<T><<<block_num,thread_num>>>(n,A_dev->ptr,A_dev->col,A_dev->val,b,x_cur,x_next);
		//r=b-Ax
		Residual(A_dev,x_cur,b,r);
		////r_norm=<r,r>
		T r_norm=Dot(r,r,n);
        if(r_norm<tolerance){std::cout<<"Jacobi converges ("<<r_norm<<" < "<<tolerance<<") after "<<i<<" iterations"<<std::endl;break;}
		////swap current and next pointers
		T* tmp=x_cur;x_cur=x_next;x_next=tmp;}
}

template<class T> void Jacobi(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance)
{
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();
	int n=A_host->m;
	SparseMatrixCuda<T>* A_dev = new SparseMatrixCuda<T>(*A_host, DEVICE);
	T* x[2] = { nullptr,nullptr }; for (int i = 0; i < 2; i++)Global_Realloc_And_Copy_Array<T>(x[i], _x, n, DEVICE, HOST);
	T* x_cur=x[0];T* x_next=x[1];
	T* b = nullptr; Global_Realloc_And_Copy_Array(b, _b, n, DEVICE, HOST);
	T* r = Global_Malloc<T>(n, DEVICE);
	Global_Memset<T>(r, 0, n, DEVICE);

	Jacobi_Device(A_dev,b,x_cur,x_next,r,iter_num,tolerance);
	Global_Copy_Array<T>(_x, x_cur, n, HOST, DEVICE);

	delete A_dev;
	Global_Free(x[0], DEVICE);
	Global_Free(x[1], DEVICE);
	Global_Free(b, DEVICE);
	Global_Free(r, DEVICE);
}

//////////////////////////////////////////////////////////////////////////
//// G-S

template<class T,int thread_n,int nnz_max> __global__ void Multi_Color_Gauss_Seidel_One_Iteration(int n,int* ptr,int* col,T* val,T* rhs,T* x,
	int start,int indir_n,int* indir_idx,int block_size,int order)
{
	int thread_idx=blockIdx.x*blockDim.x+threadIdx.x;
	if(thread_idx>=indir_n)return;
	int block_idx=indir_idx[start+thread_idx];	////1 read

	__shared__ T sval[thread_n][nnz_max];

	for(int dd=0;dd<block_size;dd++){
		int d=(order==0?dd:block_size-1-dd);
		int i=block_idx*block_size+d;
		int first=ptr[i];		////1 read
		int row_nnz=ptr[i+1]-first;	////1 read

		int si=threadIdx.x;
		T sum=rhs[i];
		T dia=(T)0;
		int num_segment=row_nnz/nnz_max+(row_nnz%nnz_max!=0);
		for(int s=0;s<num_segment;s++){
			  int offset=s*nnz_max;
			  int kmax=min(nnz_max,row_nnz-offset);
			  for(int k=0;k<kmax;k++){sval[si][k]=val[first+k+offset];}
			  for(int k=0;k<kmax;k++){
				int j=col[first+k+offset];		////1 read
				if(i==j)dia=sval[si][k];
				else sum-=sval[si][k]*x[j];}}
		sum/=dia;
		x[i]=sum;	////1 write
	}
}

template<class T,int thread_n,int nnz_max> void Multi_Color_Gauss_Seidel_Device(SparseMatrixCuda<T>* A_dev,T* b,T* x,T* r,int iter_num,T tolerance,int order)
{
	int n=A_dev->m;
	for(int i=0;i<iter_num;i++){
		for(int jj=0;jj<A_dev->color_n;jj++){int j=(order==0?jj:A_dev->color_n-1-jj);
			int start=A_dev->color_ptr[j];
			int nc=A_dev->color_ptr[j+1]-start;
			int block_num=(nc-1)/thread_n+1;
			Multi_Color_Gauss_Seidel_One_Iteration<T,thread_n,nnz_max><<<block_num,thread_n>>>
				(n,A_dev->ptr,A_dev->col,A_dev->val,b,x,start,nc,A_dev->color,A_dev->block_size,order);}}
}

template<class T,int thread_n,int nnz_max> void Multi_Color_Gauss_Seidel(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,int order)
{
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();

	int n=A_host->m;
	SparseMatrixCuda<T>* A_dev = new SparseMatrixCuda<T>(*A_host, DEVICE);
	T* b=nullptr; Global_Realloc_And_Copy_Array<T>(b,_b,n,DEVICE,HOST);
	T* x=nullptr; Global_Realloc_And_Copy_Array<T>(x,_x,n,DEVICE,HOST);
	T* r = nullptr; r = Global_Malloc<T>(n, DEVICE); Global_Memset<T>(r, 0, n, DEVICE);

	Multi_Color_Gauss_Seidel_Device<T,thread_n,nnz_max>(A_dev,b,x,r,iter_num,tolerance,order);
	Global_Copy_Array(_x,x,n,HOST,DEVICE);
    
	delete A_dev;
	Global_Free(b, DEVICE);
	Global_Free(x, DEVICE);
	Global_Free(r, DEVICE);
}

////Serialized Gauss-Seidel
template<class T> __global__ void Serialized_Gauss_Seidel_One_Iteration(int n,int* ptr,int* col,T* val,T* rhs,T* x,int order)
{
	for(int ii=0;ii<n;ii++){int i=(order==0?ii:(n-1-ii));
		int first=ptr[i];
		int row_nnz=ptr[i+1]-first;
		T sum=rhs[i];
		T dia=(T)0;
		for(int k=0;k<row_nnz;k++){
			int idx=first+k;
			int j=col[idx];
			T A_ij=val[idx];
			T x_j=x[j];
			if(i==j)dia=A_ij;
			else sum-=A_ij*x_j;}
		sum/=dia;
		x[i]=sum;}	
}

template<class T> void Serialized_Gauss_Seidel_Device(SparseMatrixCuda<T>* A_dev,T* b,T* x,T* r,int iter_num,T tolerance,int order/*=0*/)
{
	int n=A_dev->m;
	for(int i=0;i<iter_num;i++){
		Serialized_Gauss_Seidel_One_Iteration<T><<<1,1>>>(n,A_dev->ptr,A_dev->col,A_dev->val,b,x,order);
		////r=b-Ax
		Residual(A_dev,x,b,r);
		////r_norm=<r,r>
		T r_norm=Dot(r,r,n);
        if(r_norm<tolerance){break;}}
}

template<class T> void Serialized_Gauss_Seidel(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,int order/*=0*/)
{
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();

	int n=A_host->m;
	SparseMatrixCuda<T>* A_dev = new SparseMatrixCuda<T>(*A_host, DEVICE);
	T* b = nullptr; Global_Realloc_And_Copy_Array<T>(b, _b, n, DEVICE, HOST);
	T* x = nullptr; Global_Realloc_And_Copy_Array<T>(x, _x, n, DEVICE, HOST);
	T* r = Global_Malloc<T>(n, DEVICE); Global_Memset<T>(r, 0, n, DEVICE);

	Serialized_Gauss_Seidel_Device(A_dev,b,x,r,iter_num,tolerance,order);
	Global_Copy_Array<T>(_x, x, n, HOST, DEVICE);
    
	delete A_dev;
	Global_Free(b, DEVICE);
	Global_Free(x, DEVICE);
	Global_Free(r, DEVICE);
}

//////////////////////////////////////////////////////////////////////////
//// diagonal preconditioner

template<class T> __global__ void Build_Diagonal_Preconditioner(int n,int* ptr,int* col,T* val,T* _dia)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;
	if(i>=n)return;
	
	int first=ptr[i];
	int row_nnz=ptr[i+1]-first;
	T dia=(T)0;
	for(int k=0;k<row_nnz;k++){
		int j=col[first+k];		////1 read
		if(i==j){dia=val[first+k];break;}}	////1 read
	if(dia!=(T)0)_dia[i]=(T)1/dia;		//// 1 write
	else _dia[i]=(T)1;
}

template<class T> void Build_Diagonal_Preconditioner(SparseMatrixCuda<T>*& A_dev,T*& pd)
{
	int n=A_dev->m;
	const int thread_num_per_block=256;
	int block_num=(n-1)/thread_num_per_block+1;

	Build_Diagonal_Preconditioner<T><<<block_num,thread_num_per_block>>>(n,A_dev->ptr,A_dev->col,A_dev->val,pd);
}

template<class T> __global__ void Apply_Diagonal_Preconditioner_Device(int n,T* dia,T* r,T* x)
{
	int i=blockIdx.x*blockDim.x+threadIdx.x;if(i>=n)return;
	x[i]=r[i]*dia[i];
}

template<class T> void Apply_Diagonal_Preconditioner(int n,T*& dia,T* r,T* x)
{
	const int thread_num_per_block=256;int block_num=(n-1)/thread_num_per_block+1;

	Apply_Diagonal_Preconditioner_Device<T><<<block_num,thread_num_per_block>>>(n,dia,r,x);
}

//////////////////////////////////////////////////////////////////////////
////Preconditioned Conjugate Gradient
template<class MtxT,class T> void Preconditioned_Conjugate_Gradient_Device(MtxT* A,T* b,T* x,T* r,T* d,T* q,T* s,T* pd,int iter_num,T tolerance,bool use_precond/*=true*/)
{
	if(use_precond)Build_Diagonal_Preconditioner<T>(A,pd);

	bool verbose=true;
	T one=(T)1;
	const int n=A->m;
	const int restart_iter_num=50;

	Residual(A,x,b,r);													//// r=b-Ax
	if(use_precond)Apply_Diagonal_Preconditioner<T>(n,pd,r,d);			//// d=precond(r)
	else Copy(d,r,n);
   
	T delta_0=Dot(r,d,n);												//// delta=<r^H,z>
	T delta=delta_0;
	T epsilon_0=tolerance*tolerance*delta_0;
    int iter=0;
	while(delta>epsilon_0&&iter<iter_num){
		Mv(A,d,q);														//// q = Ad
        T alpha=delta/Dot(q,d,n);		
		Axpy(&alpha,d,x,n);												//// x = x+a*d
		if(iter%restart_iter_num==0)Residual(A,x,b,r);
		else{T neg_alpha=-alpha;Axpy(&neg_alpha,q,r,n);}				//// r = r-a*q

		if(use_precond)Apply_Diagonal_Preconditioner<T>(n,pd,r,s);		//// s = precond(r)
		else Copy(s,r,n);	
		
		T delta_old=delta;
		delta=Dot(r,s,n);
        T beta=delta/delta_old;
		Scale(&beta,d,n);				
		Axpy(&one,s,d,n);												//// d = s+beta*d;
		cudaThreadSynchronize();iter++;}

	if(verbose)std::cout<<"cuda pcg converges in "<<iter<<" iterations"<<std::endl;
}

template<class MtxT,class T> void Preconditioned_Conjugate_Gradient(MtxT* A_host,T* _b,T* _x,int iter_num,T tolerance,bool use_precond)
{
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();

	int n=A_host->rows();
	SparseMatrixCuda<T>* A = new SparseMatrixCuda<T>(*A_host, DEVICE);
	T* b = nullptr; Global_Realloc_And_Copy_Array(b, _b, n, DEVICE, HOST);
	T* x = nullptr; Global_Realloc_And_Copy_Array(x, _x, n, DEVICE, HOST);
	T* r = Global_Malloc<T>(n, DEVICE); Global_Memset<T>(r, 0, n, DEVICE);
	T* d = Global_Malloc<T>(n, DEVICE);
	T* q = Global_Malloc<T>(n, DEVICE);
	T* s = Global_Malloc<T>(n, DEVICE);
	T* pd = Global_Malloc<T>(n, DEVICE);

	Preconditioned_Conjugate_Gradient_Device(A,b,x,r,d,q,s,pd,iter_num,tolerance,use_precond);

	Global_Copy_Array<T>(_x, x, n, HOST, DEVICE);
	
	delete A;

	Global_Free(b, DEVICE);
	Global_Free(x, DEVICE);
	Global_Free(r, DEVICE);
	Global_Free(d, DEVICE);
	Global_Free(q, DEVICE);
	Global_Free(s, DEVICE);
	Global_Free(pd, DEVICE);
}

//////////////////////////////////////////////////////////////////////////
//// direct solver

int Direct_Solve_Device(SparseMatrixCuda<double>* A_dev,double* b,double* x,double tolerance)
{
	int singularity=0;
	checkCudaErrors(cusolverSpDcsrlsvchol(Cusolver_Handle(),A_dev->m,A_dev->nnz,
		Cusparse_Mat_Descr(),A_dev->val,A_dev->ptr,A_dev->col,b,tolerance,/*reorder*/0,x,&singularity)); 
	return singularity;
}

int Direct_Solve_Device(SparseMatrixCuda<float>* A_dev,float* b,float* x,float tolerance)
{
	int singularity=0;checkCudaErrors(cusolverSpScsrlsvchol(Cusolver_Handle(),A_dev->m,A_dev->nnz,
		Cusparse_Mat_Descr(),A_dev->val,A_dev->ptr,A_dev->col,b,tolerance,/*reorder*/0,x,&singularity)); 
	return singularity;
}

int Direct_Solve_Host(SparseMatrixCuda<double>* A_dev,double* b,double* x,double tolerance)
{
	int singularity=0;checkCudaErrors(cusolverSpDcsrlsvcholHost(Cusolver_Handle(),A_dev->m,A_dev->nnz,
		Cusparse_Mat_Descr(),A_dev->val,A_dev->ptr,A_dev->col,b,tolerance,/*reorder*/0,x,&singularity)); 
	return singularity;
}

int Direct_Solve_Host(SparseMatrixCuda<float>* A_dev,float* b,float* x,float tolerance)
{
	int singularity=0;checkCudaErrors(cusolverSpScsrlsvcholHost(Cusolver_Handle(),A_dev->m,A_dev->nnz,
		Cusparse_Mat_Descr(),A_dev->val,A_dev->ptr,A_dev->col,b,tolerance,/*reorder*/0,x,&singularity)); 
	return singularity;
}

//////////////////////////////////////////////////////////////////////////
//// multi-grid

//// GPU multigrid v-cycle *** The one that is being used ***
template<class T> void V_Cycle_Device(MultiGridSystemCuda<T>* mg_dev,int iter_num,T tolerance)
{
	T one=(T)1;
	int levels=mg_dev->levels;
	SparseMatrixCuda<T>** A_dev=mg_dev->A;
	SparseMatrixCuda<T>** P_dev=mg_dev->P;
	SparseMatrixCuda<T>** R_dev=mg_dev->R;
	T** b_dev=mg_dev->b;
	T** x_dev=mg_dev->x;
	T** r_dev=mg_dev->r;

	for(int i=0;i<=levels-2;i++){
		////smooth A_i*x_i=b_i
		if (i > 0) Global_Memset<T>(x_dev[i], 0, A_dev[i]->m, DEVICE);
		Multi_Color_Gauss_Seidel_Device<T,32,32>(A_dev[i],b_dev[i],x_dev[i],r_dev[i],iter_num,tolerance,0);
		//Serialized_Gauss_Seidel_Device(A_dev[i],b_dev[i],x_dev[i],r_dev[i],iter_num,tolerance,0);
		////r_i=b_i-A_i*x_i
		Residual(A_dev[i],x_dev[i],b_dev[i],r_dev[i]);
		////b_{i+1}=R_i*r_i
		Mv(R_dev[i],r_dev[i],b_dev[i+1]);}

	////solve A_{n-1}x_{n-1}=b_{n-1}
	if (mg_dev->coarsest_on_host) {	////direct solve for the coarsest level
		int m = mg_dev->A[levels - 1]->m;
		Global_Copy_Array<T>(mg_dev->b_coarsest_host, b_dev[levels - 1], m, HOST, DEVICE);
		Global_Copy_Array<T>(mg_dev->x_coarsest_host, x_dev[levels - 1], m, HOST, DEVICE);
		Direct_Solve_Host(mg_dev->A_coarsest_host, mg_dev->b_coarsest_host, mg_dev->x_coarsest_host, tolerance);
		Global_Copy_Array<T>(x_dev[levels - 1], mg_dev->x_coarsest_host, m, DEVICE, HOST);
	}
	else{Direct_Solve_Device(A_dev[levels-1],b_dev[levels-1],x_dev[levels-1],tolerance);}

	for(int i=levels-2;i>=0;i--){
		////x_i=x_i+P_i*x_{i+1}
		Csrmv(P_dev[i],x_dev[i+1],x_dev[i],&one,&one);
		////smooth A_i*x_i=b_i
		Multi_Color_Gauss_Seidel_Device<T,32,32>(A_dev[i],b_dev[i],x_dev[i],r_dev[i],1,tolerance,1);
		//Serialized_Gauss_Seidel_Device(A_dev[i],b_dev[i],x_dev[i],r_dev[i],iter_num,tolerance,1);
	}
}

template<class T> void Multi_Grid(MultiGridSystemCuda<T>* mg_host,int iter_num,T tolerance)
{
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();
	MultiGridSystemCuda<T>* mg_dev=nullptr;
	Malloc_And_Copy_Host_To_Device(mg_dev,mg_host);
	V_Cycle_Device(mg_dev,iter_num,tolerance);
	Global_Copy_Array<T>(mg_host->x[0], mg_dev->x[0], mg_dev->A[0]->m, HOST, DEVICE);
	delete mg_dev;
}

template<class T> void Update_A_Levels_On_Device(MultiGridSystemCuda<T>* mg_dev,T scale,bool resize/*=false*/)
{
	int levels=mg_dev->levels;
	bool realloc=resize||(mg_dev->AP==0);
	if (mg_dev->AP == 0) {
		mg_dev->AP = new SparseMatrixCuda<T> * [levels - 1];}
	if (realloc) {
		for (int i = 0; i < levels - 1; i++) {
			mg_dev->AP[i] = new SparseMatrixCuda<T>(UNKNOWN);
			mg_dev->A[i + 1] = new SparseMatrixCuda<T>(UNKNOWN);}}
	for (int i = 0; i < levels - 1; i++) {
		//std::cout<<"A "<<i<<": "<<mg_dev->A[i]->nnz<<std::endl;
		Mm(mg_dev->A[i], mg_dev->P[i], mg_dev->AP[i]);
		Mm(mg_dev->R[i], mg_dev->AP[i], mg_dev->A[i + 1]);
		Scale(&scale, mg_dev->A[i + 1]->val, mg_dev->A[i + 1]->nnz);}
}

////Multigrid Preconditioner Conjugate Gradient
template<class T> void GMGPCG_Device(MultiGridSystemCuda<T>* mg_dev,T* b,T* x,T* p,T* q,int iter_num,T tolerance,bool verbose)
{
	T one=(T)1;
	SparseMatrixCuda<T>* A=mg_dev->A[0];
	int n=A->m;
	T* s=mg_dev->x[0];
	T* r=mg_dev->b[0];

    ////r=b-Ax
	Residual(A,x,b,r);
	T b_norm=sqrt(Dot(b,b,n));
	
	T eps=b_norm*tolerance;
	T eps2=eps*eps;
	T rho1=(T)2*eps2;
	T rho2=(T)0;
	T res_norm=sqrt(Dot(r,r,n));
	T res2=res_norm*res_norm;
	
	if(rho1==(T)0){
		if(verbose)std::cout<<"cuda mgpcg stops in 0 iterations "<<std::endl;
		return;}

	for(int iter=0;iter<iter_num;iter++){
		Global_Memset<T>(s, 0, n, DEVICE);
		V_Cycle_Device(mg_dev,1,tolerance);

		rho2=rho1;
		rho1=Dot(r,s,n);

		if(iter){
			T a=rho1/rho2;
			Scale(&a,p,n);	
			Axpy(&one,s,p,n);}
		else Copy(p,s,n);
		Mv(A,p,q);
		T qp=Dot(q,p,n);
		T alpha=rho1/qp;
		Axpy(&alpha,p,x,n);

		if(iter%30!=0){
			T neg_alpha=-alpha;
			Axpy(&neg_alpha,q,r,n);}
		else{Residual(A,x,b,r);}

		res2=Dot(r,r,n);

		if(res2<=eps2){
			if(verbose)std::cout<<"cuda mgpcg converges in "<<iter+1<<" iterations "<<std::endl;
			return;}}

	if(verbose)std::cout<<"cuda mgpcg does not converge in "<<iter_num<<" iterations "<<std::endl;
}

template<class T,class T_MG> void GMGPCG(T_MG* mg_host,T* b_host,T* x_host,int iter_num,T tolerance,bool verbose)
{
	Timer<real> timer;timer.Reset();
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();
	MultiGridSystemCuda<T>* mg_dev=nullptr;
	Malloc_And_Copy_Host_To_Device(mg_dev,mg_host);
	if(verbose){AuxFunc::Seperation_Line();timer.Elapse_And_Output_And_Reset("GPU GMGPCG: cpu to gpu transfer");}

	int n=mg_dev->A[0]->m;
	T* b = nullptr; Global_Realloc_And_Copy_Array(b, b_host, n, DEVICE, HOST);
	T* x = nullptr; Global_Realloc_And_Copy_Array(x, x_host, n, DEVICE, HOST);
	T* p = Global_Malloc<T>(n, DEVICE);
	T* q = Global_Malloc<T>(n, DEVICE);

	GMGPCG_Device(mg_dev,b,x,p,q,iter_num,tolerance,true);	////always turn on pcg verbose to show the convergence iterations
	if(verbose)timer.Elapse_And_Output_And_Reset("GPU GMGPCG: gpu solving");

	Global_Copy_Array<T>(x_host, x, n, HOST, DEVICE);
	
	delete mg_dev;
	Global_Free(b, DEVICE);
	Global_Free(x, DEVICE);
	Global_Free(p, DEVICE);
	Global_Free(q, DEVICE);
	if(verbose){timer.Elapse_And_Output_And_Reset("GPU GMGPCG: gpu to cpu transfer");AuxFunc::Seperation_Line();}
}

template<class T,class T_MG> void MICPCG(T_MG* ptr_host,T* b_host,T* x_host,int iter_num,T tolerance,bool verbose/*=false*/)
{
	Timer<real> timer;timer.Reset();
	if(!Is_Cuda_Context_Initialized())Initialize_Cuda_Context();
	MultiGridSystemCuda<T>* ptr_dev=nullptr;

	//////////////////////////////////////////////////////////////////////////
	////Malloc_And_Copy_Host_To_Device(mg_dev,mg_host);
	//////////////////////////////////////////////////////////////////////////
	ptr_dev=new MultiGridSystemCuda<T>();
	int levels=1;
	////allocate and copy A

	SparseMatrixCuda<T>** A_dev=new SparseMatrixCuda<T>*[levels];
	A_dev[0] = new SparseMatrixCuda<T>(*ptr_host->A[0], DEVICE);

	////allocate b, x, r on device
	T** b_dev=new T*[levels];
	for(int i=0;i<levels;i++){ b_dev[i]= Global_Malloc<T>(ptr_host->A_size(i),DEVICE);}
	T** x_dev=new T*[levels];
	for(int i=0;i<levels;i++){ x_dev[i]= Global_Malloc<T>(ptr_host->A_size(i), DEVICE);}
	T** r_dev=new T*[levels];
	for(int i=0;i<levels;i++){ r_dev[i]= Global_Malloc<T>(ptr_host->A_size(i), DEVICE);}

	////copy b0 and x0
	T* b0_host=ptr_host->b[0]->data();T* x0_host=ptr_host->x[0]->data();
	Global_Copy_Array<T>(b_dev[0], b0_host, ptr_host->A[0]->rows(), DEVICE, HOST);
	Global_Copy_Array<T>(x_dev[0], x0_host, ptr_host->A[0]->cols(), DEVICE, HOST);
	//////////////////////////////////////////////////////////////////////////

	////sync pointers
	ptr_dev->levels=1;
	ptr_dev->A=A_dev;
	ptr_dev->b=b_dev;
	ptr_dev->x=x_dev;
	ptr_dev->r=r_dev;
	ptr_dev->on_device=true;

	if(verbose){AuxFunc::Seperation_Line();timer.Elapse_And_Output_And_Reset("GPU GMGPCG: cpu to gpu transfer");}

	int n=ptr_dev->A[0]->m;
	T* b = nullptr; Global_Realloc_And_Copy_Array(b, b_host, n, DEVICE, HOST);
	T* x = nullptr; Global_Realloc_And_Copy_Array(x, x_host, n, DEVICE, HOST);
	T* p = Global_Malloc<T>(n, DEVICE);
	T* q = Global_Malloc<T>(n, DEVICE);

	MICPCG_Device(ptr_dev,b,x,p,q,iter_num,tolerance,true);	////always turn on pcg verbose to show the convergence iterations
	if(verbose)timer.Elapse_And_Output_And_Reset("GPU GMGPCG: gpu solving");

	Global_Copy_Array<T>(x_host, x, n, HOST, DEVICE);
	
	delete ptr_dev;
	Global_Free(b, DEVICE);
	Global_Free(x, DEVICE);
	Global_Free(p, DEVICE);
	Global_Free(q, DEVICE);
	if(verbose){timer.Elapse_And_Output_And_Reset("GPU GMGPCG: gpu to cpu transfer");AuxFunc::Seperation_Line();}
}

template<class T> void MICPCG_Device(MultiGridSystemCuda<T>* mg_dev,T* b,T* x,T* p,T* q,int iter_num,T tolerance,bool verbose/*=true*/)
{
	T one=(T)1;
	SparseMatrixCuda<T>* A=mg_dev->A[0];
	int n=A->m;
	T* s=mg_dev->x[0];
	T* r=mg_dev->b[0];

    ////r=b-Ax
	Residual(A,x,b,r);
	T b_norm=sqrt(Dot(b,b,n));
	
	T eps=b_norm*tolerance;
	T eps2=eps*eps;
	T rho1=(T)2*eps2;
	T rho2=(T)0;
	T res_norm=sqrt(Dot(r,r,n));
	T res2=res_norm*res_norm;
	
	if(rho1==(T)0){
		if(verbose)std::cout<<"cuda micpcg stops in 0 iterations "<<std::endl;
		return;}

	for(int iter=0;iter<iter_num;iter++){
		Global_Memset<T>(s, 0, n, DEVICE);
		////TODO: MIC preconditioner
		Copy(s,r,n);

		rho2=rho1;
		rho1=Dot(r,s,n);

		if(iter){
			T a=rho1/rho2;
			Scale(&a,p,n);	
			Axpy(&one,s,p,n);}
		else Copy(p,s,n);
		Mv(A,p,q);
		T qp=Dot(q,p,n);
		T alpha=rho1/qp;
		Axpy(&alpha,p,x,n);

		if(iter%30!=0){
			T neg_alpha=-alpha;
			Axpy(&neg_alpha,q,r,n);}
		else{Residual(A,x,b,r);}

		res2=Dot(r,r,n);

		if(res2<=eps2){
			if(verbose)std::cout<<"cuda micpcg converges in "<<iter+1<<" iterations "<<std::endl;
			return;}}

	if(verbose)std::cout<<"cuda micpcg does not converge in "<<iter_num<<" iterations "<<std::endl;
}

#define Inst_Helper(T)	\
template void Jacobi<T>(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance);	\
template void Jacobi_Device<T>(SparseMatrixCuda<T>*& A_dev, T*& b, T*& x_cur, T*& x_next, T*& r, int iter_num, T tolerance); \
template void Serialized_Gauss_Seidel<T>(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,int order); \
template void Serialized_Gauss_Seidel_Device<T>(SparseMatrixCuda<T>* A_dev,T* b,T* x,T* r,int iter_num,T tolerance,int order); \
template void Preconditioned_Conjugate_Gradient<SparseMatrixCuda<T>,T>(SparseMatrixCuda<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,bool use_precond); \
template void Preconditioned_Conjugate_Gradient<SparseMatrix<T>,T>(SparseMatrix<T>* A_host,T* _b,T* _x,int iter_num,T tolerance,bool use_precond); \
template void Multi_Grid<T>(MultiGridSystemCuda<T>* mg_host,int iter_num,T tolerance); \
template void V_Cycle_Device<T>(MultiGridSystemCuda<T>* mg_dev,int iter_num,T tolerance); \
template void GMGPCG<T,MultiGridSystemCuda<T> >(MultiGridSystemCuda<T>* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose); \
template void GMGPCG<T,GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> > >(GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> >* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose); \
template void GMGPCG<T,GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> > >(GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> >* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose); \
template void GMGPCG_Device<T>(MultiGridSystemCuda<T>* mg_dev,T* b,T* x,T* p,T* y,int iter_num,T tolerance,bool verbose); \
template void MICPCG<T,GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> > >(GeometricMultiGridSystem<2,T,SparseMatrix<T>,VectorN<T> >* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose); \
template void MICPCG<T,GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> > >(GeometricMultiGridSystem<3,T,SparseMatrix<T>,VectorN<T> >* mg_host,T* _b,T* _x,int iter_num,T tolerance,bool verbose); \
template void MICPCG_Device<T>(MultiGridSystemCuda<T>* mg_dev,T* b,T* x,T* p,T* y,int iter_num,T tolerance,bool verbose); 
Inst_Helper(double);
Inst_Helper(float);
#undef Inst_Helper

#define Inst_Helper_Gauss_Seidel(thread_n,nnz_max) \
template void Multi_Color_Gauss_Seidel<double,thread_n,nnz_max>(SparseMatrixCuda<double>* A_host,double* _b,double* _x,int iter_num,double tolerance,int order); \
template void Multi_Color_Gauss_Seidel_Device<double,thread_n,nnz_max>(SparseMatrixCuda<double>* A_dev,double* b,double* x,double* r,int iter_num,double tolerance,int order); \
template void Multi_Color_Gauss_Seidel<float,thread_n,nnz_max>(SparseMatrixCuda<float>* A_host,float* _b,float* _x,int iter_num,float tolerance,int order); \
template void Multi_Color_Gauss_Seidel_Device<float,thread_n,nnz_max>(SparseMatrixCuda<float>* A_dev,float* b,float* x,float* r,int iter_num,float tolerance,int order);
Inst_Helper_Gauss_Seidel(16,27); ////Poisson 3D
Inst_Helper_Gauss_Seidel(16,9); ////Poisson 2D
Inst_Helper_Gauss_Seidel(16,81); ////FEM 3D
Inst_Helper_Gauss_Seidel(16,18); ////FEM 2D
Inst_Helper_Gauss_Seidel(32,27); ////Poisson 3D
Inst_Helper_Gauss_Seidel(32,9); ////Poisson 2D
Inst_Helper_Gauss_Seidel(32,81); ////FEM 3D
Inst_Helper_Gauss_Seidel(32,18); ////FEM 2D
Inst_Helper_Gauss_Seidel(64,27); ////Poisson 3D
Inst_Helper_Gauss_Seidel(64,9); ////Poisson 2D
Inst_Helper_Gauss_Seidel(64,81); ////FEM 3D
Inst_Helper_Gauss_Seidel(64,18); ////FEM 2D
#undef Inst_Helper_Gauss_Seidel



};

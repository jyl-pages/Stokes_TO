//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#include "ContextCuda.h"
#include "AuxFuncCuda.h"

namespace AuxFuncCuda
{
	
static double one_d=1.;
static double zero_d=0.;
static double neg_one_d=-1.;
static float one_f=1.;
static float zero_f=0.;
static float neg_one_f=-1.;
cublasStatus_t cublas_status;

//////////////////////////////////////////////////////////////////////////
////Basic APIs

////Data holder calculation
DataHolder Copy_Data_Holder(const DataHolder& to, const DataHolder& from) {
	if (from == UNKNOWN) {
		std::cerr << "[Error] AuxFuncCuda::Copy_Data_Holder: from=UNKNOWN\n";
		assert(false);
	}
	if (to == UNKNOWN) return from;
	else return to;
}

////Memory allocation, communication between CPU and GPU

template<class T> T* Global_Malloc(int n, const DataHolder &side) {
	T* p = nullptr;
	if (side == UNKNOWN) { std::cerr << "[Error] AuxFuncCuda::Global_Malloc: unknown data holder\n"; }
	else if (side == HOST) { p = new T[n]; }
	else { checkCudaErrors(cudaMalloc((void**)&p, n * sizeof(T))); }
	return p;
}

template<class T> T* Global_Free(T*& ptr, const DataHolder& side) {
	if (ptr == nullptr) return nullptr;
	if (side == UNKNOWN) { std::cerr << "[Error] AuxFuncCuda::Global_Malloc: unknown data holder\n"; }
	else if (side == HOST) { delete[] ptr; }
	else if (side == DEVICE) { checkCudaErrors(cudaFree((void*)ptr)); }
	return nullptr;
}

template<class T> void Global_Copy_Array(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side) {
	size_t bytenum = sizeof(T) * n;
	cudaMemcpyKind op_kind;
	if (from_side == HOST && to_side == HOST) { op_kind = cudaMemcpyHostToHost; }
	else if (from_side == HOST && to_side == DEVICE) { op_kind = cudaMemcpyHostToDevice; }
	else if (from_side == DEVICE && to_side == HOST) { op_kind = cudaMemcpyDeviceToHost; }
	else if (from_side == DEVICE && to_side == DEVICE) { op_kind = cudaMemcpyDeviceToDevice; }
	else { std::cerr << "[Error] AuxFuncCuda::Global_Copy_Array: unknown data holder\n"; }
	cudaMemcpy(ptr_to, ptr_from, bytenum, op_kind);
}

template<class T> T* Global_Malloc_And_Copy_Array(T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side) {
	T* ptr_to = Global_Malloc<T>(n, to_side);
	Global_Copy_Array<T>(ptr_to, ptr_from, n, to_side, from_side);
	return ptr_to;
}
template<class T> void Global_Realloc_And_Copy_Array(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side) {
	Global_Free<T>(ptr_to, to_side);
	ptr_to = Global_Malloc_And_Copy_Array<T>(ptr_from, n, to_side, from_side);
}

//only 1 element
template<class T> void Global_Realloc_And_Copy(T*& ptr_to, T*& ptr_from, const DataHolder& to_side, const DataHolder& from_side)
{
	Global_Realloc_And_Copy_Array(ptr_to, ptr_from, 1, to_side, from_side);
}

template<class T> void Global_Memset(T* ptr, int v, int n, const DataHolder& side) {
	size_t bytenum = n * sizeof(T);
	if (side == HOST) { memset(ptr, v, bytenum); }
	else if (side == DEVICE) { cudaMemset(ptr, v, bytenum); }
	else { std::cerr << "[Error] AuxFuncCuda::Global_Memset: unknown data type\n"; }
}

/*template<class T> void Malloc_On_Device(T*& ptr_dev,int n)
{ptr_dev = Global_Malloc<T>(n, DEVICE);}

template<class T> void Free_On_Device(T*& ptr_dev)
{Global_Free<T>(ptr_dev, DEVICE);}

template<class T> void Malloc_On_Host(T*& ptr_host, int n)
{ptr_host = Global_Malloc<T>(n, HOST);}

template<class T> void Free_On_Host(T*& ptr_host)
{Global_Free<T>(ptr_host, HOST);}

template<class T> void Copy_Array_Host_To_Device(T*& ptr_dev,T*& ptr_host,int n)
{Global_Copy_Array<T>(ptr_dev, ptr_host, n, DEVICE, HOST);}
template<class T> void Copy_Array_Device_To_Host(T*& ptr_host, T*& ptr_dev, int n)
{Global_Copy_Array<T>(ptr_host, ptr_dev, n, HOST, DEVICE);}
template<class T> void Copy_Array_Device_To_Device(T*& ptr_dev_to, T*& ptr_dev_from, int n)
{Global_Copy_Array<T>(ptr_dev_to, ptr_dev_from, n, DEVICE, DEVICE);}
template<class T> void Copy_Array_Host_To_Host(T*& ptr_host_to,T*& ptr_host_from,int n)
{Global_Copy_Array<T>(ptr_host_to, ptr_host_from, n, HOST, HOST);}


template<class T> void Malloc_And_Copy_Array_Host_To_Device(T*& ptr_dev, T*& ptr_host, int n)
{Global_Realloc_And_Copy_Array<T>(ptr_dev, ptr_host, n, DEVICE, HOST);}
template<class T> void Malloc_And_Copy_Array_Host_To_Host(T*& ptr_host_to,T*& ptr_host_from,int n)
{Global_Realloc_And_Copy_Array<T>(ptr_host_to, ptr_host_from, n, HOST, HOST);}
template<class T> void Malloc_And_Copy_Array_Device_To_Device(T*& ptr_dev_to,T*& ptr_dev_from,int n)
{Global_Realloc_And_Copy_Array<T>(ptr_dev_to, ptr_dev_from, n, DEVICE, DEVICE);}
template<class T> void Malloc_And_Copy_Array_Device_To_Host(T*& ptr_host_to, T*& ptr_dev_from, int n)
{Global_Realloc_And_Copy_Array<T>(ptr_host_to, ptr_dev_from, n, HOST, DEVICE);}



template<class T> void Malloc_And_Copy_Host_To_Device(T*& ptr_dev, T*& ptr_host)
{Global_Realloc_And_Copy<T>(ptr_dev, ptr_host, DEVICE, HOST);}
template<class T> void Malloc_And_Copy_Host_To_Host(T*& ptr_host_to,T*& ptr_host_from)
{Global_Realloc_And_Copy<T>(ptr_host_to, ptr_host_from, HOST, HOST);}
template<class T> void Malloc_And_Copy_Device_To_Device(T*& ptr_dev_to, T*& ptr_dev_from)
{Global_Realloc_And_Copy<T>(ptr_dev_to, ptr_dev_from, DEVICE, DEVICE);}

template<class T> void Memset_On_Device(T* x, int v, int n)
{Global_Memset<T>(x, v, n, DEVICE);}*/

#define Inst_Helper(T) \
template T* Global_Malloc<T>(int n, const DataHolder &side); \
template T* Global_Free<T>(T*& ptr, const DataHolder &side); \
template void Global_Copy_Array<T>(T*& ptr_to, T*& ptr_from, int n, const DataHolder &to_side, const DataHolder &from_side); \
template T* Global_Malloc_And_Copy_Array<T>(T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side); \
template void Global_Realloc_And_Copy_Array<T>(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side) ;\
template void Global_Realloc_And_Copy<T>(T*& ptr_to, T*& ptr_from, const DataHolder& to_side, const DataHolder& from_side);\
template void Global_Memset<T>(T* ptr, int v, int n, const DataHolder& side);/*\
template void Malloc_On_Device<T>(T*&,int); \
template void Free_On_Device<T>(T*&); \
template void Malloc_On_Host<T>(T*&,int); \
template void Free_On_Host<T>(T*&); \
template void Copy_Array_Host_To_Device<T>(T*&,T*&,int); \
template void Copy_Array_Device_To_Host<T>(T*&,T*&,int); \
template void Copy_Array_Device_To_Device<T>(T*&,T*&,int); \
template void Copy_Array_Host_To_Host<T>(T*&,T*&,int); \
template void Malloc_And_Copy_Array_Host_To_Device<T>(T*&,T*&,int); \
template void Malloc_And_Copy_Array_Host_To_Host<T>(T*&,T*&,int); \
template void Malloc_And_Copy_Array_Device_To_Device<T>(T*&,T*&,int); \
template void Malloc_And_Copy_Array_Device_To_Host<T>(T*& ptr_host_to,T*& ptr_dev_from,int n);\
template void Malloc_And_Copy_Host_To_Device<T>(T*&,T*&); \
template void Malloc_And_Copy_Host_To_Host<T>(T*&,T*&); \
template void Malloc_And_Copy_Device_To_Device<T>(T*&,T*&); \
template void Memset_On_Device<T>(T*,int,int);*/
Inst_Helper(int);
Inst_Helper(float);
Inst_Helper(double);
#undef Inst_Helper

template<class T>
cusparseDnVecDescr_t Create_DnVecDescr_t(T* x, int size) {//x must be on device
	cusparseDnVecDescr_t vec_t = nullptr;
	cusparseStatus_t stat = cusparseCreateDnVec(&vec_t, size, x, Cuda_Real_Type<T>());
	return vec_t;
}

////CUDA11 TOFIX [JY]: csrmv
////Linear algebra operations on device
template<class T>//T={double, float}
void Csrmv(SparseMatrixCuda<T>* A, T* x, T* y, T* alpha, T* beta)//y=alpha*A*x+beta*y
{
	//A,x,y must be on device
	//This function will internally allocate and erase temporary space dBuffer
	cusparseHandle_t handle = nullptr; cusparseCreate(&handle);
	cusparseSpMatDescr_t A_desc = A->Get_SpMatCescr_t();
	cusparseDnVecDescr_t x_desc = Create_DnVecDescr_t(x, A->n);
	cusparseDnVecDescr_t y_desc = Create_DnVecDescr_t(y, A->m);
	size_t buffersize; void* dBuffer = nullptr;
	cusparseSpMV_bufferSize(
		handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		alpha, A_desc, x_desc, beta, y_desc, A->data_type,
		CUSPARSE_MV_ALG_DEFAULT, &buffersize);
	cudaMalloc(&dBuffer, buffersize);
	cusparseSpMV(handle, CUSPARSE_OPERATION_NON_TRANSPOSE,
		alpha, A_desc, x_desc, beta, y_desc, A->data_type,
		CUSPARSE_MV_ALG_DEFAULT, dBuffer);
	cusparseDestroy(handle);
	Global_Free(dBuffer, DEVICE);
}

void Mv(SparseMatrixCuda<double>* A,double* x,double* y)//y=alpha*a
{Csrmv(A,x,y,&one_d,&zero_d);}
void Mv(SparseMatrixCuda<float>* A,float* x,float* y)
{Csrmv(A,x,y,&one_f,&zero_f);}

////CUDA11 TOFIX [JY]: csrgemm
////example: https://github.com/NVIDIA/CUDALibrarySamples/blob/master/cuSPARSE/spgemm/spgemm_example.c
template<class T>//T={double, float}
void SpGEMM(SparseMatrixCuda<T>* A, SparseMatrixCuda<T>* B, SparseMatrixCuda<T>* C) {//C=A*B
	assert(C != nullptr);
	//A,B must on device
	//This function will Re-allocate C on device
	//and internally allocate and erase temporary space dBuffer1, dBuffer2
	static cusparseOperation_t no_op = CUSPARSE_OPERATION_NON_TRANSPOSE;
	static T alpha = 1, beta = 0;
	cusparseHandle_t handle = nullptr; cusparseCreate(&handle);
	cusparseSpGEMMDescr_t spgemmDesc; cusparseSpGEMM_createDescr(&spgemmDesc);
	cusparseSpMatDescr_t A_desc = A->Get_SpMatCescr_t(), B_desc = B->Get_SpMatCescr_t(), C_desc = nullptr;
	cusparseCreateCsr(&C_desc, A->rows(), B->cols(), 0, nullptr, nullptr, nullptr,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I, CUSPARSE_INDEX_BASE_ZERO, A->data_type);
	
	size_t bufferSize1 = 0, bufferSize2 = 0;
	void* dBuffer1 = nullptr, * dBuffer2 = nullptr;
	//calculate size of buffer 1
	cusparseSpGEMM_workEstimation(handle, no_op, no_op, &alpha, A_desc, B_desc, &beta, C_desc,
		A->data_type, CUSPARSE_SPGEMM_DEFAULT, spgemmDesc, &bufferSize1, nullptr);
	cudaMalloc((void**)&dBuffer1, bufferSize1);
	//estimate with buffer 1
	cusparseSpGEMM_workEstimation(handle, no_op, no_op, &alpha, A_desc, B_desc, &beta, C_desc,
		A->data_type, CUSPARSE_SPGEMM_DEFAULT, spgemmDesc, &bufferSize1, dBuffer1);
	//calculate size of buffer 2
	cusparseSpGEMM_compute(handle, no_op, no_op, &alpha, A_desc, B_desc, &beta, C_desc,
		A->data_type, CUSPARSE_SPGEMM_DEFAULT, spgemmDesc, &bufferSize2, nullptr);
	cudaMalloc((void**)&dBuffer2, bufferSize2);
	//compute with buffer 2
	cusparseSpGEMM_compute(handle, no_op, no_op, &alpha, A_desc, B_desc, &beta, C_desc,
		A->data_type, CUSPARSE_SPGEMM_DEFAULT, spgemmDesc, &bufferSize2, dBuffer2);
	int64_t C_num_rows, C_num_cols, C_nnz;
	cusparseSpMatGetSize(C_desc, &C_num_rows, &C_num_cols, &C_nnz);
	C->resize(C_num_rows, C_num_cols, DEVICE);
	C->resizeNonZeros(C_nnz);
	cusparseCsrSetPointers(C_desc, C->ptr, C->col, C->val);
	cusparseSpGEMM_copy(handle, no_op, no_op, &alpha, A_desc, B_desc, &beta, C_desc,
		A->data_type, CUSPARSE_SPGEMM_DEFAULT, spgemmDesc);

	cusparseSpGEMM_destroyDescr(spgemmDesc);
	Global_Free(dBuffer1, DEVICE);
	Global_Free(dBuffer2, DEVICE);
}

void Mm(SparseMatrixCuda<double>* A_dev, SparseMatrixCuda<double>* B_dev, SparseMatrixCuda<double>* C_dev)
{
	SpGEMM(A_dev, B_dev, C_dev);
}
void Mm(SparseMatrixCuda<float>* A_dev, SparseMatrixCuda<float>* B_dev, SparseMatrixCuda<float>* C_dev)
{
	SpGEMM(A_dev, B_dev, C_dev);
}

void Axpy(double* alpha,double* x,double* y,int n){cublasDaxpy(Cublas_Handle(),n,alpha,x,1,y,1);}
void Axpy(float* alpha,float* x,float* y,int n){cublasSaxpy(Cublas_Handle(),n,alpha,x,1,y,1);}

void Copy(double* to,double* from,int n){cublasDcopy(Cublas_Handle(),n,from,1,to,1);}
void Copy(float* to,float* from,int n){cublasScopy(Cublas_Handle(),n,from,1,to,1);}

double Dot(double* x,double* y,int n){double x_dot_y=0.;cublas_status=cublasDdot(Cublas_Handle(),n,x,1,y,1,&x_dot_y);return x_dot_y;}
float Dot(float* x,float* y,int n){float x_dot_y=0.;cublas_status=cublasSdot(Cublas_Handle(),n,x,1,y,1,&x_dot_y);return x_dot_y;}

void Scale(double* alpha,double* x,int n){cublas_status=cublasDscal(Cublas_Handle(),n,alpha,x,1);}
void Scale(float* alpha,float* x,int n){cublas_status=cublasSscal(Cublas_Handle(),n,alpha,x,1);}
};
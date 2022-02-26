//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA
#ifndef __AuxFuncCuda_h__
#define __AuxFuncCuda_h__
#include "Common.h"
#include "SparseCuda.h"

namespace AuxFuncCuda
{
	//////////////////////////////////////////////////////////////////////////
	////Basic APIs
	////Data holder calculation
	DataHolder Copy_Data_Holder(const DataHolder& to, const DataHolder& from);

	////CUDA type selector
	template<class T> cudaDataType_t Cuda_Real_Type(void);

	////Memory allocation, communication between CPU and GPU
	//"Global" comes from __global__ in CUDA, for it can handle both device and host situations
	template<class T> T* Global_Malloc(int n, const DataHolder &side);
	template<class T> T* Global_Free(T*& ptr, const DataHolder &side);//Returns nullptr. When receiving nullptr, it does nothing
	template<class T> void Global_Copy_Array(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side);
	template<class T> T* Global_Malloc_And_Copy_Array(T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side);
	template<class T> void Global_Realloc_And_Copy_Array(T*& ptr_to, T*& ptr_from, int n, const DataHolder& to_side, const DataHolder& from_side);
	template<class T> void Global_Realloc_And_Copy(T*& ptr_to, T*& ptr_from, const DataHolder& to_side, const DataHolder& from_side);//only 1 element
	template<class T> void Global_Memset(T* ptr, int v, int n, const DataHolder& side);
	//////////////////////////////////////////////////////////////////////////
	////Linear algebra operations on device
	////sparse matrix-vector multiplication
	template<class T> void Csrmv(SparseMatrixCuda<T>* A_dev, T* x_dev, T* y_dev, T* alpha, T* beta);//T={double, float}, y=alpha*A*x+beta*y
	////matrix-vector multiplication
	void Mv(SparseMatrixCuda<double>* A_dev,double* x_dev,double* y_dev);
	void Mv(SparseMatrixCuda<float>* A_dev,float* x_dev,float* y_dev);
	////matrix-matrix multiplication
	template<class T> void SpGEMM(SparseMatrixCuda<T>* A_dev, SparseMatrixCuda<T>* B_dev, SparseMatrixCuda<T>* C_dev);//T={double, float}, C=A*B
	void Mm(SparseMatrixCuda<double>* A_dev, SparseMatrixCuda<double>* B_dev, SparseMatrixCuda<double>* C_dev);
	void Mm(SparseMatrixCuda<float>* A_dev, SparseMatrixCuda<float>* B_dev, SparseMatrixCuda<float>* C_dev);
	////Ax+y
	void Axpy(double* alpha,double* x,double* y,int n);
	void Axpy(float* alpha,float* x,float* y,int n);
	////copy
	void Copy(double* to,double* from,int n);
	void Copy(float* to,float* from,int n);
	////dot product
	double Dot(double* x,double* y,int n);
	float Dot(float* x,float* y,int n);
	////scale
	void Scale(double* alpha,double* x,int n);
	void Scale(float* alpha,float* x,int n);

	template<class T> cudaDataType_t Cuda_Real_Type(void) 
	{
		int siz = sizeof(T);
		if (siz == 4) { return CUDA_R_32F; }
		else if (siz == 8) { return CUDA_R_64F; }
		else { std::cerr << "[Error] AuxFuncCuda::Cuda_Type: Unknown data type\n"; return cudaDataType_t(); }
	}
};
#endif
#endif
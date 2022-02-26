//////////////////////////////////////////////////////////////////////////
// Sparse CUDA
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA
#ifndef __SparseCuda_h__
#define __SparseCuda_h__
#include <iostream>
#include "SparseFunc.h"
#include "cusparse.h"
//Note: this header should not include AuxFuncCuda.h

////The sparse matrix is stored in CRS format
template<class T> class SparseMatrixCuda
{public:
	enum DataHolder data_side;	//data on device/host
	bool realloc_on_shrink = true;
	////crs matrix
	int m=0;					////rows
	int n=0;					////cols
	int nnz=0;					////nonzeros
	int* ptr=nullptr;			////dev/host
	int* col=nullptr;			////dev/host
	T* val=nullptr;				////dev/host
	cudaDataType_t data_type;   ////host
	////color
	int color_n=0;				////number of colors
	int block_size=1;			////d
	int* color_ptr=nullptr;		////always on host
	int* color=nullptr;			////dev/host
	
	SparseMatrixCuda(const enum DataHolder &_side=DataHolder::HOST);
	SparseMatrixCuda(const SparseMatrixCuda<T>& copy,const enum DataHolder& _side=DataHolder::UNKNOWN);	//deep copy
	SparseMatrixCuda(const SparseMatrix<T>& copy,const enum DataHolder& _side=DataHolder::UNKNOWN);		//deep copy
	~SparseMatrixCuda(){Clear();}
	SparseMatrixCuda<T>& operator=(const SparseMatrixCuda<T>& copy);	////deep copy
	SparseMatrixCuda<T>& operator=(const SparseMatrix<T>& copy);		////deep copy
	operator SparseMatrix<T>&();
	void deepcopy(const SparseMatrixCuda<T>& copy);
	void deepcopy(const SparseMatrix<T>& copy);
	void Initialize_Color(int _color_n,int _block_size,const std::vector<int>& _color_ptr,const std::vector<int>& _color);	////copy color from the input arrays
	void Clear();
	size_t Memory_Size();
	//void Realloc_On_Device(int m, int n, int nnz);

	////Eigen SparseMatrix interfaces
	int rows() const {return m;}
	int cols() const {return n;}
	int outerSize() const {return m;}
	int nonZeros() const {return nnz;}
	T* valuePtr(){return val;}
	const T* valuePtr() const {return val;}
	int* outerIndexPtr(){return ptr;}
	const int* outIndexPtr() const {return ptr;}
	int* innerIndexPtr(){return col;}
	const int* innerIndexPtr() const {return col;}
	//similar to SparseMatrix::resize(). May realloc ptr, but not others. With default option, do not change data side.
	void resize(int _m, int _n, enum DataHolder new_side=DataHolder::UNKNOWN);
	void resizeNonZeros(int _nnz);//resize() must be called before, so we assume data side already set here

	////Cuda SparseMatrix interfaces
	cusparseSpMatDescr_t Get_SpMatCescr_t(void);
};


//////////////////////////////////////////////////////////////////////////
////Aux functions

////Interfaces for old code. Try not to use them in new projects.
/*template<class T> void Malloc_And_Copy_Host_To_Device(SparseMatrixCuda<T>*& ptr_dev,SparseMatrixCuda<T>*& ptr_host);
template<class T> void Malloc_And_Copy_Host_To_Device(SparseMatrixCuda<T>*& ptr_dev,SparseMatrix<T>*& ptr_host);
template<class T> void Malloc_And_Copy_Device_To_Host(SparseMatrixCuda<T>*& ptr_host, SparseMatrixCuda<T>*& ptr_dev);
template<class T> void Malloc_And_Copy_Host_To_Host(SparseMatrixCuda<T>*& ptr_host_to,SparseMatrix<T>*& ptr_host_from);
template<class T> void Malloc_And_Copy_Host_To_Host(SparseMatrixCuda<T>*& ptr_host_to,SparseMatrixCuda<T>*& ptr_host_from);
template<class T> void Malloc_And_Copy_Host_To_Host(SparseMatrix<T>*& ptr_host_to, SparseMatrixCuda<T>*& ptr_host_from);*/
////r=b-Ax
template<class T> void Residual(SparseMatrixCuda<T>* A, T* x, T* b, T* r);

#endif
#endif
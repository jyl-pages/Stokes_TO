//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "ContextCuda.h"
#include "AuxFuncCuda.h"
#include "SparseCuda.h"
#include <iostream>
using namespace AuxFuncCuda;

template<class T>SparseMatrixCuda<T>::SparseMatrixCuda(const enum DataHolder& _side) {
	data_side = _side;
	data_type = Cuda_Real_Type<T>();
}

template<class T> SparseMatrixCuda<T>::SparseMatrixCuda(const SparseMatrixCuda<T>& copy, const DataHolder& _side)
{
	data_side = _side;
	data_type = Cuda_Real_Type<T>();
	deepcopy(copy);
}
template<class T> SparseMatrixCuda<T>::SparseMatrixCuda(const SparseMatrix<T>& copy, const DataHolder& _side)
{
	data_side = _side;
	data_type = Cuda_Real_Type<T>();
	deepcopy(copy);
}

template<class T> SparseMatrixCuda<T>& SparseMatrixCuda<T>::operator=(const SparseMatrixCuda<T>& copy)
{deepcopy(copy); return *this;}

template<class T> SparseMatrixCuda<T>& SparseMatrixCuda<T>::operator=(const SparseMatrix<T>& copy)
{deepcopy(copy); return *this;}

template<class T> SparseMatrixCuda<T>::operator SparseMatrix<T>&()
{
	assert(data_side != UNKNOWN);
	int* ptr_to = Global_Malloc_And_Copy_Array<int>(ptr, m + 1, HOST, data_side);
	int* col_to = Global_Malloc_And_Copy_Array<int>(col, nnz, HOST, data_side);
	T* val_to = Global_Malloc_And_Copy_Array<T>(val, nnz, HOST, data_side);
	Eigen::Map<SparseMatrix<T> > sp_map(m, n, nnz, ptr_to, col_to, val_to);
	SparseMatrix<T>* ptr_host_to = new SparseMatrix<T>(sp_map.eval());
	return *ptr_host_to;
}

template<class T> void SparseMatrixCuda<T>::deepcopy(const SparseMatrixCuda<T>& copy)
{
	data_side = Copy_Data_Holder(data_side, copy.data_side);//If data_side is UNKNOWN, then assign it with copy.data_side
	m = copy.m; n = copy.n; nnz = copy.nnz;
	color_n = copy.color_n;
	block_size = copy.block_size;
	int block_n = m / block_size;

	int* copy_ptr = const_cast<int*>(copy.ptr);
	int* copy_col = const_cast<int*>(copy.col);
	T* copy_val = const_cast<T*>(copy.val);
	int* copy_color_ptr = const_cast<int*>(copy.color_ptr);
	int* copy_color = const_cast<int*>(copy.color);

	ptr = Global_Malloc_And_Copy_Array<int>(copy_ptr, m + 1, data_side, copy.data_side);
	col = Global_Malloc_And_Copy_Array<int>(copy_col, nnz, data_side, copy.data_side);
	val = Global_Malloc_And_Copy_Array<T>(copy_val, nnz, data_side, copy.data_side);
	color_ptr = Global_Malloc_And_Copy_Array<int>(copy_color_ptr, color_n + 1, HOST, copy.data_side);//color_ptr always on host
	color = Global_Malloc_And_Copy_Array<int>(copy_color, block_n, data_side, copy.data_side);
}

template<class T> void SparseMatrixCuda<T>::deepcopy(const SparseMatrix<T>& copy)
{
	data_side = Copy_Data_Holder(data_side, HOST);//SparseMatrixCuda is always on host
	m = copy.rows(); n = copy.cols(); nnz = copy.nonZeros();
	int* copy_ptr = const_cast<int*>(copy.outerIndexPtr());
	int* copy_col = const_cast<int*>(copy.innerIndexPtr());
	T* copy_val = const_cast<T*>(copy.valuePtr());

	ptr = Global_Malloc_And_Copy_Array<int>(copy_ptr, m + 1, data_side, HOST);
	col = Global_Malloc_And_Copy_Array<int>(copy_col, nnz, data_side, HOST);
	val = Global_Malloc_And_Copy_Array<T>(copy_val, nnz, data_side, HOST);

	////no color copy
}

template<class T> void SparseMatrixCuda<T>::Initialize_Color(int _color_n,int _block_size,const std::vector<int>& _color_ptr,const std::vector<int>& _color)
{
	//_color_ptr and _color are on host
	color_n=_color_n;block_size=_block_size;int* color_ptr_cast=const_cast<int*>(&_color_ptr[0]);
	color_ptr=Global_Malloc_And_Copy_Array<int>(color_ptr_cast,color_n+1,HOST,HOST);		//always on host

	int block_n=m/block_size;int* color_cast=const_cast<int*>(&_color[0]);
	data_side=Copy_Data_Holder(data_side,HOST);
	color=Global_Malloc_And_Copy_Array<int>(color_cast,block_n,data_side,HOST);
}

template<class T> void SparseMatrixCuda<T>::Clear()
{
	ptr = Global_Free<int>(ptr, data_side);//set as nullptr after freeing
	col = Global_Free<int>(col, data_side);
	val = Global_Free<T>(val, data_side);
	color_ptr = Global_Free<int>(color_ptr, HOST);//always on host
	color = Global_Free<int>(color, data_side);
}

template<class T> size_t SparseMatrixCuda<T>::Memory_Size()
{
	size_t size=0;
	//size+=4*sizeof(int);size+=1*sizeof(bool);	////m,n,nnz,color_n,on_device
	if(ptr)size+=(m+1)*sizeof(int);				////ptr
	if(col)size+=nnz*sizeof(int);				////col
	if(val)size+=nnz*sizeof(T);					////val
	if(color_ptr)size+=(color_n+1)*sizeof(int);	////color_ptr
	if(color)size+=m/block_size*sizeof(int);	////color
	return size;
}

template<class T> void SparseMatrixCuda<T>::resize(int _m, int _n, enum DataHolder new_side)
{
	assert(!(data_side == UNKNOWN && new_side == UNKNOWN));
	bool side_change = false;
	bool realloc = false;
	if (new_side == UNKNOWN) new_side = data_side;//do not suggest a new side, then data remain where it is
	else {
		//if side change happens, then realloc
		if (data_side == UNKNOWN) data_side = new_side;
		if (data_side != new_side) {
			realloc = true;
			side_change = true;
		}
	}
	//if size change happens, then realloc
	if (m != _m && (realloc_on_shrink || _m > m)) realloc = true;
	if (side_change) {
		Clear();
		nnz = color_n = 0;
		block_size = 1;
	}
	if (realloc) {
		Global_Free<int>(ptr, data_side);
		ptr = Global_Malloc<int>(_m + 1, new_side);
	}
	data_side = new_side;
	m = _m;
	n = _n;
}

template<class T> void SparseMatrixCuda<T>::resizeNonZeros(int _nnz)
{
	assert(data_side != UNKNOWN);
	if (nnz != _nnz && (realloc_on_shrink || _nnz > nnz)) {//realloc
		Global_Free<int>(col, data_side);
		col = Global_Malloc<int>(_nnz, data_side);
		Global_Free<T>(val, data_side);
		val = Global_Malloc<T>(_nnz, data_side);
	}
	nnz = _nnz;
}

template<class T> cusparseSpMatDescr_t SparseMatrixCuda<T>::Get_SpMatCescr_t()
{
	if (data_side != DEVICE) {
		std::cerr << "[Error] SparseMatrixCuda::valueDescr_t(): data not on device\n";
		return nullptr;
	}
	cusparseSpMatDescr_t matA = nullptr;
	cusparseStatus_t status = cusparseCreateCsr(&matA, m, n, nnz,
		ptr, col, val,
		CUSPARSE_INDEX_32I, CUSPARSE_INDEX_32I,
		CUSPARSE_INDEX_BASE_ZERO, data_type);
	return matA;
}

template class SparseMatrixCuda<float>;
template class SparseMatrixCuda<double>;

template<class T> void Residual(SparseMatrixCuda<T>* A, T* x, T* b, T* r) 
{
	assert(A->data_side == DEVICE);
	T neg_one = (T)-1; T one = (T)1;
	Copy(r, b, A->m);
	Csrmv(A, x, r, &neg_one, &one);
}

//New interfaces
template void Residual<double>(SparseMatrixCuda<double>*, double*, double*, double*);
template void Residual<float>(SparseMatrixCuda<float>*, float*, float*, float*);

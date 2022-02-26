//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include <cuda_runtime.h>
#include <cusparse.h>
#include <cublas_v2.h>
#include <cusolverSp.h>
#include <helper_functions.h>
#include <helper_cuda.h>
#include "ContextCuda.h"

////Auxiliary variables
cublasHandle_t cublas_handle=0;
cusparseHandle_t cusparse_handle=0;
cusparseMatDescr_t cusparse_mat_descr=0;
cusolverSpHandle_t cusolver_handle=0;
cublasStatus_t cublas_status;
bool cu_context_initialized=false;

bool Is_Cuda_Context_Initialized()
{return cu_context_initialized;}

cublasHandle_t Cublas_Handle()
{return cublas_handle;}

cusparseHandle_t Cusparse_Handle()
{return cusparse_handle;}

cusolverSpHandle_t Cusolver_Handle()
{return cusolver_handle;}

cusparseMatDescr_t Cusparse_Mat_Descr()
{return cusparse_mat_descr;}

void Initialize_Cuda_Context()
{
    cublasStatus_t cublasStatus;
    cublasStatus=cublasCreate(&cublas_handle);
    checkCudaErrors(cublasStatus);
    cusparseStatus_t cusparseStatus;
    cusparseStatus=cusparseCreate(&cusparse_handle);
    checkCudaErrors(cusparseStatus);
    cusparseStatus=cusparseCreateMatDescr(&cusparse_mat_descr);
    checkCudaErrors(cusparseStatus);
    cusparseSetMatType(cusparse_mat_descr,CUSPARSE_MATRIX_TYPE_GENERAL);
    cusparseSetMatIndexBase(cusparse_mat_descr,CUSPARSE_INDEX_BASE_ZERO);

	cusolverStatus_t cusolverStatus;
	cusolverStatus=cusolverSpCreate(&cusolver_handle);
	checkCudaErrors(cusolverStatus);

	cu_context_initialized=true;
}
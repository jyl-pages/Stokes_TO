//////////////////////////////////////////////////////////////////////////
// Auxiliary Function CUDA
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifdef USE_CUDA
#ifndef __ContextCuda_h__
#define __ContextCuda_h__
#include <cusparse.h>
#include <cusolverSp.h>
#include <cublas_api.h>

////This file can be included in *.cu only

bool Is_Cuda_Context_Initialized();
void Initialize_Cuda_Context();
cusparseHandle_t Cusparse_Handle();
cusolverSpHandle_t Cusolver_Handle();
cusparseMatDescr_t Cusparse_Mat_Descr();
cublasHandle_t Cublas_Handle();
#endif
#endif
#pragma once

#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "para.h"

#ifdef USE_DOUBLE
#define Dot cublasDdot
#define Axpy cublasDaxpy
#define Scal cublasDscal
#define DnSolve_buffer cusolverDnDpotrf_bufferSize
#define DnSolve_factorize cusolverDnDpotrf
#define DnSolve_solve cusolverDnDpotrs
#define DnSolve_buffer_lu cusolverDnDgetrf_bufferSize
#define DnSolve_factorize_lu cusolverDnDgetrf
#define DnSolve_solve_lu cusolverDnDgetrs
#define Batch_factorize_lu cublasDgetrfBatched
#define Batch_solve_lu cublasDgetrsBatched
#else
#define Dot cublasSdot
#define Axpy cublasSaxpy
#define Scal cublasSscal
#define DnSolve_buffer_chol cusolverDnSpotrf_bufferSize
#define DnSolve_factorize_chol cusolverDnSpotrf
#define DnSolve_solve_chol cusolverDnSpotrs
#define DnSolve_buffer_lu cusolverDnSgetrf_bufferSize
#define DnSolve_factorize_lu cusolverDnSgetrf
#define DnSolve_solve_lu cusolverDnSgetrs
#define Batch_factorize_lu cublasSgetrfBatched
#define Batch_solve_lu cublasSgetrsBatched

#endif


template<typename A, typename F>
__global__ void cwise_mapping(A v1, F f, int N)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= N) return;
	f(v1[i]);
}

template<typename A, typename F>
void cwise_mapping_wrapper(A v1, F f, int N)
{
	cwise_mapping << <((N + 63) >> 6), 64 >> > (v1, f, N);
}

template<typename A, typename B, typename F>
__global__ void cwise_mapping(A v1, B v2, F f, int N)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= N) return;
	f(v1[i], v2[i]);
}

template<typename A, typename B, typename F>
void cwise_mapping_wrapper(A v1, B v2, F f, int N)
{
	cwise_mapping << <((N + 63) >> 6), 64 >> > (v1, v2, f, N);
}

template<typename A, typename B, typename C, typename F>
__global__ void cwise_mapping(A v1, B v2, C v3, F f, int N)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= N) return;
	f(v1[i], v2[i], v3[i]);
}

template<typename A, typename B, typename C, typename F>
void cwise_mapping_wrapper(A v1, B v2, C v3, F f, int N)
{
	cwise_mapping << <((N + 63) >> 6), 64 >> > (v1, v2, v3, f, N);
}

template<typename A, typename F>
__global__ void cwise_mapping_with_idx(A v1, F f, int N)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= N) return;
	f(v1[i], i);
}

template<typename A, typename F>
void cwise_mapping_with_idx_wrapper(A v1, F f, int N)
{
	cwise_mapping_with_idx << <((N + 63) >> 6), 64 >> > (v1, f, N);
}

template<typename A, typename B, typename F>
__global__ void cwise_mapping_with_idx(A v1, B v2, F f, int N)
{
	int i = blockIdx.x*blockDim.x + threadIdx.x;
	if (i >= N) return;
	f(v1[i], v2[i], i);
}

template<typename A, typename B, typename F>
void cwise_mapping_with_idx_wrapper(A v1, B v2, F f, int N)
{
	cwise_mapping_with_idx << <((N + 63) >> 6), 64 >> > (v1, v2, f, N);
}

#pragma once
#include "BatchInverse.cuh"
#include "gpuUtils.h"
#include "cublas_v2.h"
#include "cuda_runtime_api.h"
#include <assert.h>

//#define DIRECT_SOLVE_BLOCK_INVERSE_DEBUG

template<int k>
class DirectSolveBlockInverse :public LinearMapping
{
public:
	LinearMapping *mapping;
	int dof, N;
	Scalar *d_A_val;
	Scalar **d_A_ptr;
	int *d_Ipiv;
	int *h_info, *d_info;
	Scalar **d_B_ptr;
	Scalar *d_buff;
	cublasHandle_t cublasHandle;

public:
	DirectSolveBlockInverse(LinearMapping *_mapping, Scalar *_buffer=nullptr)
	{
		//assert(_mapping->xDoF() == _mapping->yDoF());
		//dof = _mapping->xDoF();
		//N = dof / k;
		//cudaMalloc(&d_A_val, sizeof(Scalar)*dof*k);
		//cudaMalloc(&d_A_ptr, sizeof(Scalar*)*N);
		//cudaMalloc(&d_Ipiv, sizeof(int)*N*k);
		//h_info = (int*)malloc(sizeof(int)*N);
		//cudaMalloc(&d_info, sizeof(int)*N);
		//cudaMalloc(&d_B_ptr, sizeof(Scalar*)*N);

		//cudaMalloc(&d_buff, sizeof(Scalar)*LinearMappingToBlockMatrix_BufferSize<k, k>(_mapping));
		//LinearMappingToBlockMatrix<k, k>(d_A_val, _mapping, d_buff);

		//cublasCreate(&cublasHandle);


		//cudaDeviceSynchronize();
		//checkCudaErrors(cudaGetLastError());

		//test();

		//cudaDeviceSynchronize();
		//checkCudaErrors(cudaGetLastError());

		//Batch_factorize_lu(cublasHandle, k, d_A_ptr, k, d_Ipiv, d_info, N);
		//cudaDeviceSynchronize();
		//checkCudaErrors(cudaGetLastError());

		//cudaMemcpy(h_info, d_info, sizeof(int)*N, cudaMemcpyDeviceToHost);
		//cudaDeviceSynchronize();
	
		//printf("batch lu factorize info:\n");
		//for (int i = 0; i < N; i++)
		//{
		//	if (h_info[i] != 0)
		//		printf("\t%d %d\n",i, h_info[i]);
		//}

		mapping = _mapping;
		allocate(_buffer);
	}

	// work around for nvcc's bug:
	// seems like device lambda functions cannot be used in constructor
	void allocate(Scalar *_buffer)
	{
		assert(mapping->xDoF() == mapping->yDoF());
		dof = mapping->xDoF();
		N = dof / k;
		cudaMalloc(&d_A_val, sizeof(Scalar)*dof*k);
		cudaMalloc(&d_A_ptr, sizeof(Scalar*)*N);
		cudaMalloc(&d_Ipiv, sizeof(int)*N*k);
		h_info = (int*)malloc(sizeof(int)*N);
		cudaMalloc(&d_info, sizeof(int)*N);
		cudaMalloc(&d_B_ptr, sizeof(Scalar*)*N);

		if (_buffer == nullptr)
			cudaMalloc(&d_buff, sizeof(Scalar)*LinearMappingToBlockMatrix_BufferSize<k, k>(mapping));
		else
			d_buff = _buffer;

		Scalar *temp = d_A_val;
		auto assign = [=] __device__(Scalar*& tv, int ind) { tv = temp + ind * k*k; };
		cwise_mapping_with_idx_wrapper(d_A_ptr, assign, N);

		cublasCreate(&cublasHandle);
	}

	int xDoF() override
	{
		return dof;
	}

	int yDoF() override
	{
		return dof;
	}

	void applyMapping(Scalar *Ap, Scalar *p) override
	{
		auto assign = [=] __device__(Scalar*& tv, int ind) { tv = Ap + ind * k; };
		cwise_mapping_with_idx_wrapper(d_B_ptr, assign, N);
		cudaMemcpy(Ap, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
		int info;
		Batch_solve_lu(cublasHandle, CUBLAS_OP_N, k, 1, d_A_ptr, k, d_Ipiv, d_B_ptr, k, &info, N);
		//cudaDeviceSynchronize();
		//printf("batch lu solve info: %d\n", info);
	}

	void update() override
	{
		LinearMappingToBlockMatrix<k, k>(d_A_val, mapping, d_buff);

		Batch_factorize_lu(cublasHandle, k, d_A_ptr, k, d_Ipiv, d_info, N);
#ifdef DIRECT_SOLVE_BLOCK_INVERSE_DEBUG
		cudaMemcpy(h_info, d_info, sizeof(int)*N, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();

		printf("batch lu factorize info:\n");
		for (int i = 0; i < N; i++)
		{
			if (h_info[i] != 0)
				printf("\t%d %d\n", i, h_info[i]);
		}
#endif
	}
};

#pragma once
#include "linearMapping.h"
#include "cuda_runtime_api.h"
#include "gpuUtils.h"

//#define BATCH_INVERSE_DEBUG

#ifdef BATCH_INVERSE_DEBUG
#include <cstdio>
#endif

// only for size<32
template<int kx, int ky>
__global__ void BlockMatrixReorderKernel(Scalar *out, Scalar *in, int N)
{
	__shared__ Scalar shared[kx * 64];
	__shared__ Scalar start, end;
	const int block = blockIdx.x;
	const int idx = threadIdx.x;
	if (idx == 0) start = 0;
	__syncthreads();
	for (int i = 0; i < ky; i++)
	{
		Scalar * const buf = shared + (i & 1)*(kx * 32);
		for (int j = 0; j < kx; j++)
			if (block * 32 * ky + i * 32 + idx < N)
				buf[32 * j + idx] = in[j*N + block * 32 * ky + i * 32 + idx];
		if (idx == 0)
		{
			end = start + kx;
			if (i == 0 && 32 % ky != 0) end--;
			if (i == ky - 1 && 32 % ky != 0) end++;
		}
		__syncthreads();
		for (int j = start; j < end; j++)
		{
			int ind = j * 32 + idx;
			int mat_ind = ind / (ky*kx);
			int col_ind = (ind / ky) % kx;
			int row_ind = ind % ky;
			int buf_row_ind = (mat_ind*ky + row_ind) % 64;
			if (block * 32 * kx*ky + ind < N*kx)
				out[block * 32 * kx*ky + ind] = shared[buf_row_ind / 32 * (32 * kx) + buf_row_ind % 32 + 32 * col_ind];
		}
		if (idx == 0) start = end;
	}
}

template<int kx, int ky>
int LinearMappingToBlockMatrix_BufferSize(LinearMapping *mapping)
{
	return mapping->xDoF() + kx * mapping->yDoF();
}

template<int kx, int ky>
int LinearMappingToBlockMatrix_ResultSize(LinearMapping *mapping)
{
	return kx * mapping->yDoF();
}

template<int kx, int ky>
void LinearMappingToBlockMatrix(Scalar *result, LinearMapping *mapping, Scalar *buffer)
{
	int xdof = mapping->xDoF(), ydof = mapping->yDoF();

	for (int i = 0; i < kx; i++)
	{
		auto assign = [=] __device__(Scalar& tv, int ind) { tv = (ind%kx == i) ? (Scalar)1 : (Scalar)0; };
		cwise_mapping_with_idx_wrapper(buffer, assign, xdof);
		mapping->applyMapping(buffer + xdof + i * ydof, buffer);
	}

	//LinearMappingToBlockMatrixHelper<kx - 1, kx>(result, mapping, buffer);

	BlockMatrixReorderKernel<kx, ky> << <(ydof + 31) / 32, 32 >> > (result, buffer + xdof, ydof);
#ifdef BATCH_INVERSE_DEBUG
	Scalar *watch = (Scalar*)malloc(sizeof(Scalar)*LinearMappingToBlockMatrix_ResultSize<kx,ky>(mapping));
	Scalar *mat = (Scalar*)malloc(sizeof(Scalar)*kx*ky);
	cudaMemcpy(watch, result, sizeof(Scalar)*LinearMappingToBlockMatrix_ResultSize<kx,ky>(mapping), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	FILE *f = fopen("BATCH_INVERSE_DEBUG.txt", "w");

	for (int i = 0; i < xdof / kx; i++)
	{
		for (int j = 0; j < kx; j++)
			for (int k = 0; k < ky; k++)
				mat[k*kx + j] = watch[i*kx*ky + j * ky + k];
		for (int k = 0; k < ky; k++)
		{
			for (int j = 0; j < kx; j++)
				fprintf(f, "%lf ", mat[k*kx + j]);
			fprintf(f, "\n");
		}
		fprintf(f, "\n");
	}

	fclose(f);
	free(watch);
	free(mat);
#endif
}

//template<int T, int kx>
//void LinearMappingToBlockMatrixHelper(Scalar *result, LinearMapping *mapping, Scalar *buffer)
//{
//	int xdof = mapping->xDoF(), ydof = mapping->yDoF();
//	auto assign = [=] __device__(Scalar& tv, int ind) { tv = (ind%kx == T) ? (Scalar)1 : (Scalar)0; };
//	//auto assign = [=] __device__(Scalar& tv, int ind) { tv = (Scalar)ind; };
//	cwise_mapping_with_idx_wrapper(buffer, assign, xdof);
//	mapping->applyMapping(result + T * ydof, buffer);
//	LinearMappingToBlockMatrixHelper<T - 1, kx>(result, mapping, buffer);
//}
//
//
//template<int kx>
//void LinearMappingToBlockMatrixHelper<0, kx>(Scalar *result, LinearMapping *mapping, Scalar *buffer)
//{
//	int xdof = mapping->xDoF(), ydof = mapping->yDoF();
//	auto assign = [=] __device__(Scalar& tv, int ind) { tv = (ind%kx == 0) ? (Scalar)1 : (Scalar)0; };
//	//auto assign = [=] __device__(Scalar& tv, int ind) { tv = (Scalar)ind; };
//	cwise_mapping_with_idx_wrapper(buffer, assign, xdof);
//	mapping->applyMapping(result + 0 * ydof, buffer);
//}

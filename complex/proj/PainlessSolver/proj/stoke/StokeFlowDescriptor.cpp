#include "StokeFlowDescriptor.h"
#include "cuda_runtime_api.h"

void StokeFlowDescriptor::init(int _Nx, int _Ny)
{
	Nx = _Nx;
	Ny = _Ny;
	grid.init(_Nx, _Ny);
	size = grid.cell_size() + grid.face_size(0) + grid.face_size(1);
	h_vol = (Scalar*)malloc(sizeof(Scalar)*size);
	h_fixed = (bool*)malloc(sizeof(bool)*size);
	h_penalty = (Scalar*)malloc(sizeof(Scalar)*size);
	cudaMalloc(&d_vol, sizeof(Scalar)*size);
	cudaMalloc(&d_fixed, sizeof(bool)*size);
	cudaMalloc(&d_penalty, sizeof(Scalar)*size);
}

void StokeFlowDescriptor::setFixed(bool *fixed)
{
	memcpy(h_fixed, fixed, sizeof(bool)*size);
}

void StokeFlowDescriptor::setVol(Scalar *vol)
{
	memcpy(h_vol, vol, sizeof(Scalar)*size);
}

void StokeFlowDescriptor::setPenalty(Scalar *penalty)
{
	memcpy(h_penalty, penalty, sizeof(Scalar)*size);
}

void StokeFlowDescriptor::finish()
{
	grid2DOperator::ApplyFaceRedudantFixed(grid, d_fixed + grid.cell_size(), 0);
	grid2DOperator::ApplyFaceRedudantFixed(grid, d_fixed + grid.cell_size() + grid.face_size(0), 1);
	grid2DOperator::ApplyMask(grid, d_vol, d_fixed, size);
	grid2DOperator::ApplyMask(grid, d_penalty , d_fixed, size);
}

void StokeFlowDescriptor::toHost()
{
	cudaMemcpy(h_vol, d_vol, sizeof(Scalar)*size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fixed, d_fixed, sizeof(bool)*size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_penalty, d_penalty, sizeof(Scalar)*size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void StokeFlowDescriptor::toDevice()
{
	cudaMemcpy(d_vol, h_vol, sizeof(Scalar)*size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_fixed, h_fixed, sizeof(bool)*size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_penalty, h_penalty, sizeof(Scalar)*size, cudaMemcpyHostToDevice);
}

#include "PoissonDescriptor.h"
#include "cuda_runtime_api.h"

void PoissonDescriptor::init(int _Nx, int _Ny)
{
	Nx = _Nx;
	Ny = _Ny;
	grid.init(_Nx, _Ny);
	size = grid.cell_size();
	h_vol = (Scalar*)malloc(sizeof(Scalar)*(grid.face_size(0)+grid.face_size(1)));
	h_fixed = (bool*)malloc(sizeof(bool)*size);
	cudaMalloc(&d_vol, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
	cudaMalloc(&d_fixed, sizeof(bool)*size);
}

void PoissonDescriptor::setFixed(bool * fixed)
{
	memcpy(h_fixed, fixed, sizeof(bool)*size);
}

void PoissonDescriptor::setVol(Scalar * vol)
{
	memcpy(h_vol, vol, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
}

void PoissonDescriptor::finish()
{
}

void PoissonDescriptor::toHost()
{
	cudaMemcpy(h_vol, d_vol, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)), cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fixed, d_fixed, sizeof(bool)*size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void PoissonDescriptor::toDevice()
{
	cudaMemcpy(d_vol, h_vol, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)), cudaMemcpyHostToDevice);
	cudaMemcpy(d_fixed, h_fixed, sizeof(bool)*size, cudaMemcpyHostToDevice);
}

#include "PoissonDescriptor3D.h"
#include "cuda_runtime_api.h"

void PoissonDescriptor3D::init(int _Nx, int _Ny, int _Nz)
{
	Nx = _Nx;
	Ny = _Ny;
	Nz = _Nz;
	grid.init(_Nx, _Ny, _Nz);
	size = grid.cell_size();
	fsize = grid.face_size(0) + grid.face_size(1) + grid.face_size(2);
	h_vol = (Scalar*)malloc(sizeof(Scalar)*fsize);
	h_fixed = (bool*)malloc(sizeof(bool)*size);
	cudaMalloc(&d_vol, sizeof(Scalar)*fsize);
	cudaMalloc(&d_fixed, sizeof(bool)*size);
}

void PoissonDescriptor3D::setFixed(bool * fixed)
{
	memcpy(h_fixed, fixed, sizeof(bool)*size);
}

void PoissonDescriptor3D::setVol(Scalar * vol)
{
	memcpy(h_vol, vol, sizeof(Scalar)*fsize);
}

void PoissonDescriptor3D::finish()
{
}

void PoissonDescriptor3D::toHost()
{
	cudaMemcpy(h_vol, d_vol, sizeof(Scalar)*fsize, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fixed, d_fixed, sizeof(bool)*size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void PoissonDescriptor3D::toDevice()
{
	cudaMemcpy(d_vol, h_vol, sizeof(Scalar)*fsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_fixed, h_fixed, sizeof(bool)*size, cudaMemcpyHostToDevice);
}

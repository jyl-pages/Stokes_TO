#include "StokeFlowDescriptor3D.h"
#include "cuda_runtime_api.h"
#include "gpuUtils.h"
#include <memory>

void StokeFlowDescriptor3D::init(int _Nx, int _Ny, int _Nz)
{
	Nx = _Nx;
	Ny = _Ny;
	Nz = _Nz;
	grid.init(_Nx, _Ny, _Nz);
	size = grid.cell_size() + grid.face_size(0) + grid.face_size(1) + grid.face_size(2);
	tsize=size+ grid.edge_size(0) + grid.edge_size(1) + grid.edge_size(2);
	h_vol = (Scalar*)malloc(sizeof(Scalar)*size);
	h_fixed = (bool*)malloc(sizeof(bool)*size);
	h_penalty = (Scalar*)malloc(sizeof(Scalar)*size);
	cudaMalloc(&d_vol, sizeof(Scalar)*size);
	cudaMalloc(&d_fixed, sizeof(bool)*size);
	cudaMalloc(&d_penalty, sizeof(Scalar)*size);
}

void StokeFlowDescriptor3D::setFixed(bool *fixed)
{
	memcpy(h_fixed, fixed, sizeof(bool)*size);
}

void StokeFlowDescriptor3D::setVol(Scalar *vol)
{
	memcpy(h_vol, vol, sizeof(Scalar)*size);
}

void StokeFlowDescriptor3D::setPenalty(Scalar *penalty)
{
	memcpy(h_penalty, penalty, sizeof(Scalar)*size);
}

void StokeFlowDescriptor3D::finish()
{
	int offset = grid.cell_size();
	for (int d = 0; d < 3; d++)
	{
		grid3D grid_temp = grid;
		auto F = [=]__device__(bool& tfixed, int i) { int x, y, z; grid_temp.ind_face(i, d, x, y, z); tfixed |= !grid_temp.face_valid(x, y, z, d); };
		cwise_mapping_with_idx_wrapper(d_fixed + offset, F, grid.face_size(d));
		offset += grid.face_size(d);
	}
	auto fixed = [] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	cwise_mapping_wrapper(d_vol, d_fixed, fixed, size);
	cwise_mapping_wrapper(d_penalty, d_fixed, fixed, size);
}

void StokeFlowDescriptor3D::toHost()
{
	cudaMemcpy(h_vol, d_vol, sizeof(Scalar)*size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fixed, d_fixed, sizeof(bool)*size, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_penalty, d_penalty, sizeof(Scalar)*size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

void StokeFlowDescriptor3D::toDevice()
{
	cudaMemcpy(d_vol, h_vol, sizeof(Scalar)*size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_fixed, h_fixed, sizeof(bool)*size, cudaMemcpyHostToDevice);
	cudaMemcpy(d_penalty, h_penalty, sizeof(Scalar)*size, cudaMemcpyHostToDevice);
}

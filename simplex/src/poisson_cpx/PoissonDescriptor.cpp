#include "PoissonDescriptor.h"
#include "cuda_runtime_api.h"

template<int d>
void PoissonDescriptor<d>::init(const VectorDi& _N)
{
	N = _N;
	if constexpr (d == 2) {
		grid.init(N[0], N[1]);
	}
	else {
		grid.init(N[0], N[1], N[2]);
	}
	size = grid.cell_size();
	fsize = 0; for (int i = 0; i < d; i++) fsize += grid.face_size(i);
	h_vol = (Scalar*)malloc(sizeof(Scalar) * fsize);
	h_fixed = (bool*)malloc(sizeof(bool) * size);
	cudaMalloc(&d_vol, sizeof(Scalar) * fsize);
	cudaMalloc(&d_fixed, sizeof(bool) * size);
}

template<int d>
void PoissonDescriptor<d>::setFixed(bool * fixed)
{
	memcpy(h_fixed, fixed, sizeof(bool)*size);
}

template<int d>
void PoissonDescriptor<d>::setVol(Scalar * vol)
{
	memcpy(h_vol, vol, sizeof(Scalar) * fsize);
}

template<int d>
void PoissonDescriptor<d>::finish()
{
}

template<int d>
void PoissonDescriptor<d>::toHost()
{
	cudaMemcpy(h_vol, d_vol, sizeof(Scalar) * fsize, cudaMemcpyDeviceToHost);
	cudaMemcpy(h_fixed, d_fixed, sizeof(bool) * size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

template<int d>
void PoissonDescriptor<d>::toDevice()
{
	cudaMemcpy(d_vol, h_vol, sizeof(Scalar) * fsize, cudaMemcpyHostToDevice);
	cudaMemcpy(d_fixed, h_fixed, sizeof(bool) * size, cudaMemcpyHostToDevice);
}

template class PoissonDescriptor<2>;
template class PoissonDescriptor<3>;
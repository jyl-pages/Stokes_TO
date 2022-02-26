#include "FluidEnergyGrad.h"
#include "gpuUtils.h"
#include <cstdlib>
#include <memory.h>

void FluidEnergyGrad::init(StokeFlowDescriptor *_descr)
{
	descr = _descr;
	Nx = descr->Nx;
	Ny = descr->Ny;

	grid.init(Nx, Ny);

	dof = grid.cell_size() + grid.face_size(0) + grid.face_size(1);
	cudaMalloc(&d_temp, sizeof(Scalar)*dof);

	cublasCreate(&cublasHandle);
}

int FluidEnergyGrad::xDoF()
{
	return dof;
}

int FluidEnergyGrad::yDoF()
{
	return dof;
}

void FluidEnergyGrad::applyMapping(Scalar *Ap, Scalar *p)
{
	int cell_off = 0;
	int face_x_off = cell_off + grid.cell_size();
	int face_y_off = face_x_off + grid.face_size(0);

	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;
	Scalar mu = descr->lap_coeff*(Scalar)1;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [=] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };

	// lap v
	grid2DOperator::faceLapMapping(grid, d_temp + face_x_off, d_temp + face_y_off, p + face_x_off, p + face_y_off);

	//mapping.ApplyMask(d_temp + face_x_off, d_fixed + face_x_off, mapping.face_size(0));
	//mapping.ApplyMask(d_temp + face_y_off, d_fixed + face_y_off, mapping.face_size(1));
	cwise_mapping_wrapper(d_temp + face_x_off, descr->d_fixed + face_x_off, fix, grid.face_size(0));
	cwise_mapping_wrapper(d_temp + face_y_off, descr->d_fixed + face_y_off, fix, grid.face_size(1));

	Axpy(cublasHandle, grid.face_size(0), &mu, d_temp + face_x_off, 1, Ap + face_x_off, 1);
	Axpy(cublasHandle, grid.face_size(1), &mu, d_temp + face_y_off, 1, Ap + face_y_off, 1);

	// penalty
	cudaMemcpy(d_temp, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

	cwise_mapping_wrapper(d_temp, descr->d_penalty, multi, dof);

	Axpy(cublasHandle, dof, &one, d_temp, 1, Ap, 1);
}

#include "PoissonMapping.h"
#include "gpuUtils.h"
#include <cstdlib>
#include <memory.h>

void PoissonMapping::init(PoissonDescriptor<2> *_descr)
{
	descr = _descr;
	Nx = descr->N[0];
	Ny = descr->N[1];

	grid.init(Nx, Ny);

	dof = grid.cell_size();

	cudaMalloc(&d_temp, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1)));
}

int PoissonMapping::xDoF()
{
	return dof;
}

int PoissonMapping::yDoF()
{
	return dof;
}

void PoissonMapping::applyMapping(Scalar * Ap, Scalar * p)
{
	int face_x_off = 0;
	int face_y_off = face_x_off + grid.face_size(0);

	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [=] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };
	auto inv = [=]__device__(Scalar& tv, Scalar tvol) { tv /= tvol; };
	auto neg = [=]__device__(Scalar& tv) { tv = -tv; };

	grid2DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, p);

	cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));
	cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
	cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

	grid2DOperator::D1Mapping(grid, Ap, d_temp + face_x_off, d_temp + face_y_off);

	cwise_mapping_wrapper(Ap, descr->d_fixed, fix, grid.cell_size());

	cwise_mapping_wrapper(Ap, neg, grid.cell_size());
}

void PoissonMappingFixed::init(PoissonDescriptor<2> *_descr)
{
	descr = _descr;
	Nx = descr->N[0];
	Ny = descr->N[1];

	grid.init(Nx, Ny);

	dof = grid.cell_size();

	cudaMalloc(&d_temp, sizeof(Scalar)*(dof + grid.face_size(0) + grid.face_size(1)));
}

int PoissonMappingFixed::xDoF()
{
	return dof;
}

int PoissonMappingFixed::yDoF()
{
	return dof;
}

void PoissonMappingFixed::applyMapping(Scalar * Ap, Scalar * p)
{
	int face_x_off = 0;
	int face_y_off = face_x_off + grid.face_size(0);
	int cell_off = face_y_off + grid.face_size(1);

	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [=] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };
	auto neg = [=]__device__(Scalar& tv) { tv = -tv; };

	cudaMemcpy(d_temp + cell_off, p, sizeof(Scalar)*grid.cell_size(), cudaMemcpyDeviceToDevice);
	cwise_mapping_wrapper(d_temp + cell_off, descr->d_fixed, fix, grid.cell_size());

	grid2DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + cell_off);

	cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));
	cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
	cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

	grid2DOperator::D1Mapping(grid, Ap, d_temp + face_x_off, d_temp + face_y_off);

	cwise_mapping_wrapper(Ap, descr->d_fixed, fix, grid.cell_size());

	cwise_mapping_wrapper(Ap, neg, grid.cell_size());

	auto cond_add = [=]__device__(Scalar &tv1, Scalar tv2, bool tfixed) { if (tfixed) tv1 += tv2; };
	cwise_mapping_wrapper(Ap, p, descr->d_fixed, cond_add, grid.cell_size());
}

#include "PoissonMapping3D.h"
#include "gpuUtils.h"
#include <cstdlib>
#include <memory.h>

void PoissonMapping3D::init(PoissonDescriptor<3> *_descr)
{
	descr = _descr;
	Nx = descr->N[0];
	Ny = descr->N[1];
	Nz = descr->N[2];

	grid.init(Nx, Ny, Nz);

	dof = grid.cell_size();

	cudaMalloc(&d_temp, sizeof(Scalar)*(grid.face_size(0) + grid.face_size(1) + grid.face_size(2)));
}

int PoissonMapping3D::xDoF()
{
	return dof;
}

int PoissonMapping3D::yDoF()
{
	return dof;
}

void PoissonMapping3D::applyMapping(Scalar * Ap, Scalar * p)
{
	int cell_off = 0;
	int face_x_off = 0;
	int face_y_off = face_x_off + grid.face_size(0);
	int face_z_off = face_y_off + grid.face_size(1);

	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [=] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };
	auto inv = [=]__device__(Scalar& tv, Scalar tvol) { tv /= tvol; };
	auto neg = [=]__device__(Scalar& tv) { tv = -tv; };

	grid3DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off, p + cell_off);

	cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
	cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));
	cwise_mapping_wrapper(d_temp + face_z_off, descr->d_vol + face_z_off, multi, grid.face_size(2));

	grid3DOperator::D2Mapping(grid, Ap + cell_off, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off);

	cwise_mapping_wrapper(Ap + cell_off, descr->d_fixed + cell_off, fix, grid.cell_size());

	cwise_mapping_wrapper(Ap + cell_off, neg, grid.cell_size());
}


void PoissonMapping3DFixed::init(PoissonDescriptor<3> *_descr)
{
	descr = _descr;
	Nx = descr->N[0];
	Ny = descr->N[1];
	Nz = descr->N[2];

	grid.init(Nx, Ny, Nz);

	dof = grid.cell_size();

	cudaMalloc(&d_temp, sizeof(Scalar)*(dof + grid.face_size(0) + grid.face_size(1) + grid.face_size(2)));
}

int PoissonMapping3DFixed::xDoF()
{
	return dof;
}

int PoissonMapping3DFixed::yDoF()
{
	return dof;
}

void PoissonMapping3DFixed::applyMapping(Scalar * Ap, Scalar * p)
{
	int face_x_off = 0;
	int face_y_off = face_x_off + grid.face_size(0);
	int face_z_off = face_y_off + grid.face_size(1);
	int cell_off = face_z_off + grid.face_size(2);


	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [=] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };
	auto inv = [=]__device__(Scalar& tv, Scalar tvol) { tv /= tvol; };
	auto neg = [=]__device__(Scalar& tv) { tv = -tv; };

	// -div grad
	{
		cudaMemcpy(d_temp + cell_off, p, sizeof(Scalar)*grid.cell_size(), cudaMemcpyDeviceToDevice);
		cwise_mapping_wrapper(d_temp + cell_off, descr->d_fixed, fix, grid.cell_size());

		grid3DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off, d_temp + cell_off);

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_z_off, descr->d_vol + face_z_off, multi, grid.face_size(2));

		grid3DOperator::D2Mapping(grid, Ap, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off);

		cwise_mapping_wrapper(Ap, descr->d_fixed, fix, grid.cell_size());

		cwise_mapping_wrapper(Ap, neg, grid.cell_size());
	}

	// definite fix
	{
		auto cond_add = [=]__device__(Scalar &tv1, Scalar tv2, bool tfixed) { if (tfixed) tv1 += tv2; };
		cwise_mapping_wrapper(Ap, p, descr->d_fixed, cond_add, dof);
	}
}

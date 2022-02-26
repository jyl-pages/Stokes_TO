#include "PoissonSubSystem.h"
#include "SubSystemUtils.cuh"
#include "cuda_runtime.h"
#include "gpuUtils.h"

PoissonDownSample::PoissonDownSample(PoissonDescriptor *_f_desr, PoissonDescriptor *_c_desr)
{
	f_descr = _f_desr;
	c_descr = _c_desr;

	f_Nx = f_descr->Nx;
	f_Ny = f_descr->Ny;
	c_Nx = c_descr->Nx;
	c_Ny = c_descr->Ny;

	xdof = f_descr->size;
	ydof = c_descr->size;
}

int PoissonDownSample::xDoF()
{
	return xdof;
}

int PoissonDownSample::yDoF()
{
	return ydof;
}

void PoissonDownSample::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid2D& c_grid = c_descr->grid;
	const grid2D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };

	CellDownSampleKernel << <dim3(c_Nx / 8, c_Ny / 8), dim3(8, 8) >> > \
		(Ap, p, c_grid, f_grid, cell_kernel);

	cwise_mapping_wrapper(Ap, c_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

PoissonUpSample::PoissonUpSample(PoissonDescriptor *_f_desr, PoissonDescriptor *_c_desr)
{
	f_descr = _f_desr;
	c_descr = _c_desr;

	f_Nx = f_descr->Nx;
	f_Ny = f_descr->Ny;
	c_Nx = c_descr->Nx;
	c_Ny = c_descr->Ny;

	xdof = c_descr->size;
	ydof = f_descr->size;

}

int PoissonUpSample::xDoF()
{
	return xdof;
}

int PoissonUpSample::yDoF()
{
	return ydof;
}

void PoissonUpSample::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid2D& c_grid = c_descr->grid;
	const grid2D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar v, int x, int y)->Scalar { return v; };

	CellUpSampleKernel << <dim3(f_Nx / 16, f_Ny / 16), dim3(16, 16) >> > \
		(Ap, p, f_grid, c_grid, cell_kernel);

	cwise_mapping_wrapper(Ap, f_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

void updateSubSystem(PoissonDescriptor& c_descr, const PoissonDescriptor& f_descr)
{
	const grid2D& c_grid = c_descr.grid;
	const grid2D& f_grid = f_descr.grid;

	//auto fixed_cell_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 && v01 && v10 && v11; };
	//auto fixed_face_x_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 && v01 && v10 && v11; };
	//auto fixed_face_y_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 && v01 && v10 && v11; };

	auto fixed_cell_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 || v01 || v10 || v11; };

	cudaMemset(c_descr.d_fixed, 0, sizeof(bool)*c_grid.cell_size());

	CellDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_fixed, f_descr.d_fixed, c_grid, f_grid, fixed_cell_kernel);

	auto vol_face_x_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v10; };
	auto vol_face_y_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01; };

	cudaMemset(c_descr.d_vol, 0, sizeof(Scalar)*(c_grid.face_size(0) + c_grid.face_size(1)));

	FaceXDownSampleKernel << <dim3(c_grid.Nx / 8 + 1, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_vol, f_descr.d_vol, c_grid, f_grid, vol_face_x_kernel, (Scalar)0);
	FaceYDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8 + 1), dim3(8, 8) >> > \
		(c_descr.d_vol + c_grid.face_size(0), f_descr.d_vol + f_grid.face_size(0), c_grid, f_grid, vol_face_y_kernel, (Scalar)0);

	//cwise_mapping_wrapper(c_descr.d_vol, c_descr.d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, c_descr.size);
}

PoissonDescriptor createSubSystem(const PoissonDescriptor& f_descr)
{
	PoissonDescriptor c_descr;
	c_descr.init(f_descr.Nx / 2, f_descr.Ny / 2);

	updateSubSystem(c_descr, f_descr);

	return c_descr;
}

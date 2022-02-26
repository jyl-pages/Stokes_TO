#include "StokeFlowSubSystem.h"
#include "SubSystemUtils.cuh"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"
#include "gpuUtils.h"

void StokeFlowDownSample::init(StokeFlowDescriptor *_f_desr, StokeFlowDescriptor *_c_desr)
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

int StokeFlowDownSample::xDoF()
{
	return xdof;
}

int StokeFlowDownSample::yDoF()
{
	return ydof;
}

void StokeFlowDownSample::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid2D& c_grid = c_descr->grid;
	const grid2D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__ (Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };
	auto face_x_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };
	auto face_y_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };

	CellDownSampleKernel << <dim3(c_Nx / 8, c_Ny / 8), dim3(8, 8) >> > \
		(Ap, p, c_grid, f_grid, cell_kernel);
	FaceXDownSampleKernel << <dim3(c_Nx / 8 + 1, c_Ny / 8), dim3(8, 8) >> > \
		(Ap + c_grid.cell_size(), p + f_grid.cell_size(), c_grid, f_grid, face_x_kernel, (Scalar)0);
	FaceYDownSampleKernel << <dim3(c_Nx / 8, c_Ny / 8 + 1), dim3(8, 8) >> > \
		(Ap + c_grid.cell_size() + c_grid.face_size(0), p + f_grid.cell_size() + f_grid.face_size(0), c_grid, f_grid, face_y_kernel, (Scalar)0);

	cwise_mapping_wrapper(Ap, c_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

void StokeFlowUpSample::init(StokeFlowDescriptor *_f_desr, StokeFlowDescriptor *_c_desr)
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

int StokeFlowUpSample::xDoF()
{
	return xdof;
}

int StokeFlowUpSample::yDoF()
{
	return ydof;
}

void StokeFlowUpSample::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid2D& c_grid = c_descr->grid;
	const grid2D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar v, int x, int y)->Scalar { return v; };
	auto face_x_kernel = [=] __device__(Scalar v, int x, int y)->Scalar { return v; };
	auto face_y_kernel = [=] __device__(Scalar v, int x, int y)->Scalar { return v; };

	CellUpSampleKernel << <dim3(f_Nx / 16, f_Ny / 16), dim3(16, 16) >> > \
		(Ap, p, f_grid, c_grid, cell_kernel);
	FaceXUpSampleKernel << <dim3(f_Nx / 16 + 1, f_Ny / 16), dim3(16, 16) >> > \
		(Ap + f_grid.cell_size(), p + c_grid.cell_size(), f_grid, c_grid, face_x_kernel);
	FaceYUpSampleKernel << <dim3(f_Nx / 16, f_Ny / 16 + 1), dim3(16, 16) >> > \
		(Ap + f_grid.cell_size() + f_grid.face_size(0), p + c_grid.cell_size() + c_grid.face_size(0), f_grid, c_grid, face_y_kernel);

	cwise_mapping_wrapper(Ap, f_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

void updateSubSystem(StokeFlowDescriptor& c_descr, const StokeFlowDescriptor& f_descr)
{
	const grid2D& c_grid = c_descr.grid;
	const grid2D& f_grid = f_descr.grid;

	//auto fixed_cell_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 && v01 && v10 && v11; };
	//auto fixed_face_x_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 && v01 && v10 && v11; };
	//auto fixed_face_y_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 && v01 && v10 && v11; };

	auto fixed_cell_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 || v01 || v10 || v11; };
	auto fixed_face_x_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 || v01 || v10 || v11; };
	auto fixed_face_y_kernel = [=] __device__(bool v00, bool v01, bool v10, bool v11)->bool { return v00 || v01 || v10 || v11; };

	cudaMemset(c_descr.d_fixed, 0, sizeof(bool)*c_descr.size);

	CellDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_fixed, f_descr.d_fixed, c_grid, f_grid, fixed_cell_kernel);
	FaceXDownSampleKernel << <dim3(c_grid.Nx / 8 + 1, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_fixed + c_grid.cell_size(), f_descr.d_fixed + f_grid.cell_size(), c_grid, f_grid, fixed_face_x_kernel, false);
	FaceYDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8 + 1), dim3(8, 8) >> > \
		(c_descr.d_fixed + c_grid.cell_size() + c_grid.face_size(0), f_descr.d_fixed + f_grid.cell_size() + f_grid.face_size(0), c_grid, f_grid, fixed_face_y_kernel, false);

	auto vol_cell_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };
	auto vol_face_x_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v10; };
	auto vol_face_y_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01; };

	cudaMemset(c_descr.d_vol, 0, sizeof(Scalar)*c_descr.size);

	CellDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_vol, f_descr.d_vol, c_grid, f_grid, vol_cell_kernel);
	FaceXDownSampleKernel << <dim3(c_grid.Nx / 8 + 1, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_vol + c_grid.cell_size(), f_descr.d_vol + f_grid.cell_size(), c_grid, f_grid, vol_face_x_kernel, (Scalar)0);
	FaceYDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8 + 1), dim3(8, 8) >> > \
		(c_descr.d_vol + c_grid.cell_size() + c_grid.face_size(0), f_descr.d_vol + f_grid.cell_size() + f_grid.face_size(0), c_grid, f_grid, vol_face_y_kernel, (Scalar)0);

	//cwise_mapping_wrapper(c_descr.d_vol, c_descr.d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, c_descr.size);

	auto penalty_cell_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };
	auto penalty_face_x_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };
	auto penalty_face_y_kernel = [=] __device__(Scalar v00, Scalar v01, Scalar v10, Scalar v11)->Scalar { return v00 + v01 + v10 + v11; };

	cudaMemset(c_descr.d_penalty, 0, sizeof(Scalar)*c_descr.size);

	CellDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_penalty, f_descr.d_penalty, c_grid, f_grid, penalty_cell_kernel);
	FaceXDownSampleKernel << <dim3(c_grid.Nx / 8 + 1, c_grid.Ny / 8), dim3(8, 8) >> > \
		(c_descr.d_penalty + c_grid.cell_size(), f_descr.d_penalty + f_grid.cell_size(), c_grid, f_grid, penalty_face_x_kernel, (Scalar)0);
	FaceYDownSampleKernel << <dim3(c_grid.Nx / 8, c_grid.Ny / 8 + 1), dim3(8, 8) >> > \
		(c_descr.d_penalty + c_grid.cell_size() + c_grid.face_size(0), f_descr.d_penalty + f_grid.cell_size() + f_grid.face_size(0), c_grid, f_grid, penalty_face_y_kernel, (Scalar)0);

	//cwise_mapping_wrapper(c_descr.d_penalty, c_descr.d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, c_descr.size);

	c_descr.finish();
}

StokeFlowDescriptor createSubSystem(const StokeFlowDescriptor& f_descr)
{
	StokeFlowDescriptor c_descr;
	c_descr.init(f_descr.Nx / 2, f_descr.Ny / 2);
	c_descr.lap_coeff = f_descr.lap_coeff*(Scalar)2;

	updateSubSystem(c_descr, f_descr);

	return c_descr;
}

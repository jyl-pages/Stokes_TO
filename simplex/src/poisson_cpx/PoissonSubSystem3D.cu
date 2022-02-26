#include "PoissonSubSystem3D.h"
#include "SubSystemUtils3D.cuh"
#include "cuda_runtime.h"
#include "gpuUtils.h"

PoissonDownSample3D::PoissonDownSample3D(PoissonDescriptor<3> *_f_desr, PoissonDescriptor<3> *_c_desr)
{
	f_descr = _f_desr;
	c_descr = _c_desr;

	xdof = f_descr->size;
	ydof = c_descr->size;
}

int PoissonDownSample3D::xDoF()
{
	return xdof;
}

int PoissonDownSample3D::yDoF()
{
	return ydof;
}

void PoissonDownSample3D::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid3D& c_grid = c_descr->grid;
	const grid3D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };

	grid3DSubsystem::CellDownSample(Ap, p, c_grid, f_grid, cell_kernel);

	cwise_mapping_wrapper(Ap, c_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

PoissonUpSample3D::PoissonUpSample3D(PoissonDescriptor<3> *_f_desr, PoissonDescriptor<3> *_c_desr)
{
	f_descr = _f_desr;
	c_descr = _c_desr;

	xdof = c_descr->size;
	ydof = f_descr->size;

}

int PoissonUpSample3D::xDoF()
{
	return xdof;
}

int PoissonUpSample3D::yDoF()
{
	return ydof;
}

void PoissonUpSample3D::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid3D& c_grid = c_descr->grid;
	const grid3D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar v, int x, int y, int z)->Scalar { return v; };

	grid3DSubsystem::CellUpSample(Ap, p, f_grid, c_grid, cell_kernel);

	cwise_mapping_wrapper(Ap, f_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

void updateSubSystem(PoissonDescriptor<3>& c_descr, const PoissonDescriptor<3>& f_descr)
{
	const grid3D& c_grid = c_descr.grid;
	const grid3D& f_grid = f_descr.grid;

	auto fixed_cell_kernel = [=] __device__(bool *v)->bool { bool sum = false; for (int i = 0; i < 8; i++) sum |= v[i]; return sum; };

	cudaMemset(c_descr.d_fixed, 0, sizeof(bool)*c_grid.cell_size());

	grid3DSubsystem::CellDownSample(c_descr.d_fixed, f_descr.d_fixed, c_grid, f_grid, fixed_cell_kernel);

	auto vol_face_x_kernel = [=] __device__(Scalar *v)->Scalar { return v[0b000] + v[0b010] + v[0b100] + v[0b110]; };
	auto vol_face_y_kernel = [=] __device__(Scalar *v)->Scalar { return v[0b000] + v[0b001] + v[0b100] + v[0b101]; };
	auto vol_face_z_kernel = [=] __device__(Scalar *v)->Scalar { return v[0b000] + v[0b001] + v[0b010] + v[0b011]; };

	cudaMemset(c_descr.d_vol, 0, sizeof(Scalar)*(c_grid.face_size(0) + c_grid.face_size(1) + c_grid.face_size(2)));

	grid3DSubsystem::FaceXDownSample(c_descr.d_vol, f_descr.d_vol, c_grid, f_grid, vol_face_x_kernel, (Scalar)0);
	grid3DSubsystem::FaceYDownSample(c_descr.d_vol + c_grid.face_size(0), f_descr.d_vol + f_grid.face_size(0), c_grid, f_grid, vol_face_y_kernel, (Scalar)0);
	grid3DSubsystem::FaceZDownSample\
		(c_descr.d_vol + c_grid.face_size(0) + c_grid.face_size(1), f_descr.d_vol + f_grid.face_size(0) + f_grid.face_size(1), c_grid, f_grid, vol_face_z_kernel, (Scalar)0);

	//cwise_mapping_wrapper(c_descr.d_vol, c_descr.d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, c_descr.size);
}

PoissonDescriptor<3> createSubSystem(const PoissonDescriptor<3>& f_descr)
{
	PoissonDescriptor<3> c_descr;
	c_descr.init(f_descr.N / 2);

	updateSubSystem(c_descr, f_descr);

	return c_descr;
}

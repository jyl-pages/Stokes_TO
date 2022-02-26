#include "StokeFlowSubSystem3D.h"
#include "SubSystemUtils3D.cuh"
#include "cuda_runtime.h"
#include "gpuUtils.h"

void StokeFlowDownSample3D::init(StokeFlowDescriptor3D *_f_desr, StokeFlowDescriptor3D *_c_desr)
{
	f_descr = _f_desr;
	c_descr = _c_desr;

	f_Nx = f_descr->Nx;
	f_Ny = f_descr->Ny;
	f_Nz = f_descr->Nz;
	c_Nx = c_descr->Nx;
	c_Ny = c_descr->Ny;
	c_Nz = c_descr->Nz;

	xdof = f_descr->size;
	ydof = c_descr->size;
}

int StokeFlowDownSample3D::xDoF()
{
	return xdof;
}

int StokeFlowDownSample3D::yDoF()
{
	return ydof;
}

void StokeFlowDownSample3D::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid3D& c_grid = c_descr->grid;
	const grid3D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto face_x_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto face_y_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto face_z_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };

	grid3DSubsystem::CellDownSample(Ap, p, c_grid, f_grid, cell_kernel);
	grid3DSubsystem::FaceXDownSample(Ap + c_grid.cell_size(), p + f_grid.cell_size(), c_grid, f_grid, face_x_kernel, (Scalar)0);
	grid3DSubsystem::FaceYDownSample(Ap + c_grid.cell_size()+c_grid.face_size(0), p + f_grid.cell_size()+f_grid.face_size(0),\
		c_grid, f_grid, face_y_kernel, (Scalar)0);
	grid3DSubsystem::FaceZDownSample(Ap + c_grid.cell_size() + c_grid.face_size(0) + c_grid.face_size(1), \
		p + f_grid.cell_size() + f_grid.face_size(0) + f_grid.face_size(1),\
		c_grid, f_grid, face_z_kernel, (Scalar)0);

	cwise_mapping_wrapper(Ap, c_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

void StokeFlowUpSample3D::init(StokeFlowDescriptor3D *_f_desr, StokeFlowDescriptor3D *_c_desr)
{
	f_descr = _f_desr;
	c_descr = _c_desr;

	f_Nx = f_descr->Nx;
	f_Ny = f_descr->Ny;
	f_Nz = f_descr->Nz;
	c_Nx = c_descr->Nx;
	c_Ny = c_descr->Ny;
	c_Nz = c_descr->Nz;

	xdof = c_descr->size;
	ydof = f_descr->size;

}

int StokeFlowUpSample3D::xDoF()
{
	return xdof;
}

int StokeFlowUpSample3D::yDoF()
{
	return ydof;
}

void StokeFlowUpSample3D::applyMapping(Scalar *Ap, Scalar *p)
{
	const grid3D& c_grid = c_descr->grid;
	const grid3D& f_grid = f_descr->grid;
	cudaMemset(Ap, 0, sizeof(Scalar)*ydof);

	auto cell_kernel = [=] __device__(Scalar v, int x, int y, int z)->Scalar { return v; };
	auto face_x_kernel = [=] __device__(Scalar v, int x, int y, int z)->Scalar { return v; };
	auto face_y_kernel = [=] __device__(Scalar v, int x, int y,int z)->Scalar { return v; };
	auto face_z_kernel = [=] __device__(Scalar v, int x, int y, int z)->Scalar { return v; };

	grid3DSubsystem::CellUpSample(Ap, p, f_grid, c_grid, cell_kernel);
	grid3DSubsystem::FaceXUpSample(Ap + f_grid.cell_size(), p + c_grid.cell_size(), f_grid, c_grid, face_x_kernel, (Scalar)0);
	grid3DSubsystem::FaceYUpSample(Ap + f_grid.cell_size() + f_grid.face_size(0), p + c_grid.cell_size() + c_grid.face_size(0), \
		f_grid, c_grid, face_y_kernel, (Scalar)0);
	grid3DSubsystem::FaceZUpSample(Ap + f_grid.cell_size() + f_grid.face_size(0) + f_grid.face_size(1), \
		p + c_grid.cell_size() + c_grid.face_size(0) + c_grid.face_size(1), \
		f_grid, c_grid, face_z_kernel, (Scalar)0);

	cwise_mapping_wrapper(Ap, f_descr->d_fixed, [=] __device__(Scalar& v, bool fixed) { if (fixed) v = (Scalar)0; }, ydof);
}

void updateSubSystem(StokeFlowDescriptor3D& c_descr, const StokeFlowDescriptor3D& f_descr)
{
	const grid3D& c_grid = c_descr.grid;
	const grid3D& f_grid = f_descr.grid;

	auto fixed_cell_kernel = [=] __device__(bool *v)->bool { bool sum = false; for (int i = 0; i < 8; i++) sum |= v[i]; return sum; };
	auto fixed_face_x_kernel = [=] __device__(bool *v)->bool { bool sum = false; for (int i = 0; i < 8; i++) sum |= v[i]; return sum; };
	auto fixed_face_y_kernel = [=] __device__(bool *v)->bool { bool sum = false; for (int i = 0; i < 8; i++) sum |= v[i]; return sum; };
	auto fixed_face_z_kernel = [=] __device__(bool *v)->bool { bool sum = false; for (int i = 0; i < 8; i++) sum |= v[i]; return sum; };

	cudaMemset(c_descr.d_fixed, 0, sizeof(bool)*c_descr.size);

	grid3DSubsystem::CellDownSample(c_descr.d_fixed, f_descr.d_fixed, c_grid, f_grid, fixed_cell_kernel);
	grid3DSubsystem::FaceXDownSample(c_descr.d_fixed + c_grid.cell_size(), f_descr.d_fixed + f_grid.cell_size(), c_grid, f_grid, fixed_face_x_kernel, false);
	grid3DSubsystem::FaceYDownSample(c_descr.d_fixed + c_grid.cell_size() + c_grid.face_size(0), f_descr.d_fixed + f_grid.cell_size() + f_grid.face_size(0), \
		c_grid, f_grid, fixed_face_y_kernel, false);
	grid3DSubsystem::FaceZDownSample(c_descr.d_fixed + c_grid.cell_size() + c_grid.face_size(0) + c_grid.face_size(1), \
		f_descr.d_fixed + f_grid.cell_size() + f_grid.face_size(0) + f_grid.face_size(1), \
		c_grid, f_grid, fixed_face_z_kernel, false);

	auto vol_cell_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto vol_face_x_kernel = [=] __device__(Scalar *v)->Scalar { return v[0b000] + v[0b010] + v[0b100] + v[0b110]; };
	auto vol_face_y_kernel = [=] __device__(Scalar *v)->Scalar { return v[0b000] + v[0b001] + v[0b100] + v[0b101]; };
	auto vol_face_z_kernel = [=] __device__(Scalar *v)->Scalar { return v[0b000] + v[0b001] + v[0b010] + v[0b011]; };

	cudaMemset(c_descr.d_vol, 0, sizeof(Scalar)*c_descr.size);

	grid3DSubsystem::CellDownSample(c_descr.d_vol, f_descr.d_vol, c_grid, f_grid, vol_cell_kernel);
	grid3DSubsystem::FaceXDownSample(c_descr.d_vol + c_grid.cell_size(), f_descr.d_vol + f_grid.cell_size(), c_grid, f_grid, vol_face_x_kernel, (Scalar)0);
	grid3DSubsystem::FaceYDownSample(c_descr.d_vol + c_grid.cell_size() + c_grid.face_size(0), f_descr.d_vol + f_grid.cell_size() + f_grid.face_size(0), \
		c_grid, f_grid, vol_face_y_kernel, (Scalar)0);
	grid3DSubsystem::FaceZDownSample(c_descr.d_vol + c_grid.cell_size() + c_grid.face_size(0) + c_grid.face_size(1), \
		f_descr.d_vol + f_grid.cell_size() + f_grid.face_size(0) + f_grid.face_size(1), \
		c_grid, f_grid, vol_face_z_kernel, (Scalar)0);

	auto penalty_cell_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto penalty_face_x_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto penalty_face_y_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };
	auto penalty_face_z_kernel = [=] __device__(Scalar *v)->Scalar { Scalar sum(0); for (int i = 0; i < 8; i++) sum += v[i]; return sum; };

	cudaMemset(c_descr.d_penalty, 0, sizeof(Scalar)*c_descr.size);

	grid3DSubsystem::CellDownSample(c_descr.d_penalty, f_descr.d_penalty, c_grid, f_grid, penalty_cell_kernel);
	grid3DSubsystem::FaceXDownSample(c_descr.d_penalty + c_grid.cell_size(), f_descr.d_penalty + f_grid.cell_size(), c_grid, f_grid, penalty_face_x_kernel, (Scalar)0);
	grid3DSubsystem::FaceYDownSample(c_descr.d_penalty + c_grid.cell_size() + c_grid.face_size(0), f_descr.d_penalty + f_grid.cell_size() + f_grid.face_size(0), \
		c_grid, f_grid, penalty_face_y_kernel, (Scalar)0);
	grid3DSubsystem::FaceZDownSample(c_descr.d_penalty + c_grid.cell_size() + c_grid.face_size(0) + c_grid.face_size(1), \
		f_descr.d_penalty + f_grid.cell_size() + f_grid.face_size(0) + f_grid.face_size(1), \
		c_grid, f_grid, penalty_face_z_kernel, (Scalar)0);

	c_descr.finish();
}

StokeFlowDescriptor3D createSubSystem(const StokeFlowDescriptor3D& f_descr)
{
	StokeFlowDescriptor3D c_descr;
	c_descr.init(f_descr.Nx / 2, f_descr.Ny / 2, f_descr.Nz / 2);
	//c_descr.lap_coeff = f_descr.lap_coeff*(Scalar)2;

	updateSubSystem(c_descr, f_descr);

	return c_descr;
}

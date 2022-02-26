#include "FluidEnergy3D.h"
#include "SimplexSupport.h"
#include "cuda_runtime_api.h"

void FluidEnergy3D::init(int _grid_size)
{
	grid_size = _grid_size;
	Vector3i sz = Vector3i(grid_size, grid_size, grid_size);

	grid.Initialize(sz, 1.);
	mac_grid.Initialize(grid);

	descr.init(grid_size, grid_size, grid_size);
	grad.init(&descr);

	h_x = (Scalar*)malloc(sizeof(Scalar)*descr.size);
	cudaMalloc(&d_x, sizeof(Scalar)*descr.size);

	h_Ax = (Scalar*)malloc(sizeof(Scalar)*descr.size);
	cudaMalloc(&d_Ax, sizeof(Scalar)*descr.size);

	//h_penalty = (Scalar*)malloc(sizeof(Scalar)*descr.size);
	//cudaMalloc(&d_penalty, sizeof(Scalar)*descr.size);

	vel_grad.Resize(sz);
	p_grad.Resize(sz);

	cell_penalty_grad.Resize(sz);
	face_penalty_grad.Resize(sz);

}

void FluidEnergy3D::init_boundary(const Field<Scalar, 3>& cell_vol, const Field<int, 3>& cell_fixed, \
	const FaceField<Scalar, 3>& face_vol, const FaceField<int, 3>& face_fixed)
{
	int size = descr.size;
	bool *fixed = new bool[size];
	Scalar *vol = new Scalar[size];
	memset(fixed, 0, sizeof(bool)*descr.size);
	memset(vol, 0, sizeof(Scalar)*descr.size);

	simplexSupport::simplex2solver(cell_vol, descr.grid, vol);
	simplexSupport::simplex2solver(cell_fixed, descr.grid, fixed, [=](int v)->bool { if (v) return true; else return false; });

	simplexSupport::simplex2solver(face_vol, descr.grid, vol + descr.grid.cell_size());
	simplexSupport::simplex2solver(face_fixed, descr.grid, fixed + descr.grid.cell_size(), [=](int v)->bool { if (v) return true; else return false; });

	descr.setFixed(fixed);
	descr.setVol(vol);

	delete[] fixed;
	delete[] vol;
}

void FluidEnergy3D::update_penalty(const Field<Scalar, 3>& cell_penalty, const FaceField<Scalar, 3>& face_penalty)
{
	int size = descr.size;
	Scalar *penalty = new Scalar[size];
	memset(penalty, 0, sizeof(Scalar)*descr.size);

	simplexSupport::simplex2solver(cell_penalty, descr.grid, penalty);
	simplexSupport::simplex2solver(face_penalty, descr.grid, penalty + descr.grid.cell_size());

	descr.setPenalty(penalty);
	descr.toDevice();
	descr.finish();

	delete[] penalty;
}

void FluidEnergy3D::compute_gradient_x(const FaceField<Scalar, 3>& vel, const Field<Scalar, 3>& p)
{
	simplexSupport::simplex2solver(p, descr.grid, h_x);
	simplexSupport::simplex2solver(vel, descr.grid, h_x + descr.grid.cell_size());
	cudaMemcpy(d_x, h_x, sizeof(Scalar)*descr.size, cudaMemcpyHostToDevice);
	grad.applyMapping(d_Ax, d_x);
	cudaMemcpy(h_Ax, d_Ax, sizeof(Scalar)*descr.size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	simplexSupport::solver2simplex(descr.grid, h_Ax, p_grad);
	simplexSupport::solver2simplex(descr.grid, h_Ax + descr.grid.cell_size(), vel_grad);
}

void FluidEnergy3D::compute_gradient_q(const FaceField<Scalar, 3>& vel, const Field<Scalar, 3>& p)
{
	simplexSupport::simplex2solver(p, descr.grid, h_x);
	simplexSupport::simplex2solver(vel, descr.grid, h_x + descr.grid.cell_size());

	for (int i = 0; i < descr.size; i++) h_Ax[i] = (Scalar).5*h_x[i] * h_x[i];

	simplexSupport::solver2simplex(descr.grid, h_Ax, cell_penalty_grad);
	simplexSupport::solver2simplex(descr.grid, h_Ax + descr.grid.cell_size(), face_penalty_grad);
}

Scalar FluidEnergy3D::compute_energy(const FaceField<Scalar, 3>& vel, const Field<Scalar, 3>& p)
{
	Scalar obj = (Scalar)0;

	simplexSupport::simplex2solver(p, descr.grid, h_x);
	simplexSupport::simplex2solver(vel, descr.grid, h_x + descr.grid.cell_size());

	cudaMemcpy(d_x, h_x, sizeof(Scalar)*descr.size, cudaMemcpyHostToDevice);
	grad.applyMapping(d_Ax, d_x);
	cudaMemcpy(h_Ax, d_Ax, sizeof(Scalar)*descr.size, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	for (int i = 0; i < descr.size; i++) obj += h_Ax[i] * h_x[i];

	return (Scalar).5*obj;
}

#include "StokeFlowDirectSolve.h"
#include "cuda_runtime_api.h"
#include "gpuUtils.h"
#include <cstdlib>

void calcA(Scalar *A, const StokeFlowDescriptor& descr, bool closed)
{
	const grid2D& grid = descr.grid;
	Scalar *vol = descr.h_vol;
	bool *fixed = descr.h_fixed;
	int cell_off = 0;
	int face_x_off = cell_off + grid.cell_size();
	int face_y_off = face_x_off + grid.face_size(0);
	int size = descr.size;
	Scalar lap_coef = descr.lap_coeff;

	memset(A, 0, sizeof(Scalar)*size*size);

	for (int y = 0; y < grid.Ny; y++) for (int x = 0; x < grid.Nx; x++)
	{
		int cell_ind = cell_off + grid.cell_ind(x, y);
		int face_ind;
		Scalar temp;

		face_ind = face_x_off + grid.face_ind(x, y, 0);
		temp = vol[face_ind];
		A[cell_ind + face_ind*size] = -temp;
		A[face_ind + cell_ind*size] = -temp;


		face_ind = face_x_off + grid.face_ind(x + 1, y, 0);
		temp = vol[face_ind];
		A[cell_ind + face_ind * size] = temp;
		A[face_ind + cell_ind * size] = temp;


		face_ind = face_y_off + grid.face_ind(x, y, 1);
		temp = vol[face_ind];
		A[cell_ind + face_ind * size] = -temp;
		A[face_ind + cell_ind * size] = -temp;


		face_ind = face_y_off + grid.face_ind(x, y + 1, 1);
		temp = vol[face_ind];
		A[cell_ind + face_ind * size] = temp;
		A[face_ind + cell_ind * size] = temp;

	}

	if (closed)
	{
		for (int j = 0; j < grid.cell_size(); j++) for (int i = 0; i < grid.cell_size(); i++)
			A[i + j * size] = (Scalar)1;
	}

	for (int x = 0; x <= grid.Nx; x++) for (int y = 0; y < grid.Ny; y++)
	{
		int face_ind = face_x_off + grid.face_ind(x, y, 0);
		int face_ind_2;

		A[face_ind + face_ind * size] = (Scalar)4 * lap_coef;

		if (x != 0)
		{
			face_ind_2 = face_x_off + grid.face_ind(x - 1, y, 0);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}

		if (x != grid.Nx)
		{
			face_ind_2 = face_x_off + grid.face_ind(x + 1, y, 0);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}

		if (y != 0)
		{
			face_ind_2 = face_x_off + grid.face_ind(x, y - 1, 0);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}

		if (y != grid.Ny - 1)
		{
			face_ind_2 = face_x_off + grid.face_ind(x, y + 1, 0);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}
	}
	
	for (int y = 0; y <= grid.Ny; y++) for (int x = 0; x < grid.Nx; x++)
	{
		int face_ind = face_y_off + grid.face_ind(x, y, 1);
		int face_ind_2;

		A[face_ind + face_ind * size] = (Scalar)4 * lap_coef;

		if (y != 0)
		{
			face_ind_2 = face_y_off + grid.face_ind(x, y - 1, 1);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}

		if (y != grid.Ny)
		{
			face_ind_2 = face_y_off + grid.face_ind(x, y + 1, 1);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}

		if (x != 0)
		{
			face_ind_2 = face_y_off + grid.face_ind(x - 1, y, 1);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}

		if (x != grid.Nx - 1)
		{
			face_ind_2 = face_y_off + grid.face_ind(x + 1, y, 1);
			A[face_ind_2 + face_ind * size] = (Scalar)-1 * lap_coef;
		}
	}

	for (int i = 0; i < size; i++)
		A[i + i * size] += descr.h_penalty[i];

	for (int i = 0; i < size; i++)
		if (fixed[i])
			for (int j = 0; j < size; j++) A[i + j * size] = A[j + i * size] = (Scalar)0;

	for (int i = 0; i < size; i++)
		if (fixed[i]) A[i + i * size] = (Scalar)1;
}

void StokeFlowDirectSolve::init(StokeFlowDescriptor *_descr)
{
	descr = _descr;
	Nx = descr->Nx;
	Ny = descr->Ny;
	dof = descr->size;
	h_A = (Scalar*)malloc(sizeof(Scalar)*dof*dof);
	cudaMalloc(&d_A, sizeof(Scalar)*dof*dof);
	cudaMalloc(&d_piv, sizeof(int)*dof);
	cudaMalloc(&d_info, sizeof(int));
	cusolverDnCreate(&cusolverDnHandle);
}

void StokeFlowDirectSolve::finish()
{
	descr->toHost();
	calcA(h_A, *descr, closed);
	//{
	//	FILE *f = fopen("rua.txt", "w");
	//	fprintf(f, "%d\n%d\n%d\n", descr->size, descr->size, descr->size*descr->size);
	//	for (int i = 0; i < descr->size; i++)
	//		for (int j = 0; j < descr->size; j++)
	//			fprintf(f, "%d %d %f\n", i, j, h_A[i + j * descr->size]);
	//}
	cudaMemcpy(d_A, h_A, sizeof(Scalar)*dof*dof, cudaMemcpyHostToDevice);
	int buffer_size;
	DnSolve_buffer_lu(cusolverDnHandle, dof, dof, d_A, dof, &buffer_size);
	cudaMalloc(&d_buffer, sizeof(Scalar)*buffer_size);
	DnSolve_factorize_lu(cusolverDnHandle, dof, dof, d_A, dof, d_buffer, d_piv, d_info);
	cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
	printf("LU factorize info:%d\n", h_info);
}

int StokeFlowDirectSolve::xDoF()
{
	return dof;
}

int StokeFlowDirectSolve::yDoF()
{
	return dof;
}

void StokeFlowDirectSolve::applyMapping(Scalar *Ap, Scalar *p)
{
	cudaMemcpy(Ap, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
	DnSolve_solve_lu(cusolverDnHandle, CUBLAS_OP_N, dof, 1, d_A, dof, d_piv, Ap, dof, d_info);
}
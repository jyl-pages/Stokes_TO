#include "PoissonLike.h"
#include "gpuUtils.h"
#include "PartialMapping.h"
#include "DirectSolveBlockInverse.cuh"
#include "cusolverDn.h"
#include <assert.h>

namespace PoissonLike
{
	class WhiteMask2D
	{
	public:
		grid2D grid;
		__host__ __device__ bool operator[](int i) const
		{
			int x, y;
			grid.ind_cell(i, x, y);
			return (x^y) & 0b1;
		}
	};

	class BlackMask2D
	{
	public:
		grid2D grid;
		__host__ __device__ bool operator[](int i) const
		{
			int x, y;
			grid.ind_cell(i, x, y);
			return (x^y ^ 0b1) & 0b1;

		}
	};

	class DiagMapping2D :public LinearMapping
	{
	public:
		int dof;
		Scalar *d_temp;

		LinearMapping *mapping1;
		LinearMapping *mapping2;

	public:
		DiagMapping2D(const grid2D& _grid, LinearMapping *_mapping)
		{
			assert(_mapping->xDoF() == _mapping->yDoF());
			dof = _mapping->xDoF();

			cudaMalloc(&d_temp, sizeof(Scalar)*dof);

			WhiteMask2D *wm = new WhiteMask2D();
			wm->grid = _grid;
			mapping1 = new PartialMapping(wm, _mapping, wm);

			BlackMask2D *bm = new BlackMask2D();
			bm->grid = _grid;
			mapping2 = new PartialMapping(bm, _mapping, bm);

		}

		int xDoF() override
		{
			return dof;

		}

		int yDoF() override
		{
			return dof;
		}

		void applyMapping(Scalar *Ap, Scalar *p) override
		{
			auto add = [=] __device__(Scalar& tv1, Scalar tv2) { tv1 += tv2; };
			mapping1->applyMapping(Ap, p);
			mapping2->applyMapping(d_temp, p);
			cwise_mapping_wrapper(Ap, d_temp, add, dof);
		}
	};

	LinearMapping *GenerateJacobiPreconditioner(const grid2D & grid, LinearMapping * mapping)
	{
		DiagMapping2D *diag = new DiagMapping2D(grid, mapping);
		return new DirectSolveBlockInverse<1>(diag);
	}

	class DirectSolve2D :public LinearMapping
	{
	public:
		grid2D grid;
		LinearMapping *mapping;
		bool closed;
		int dof;
		Scalar *h_A, *d_A;
		Scalar *h_Ap, *d_Ap, *h_p, *d_p;
		int *d_piv;
		Scalar *d_buffer;
		int h_info, *d_info;
		cusolverDnHandle_t cusolverDnHandle;

	public:
		DirectSolve2D(const grid2D& _grid, LinearMapping* _mapping, bool _closed)
		{
			grid = _grid;
			mapping = _mapping;
			closed = _closed;
			assert(mapping->xDoF() == mapping->yDoF());
			dof = mapping->xDoF();
			h_A = (Scalar*)malloc(sizeof(Scalar)*dof*dof);
			h_Ap = (Scalar*)malloc(sizeof(Scalar)*dof);
			h_p = (Scalar*)malloc(sizeof(Scalar)*dof);
			cudaMalloc(&d_A, sizeof(Scalar)*dof*dof);
			cudaMalloc(&d_Ap, sizeof(Scalar)*dof);
			cudaMalloc(&d_p, sizeof(Scalar)*dof);
			cudaMalloc(&d_piv, sizeof(int)*dof);
			cudaMalloc(&d_info, sizeof(int));
			cusolverDnCreate(&cusolverDnHandle);

			int buffer_size;
			DnSolve_buffer_lu(cusolverDnHandle, dof, dof, d_A, dof, &buffer_size);
			cudaMalloc(&d_buffer, sizeof(Scalar)*buffer_size);

		}

		int xDoF() override
		{
			return dof;
		}

		int yDoF() override
		{
			return dof;
		}

		void applyMapping(Scalar *Ap, Scalar *p) override
		{
			cudaMemcpy(Ap, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
			DnSolve_solve_lu(cusolverDnHandle, CUBLAS_OP_N, dof, 1, d_A, dof, d_piv, Ap, dof, d_info);
		}

		void update() override
		{
			memset(h_A, 0, sizeof(Scalar)*dof*dof);
			for (int flag = 0; flag < 5; flag++)
			{
				for (int i = 0; i < dof; i++)
				{
					int x, y;
					grid.ind_cell(i, x, y);
					h_p[i] = ((x + y * 2 + flag) % 5 == 0) ? (Scalar)1 : (Scalar)0;
				}
				cudaMemcpy(d_p, h_p, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
				mapping->applyMapping(d_Ap, d_p);
				cudaMemcpy(h_Ap, d_Ap, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
				cudaDeviceSynchronize();
				for (int i = 0; i < dof; i++)
				{
					int x, y;
					grid.ind_cell(i, x, y);
					int x0 = x, y0 = y;
					int dflag = (x + y * 2 + flag) % 5;
					switch (dflag)
					{
					case 0:
						break;
					case 1:
						x0--;
						break;
					case 2:
						y0--;
						break;
					case 3:
						y0++;
						break;
					case 4:
						x0++;
						break;
					default:
						break;
					}
					if (x0 < grid.Nx && x0 >= 0 && y0 < grid.Ny && y0 >= 0)
					{
						int ni = grid.cell_ind(x0, y0);
						h_A[i + ni * dof] = h_Ap[i];
					}
				}
			}

			if (closed)
				for (int i = 0; i < dof*dof; i++) h_A[i] += (Scalar)1;

			cudaMemcpy(d_A, h_A, sizeof(Scalar)*dof*dof, cudaMemcpyHostToDevice);
			DnSolve_factorize_lu(cusolverDnHandle, dof, dof, d_A, dof, d_buffer, d_piv, d_info);
			cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			printf("LU factorize info:%d\n", h_info);
		}

	};

	LinearMapping *GenerateDirectSolve(const grid2D & grid, LinearMapping * mapping, bool closed)
	{
		return new DirectSolve2D(grid, mapping, closed);
	}
}

namespace PoissonLike
{
	class WhiteMask3D
	{
	public:
		grid3D grid;
		__host__ __device__ bool operator[](int i) const
		{
			int x, y, z;
			grid.ind_cell(i, x, y, z);
			return (x^y^z) & 0b1;
		}
	};

	class BlackMask3D
	{
	public:
		grid3D grid;
		__host__ __device__ bool operator[](int i) const
		{
			int x, y,z;
			grid.ind_cell(i, x, y, z);
			return (x^y^z ^ 0b1) & 0b1;

		}
	};

	class DiagMapping3D :public LinearMapping
	{
	public:
		int dof;
		Scalar *d_temp;

		LinearMapping *mapping1;
		LinearMapping *mapping2;

	public:
		DiagMapping3D(const grid3D& _grid, LinearMapping *_mapping)
		{
			assert(_mapping->xDoF() == _mapping->yDoF());
			dof = _mapping->xDoF();

			cudaMalloc(&d_temp, sizeof(Scalar)*dof);

			WhiteMask3D *wm = new WhiteMask3D();
			wm->grid = _grid;
			mapping1 = new PartialMapping(wm, _mapping, wm);

			BlackMask3D *bm = new BlackMask3D();
			bm->grid = _grid;
			mapping2 = new PartialMapping(bm, _mapping, bm);

		}

		int xDoF() override
		{
			return dof;

		}

		int yDoF() override
		{
			return dof;
		}

		void applyMapping(Scalar *Ap, Scalar *p) override
		{
			auto add = [=] __device__(Scalar& tv1, Scalar tv2) { tv1 += tv2; };
			mapping1->applyMapping(Ap, p);
			mapping2->applyMapping(d_temp, p);
			cwise_mapping_wrapper(Ap, d_temp, add, dof);
		}
	};

	LinearMapping *GenerateJacobiPreconditioner(const grid3D & grid, LinearMapping * mapping)
	{
		DiagMapping3D *diag = new DiagMapping3D(grid, mapping);
		return new DirectSolveBlockInverse<1>(diag);
	}

	class DirectSolve3D :public LinearMapping
	{
	public:
		grid3D grid;
		LinearMapping *mapping;
		bool closed;
		int dof;
		Scalar *h_A, *d_A;
		Scalar *h_Ap, *d_Ap, *h_p, *d_p;
		int *d_piv;
		Scalar *d_buffer;
		int h_info, *d_info;
		cusolverDnHandle_t cusolverDnHandle;

	public:
		DirectSolve3D(const grid3D& _grid, LinearMapping* _mapping, bool _closed)
		{
			grid = _grid;
			mapping = _mapping;
			closed = _closed;
			assert(mapping->xDoF() == mapping->yDoF());
			dof = mapping->xDoF();
			h_A = (Scalar*)malloc(sizeof(Scalar)*dof*dof);
			h_Ap = (Scalar*)malloc(sizeof(Scalar)*dof);
			h_p = (Scalar*)malloc(sizeof(Scalar)*dof);
			cudaMalloc(&d_A, sizeof(Scalar)*dof*dof);
			cudaMalloc(&d_Ap, sizeof(Scalar)*dof);
			cudaMalloc(&d_p, sizeof(Scalar)*dof);
			cudaMalloc(&d_piv, sizeof(int)*dof);
			cudaMalloc(&d_info, sizeof(int));
			cusolverDnCreate(&cusolverDnHandle);

			int buffer_size;
			DnSolve_buffer_lu(cusolverDnHandle, dof, dof, d_A, dof, &buffer_size);
			cudaMalloc(&d_buffer, sizeof(Scalar)*buffer_size);
		}

		int xDoF() override
		{
			return dof;
		}

		int yDoF() override
		{
			return dof;
		}

		void applyMapping(Scalar *Ap, Scalar *p) override
		{
			cudaMemcpy(Ap, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
			DnSolve_solve_lu(cusolverDnHandle, CUBLAS_OP_N, dof, 1, d_A, dof, d_piv, Ap, dof, d_info);
		}

		void update() override
		{
			memset(h_A, 0, sizeof(Scalar)*dof*dof);
			for (int flag = 0; flag < 7; flag++)
			{
				for (int i = 0; i < dof; i++)
				{
					int x, y, z;
					grid.ind_cell(i, x, y, z);
					h_p[i] = ((x + y * 2 + z * 3 + flag) % 7 == 0) ? (Scalar)1 : (Scalar)0;
				}
				cudaMemcpy(d_p, h_p, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
				mapping->applyMapping(d_Ap, d_p);
				cudaMemcpy(h_Ap, d_Ap, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
				cudaDeviceSynchronize();
				for (int i = 0; i < dof; i++)
				{
					int x, y, z;
					grid.ind_cell(i, x, y, z);
					int x0 = x, y0 = y, z0 = z;
					int dflag = (x + y * 2 + z * 3 + flag) % 7;
					switch (dflag)
					{
					case 0:
						break;
					case 1:
						x0--;
						break;
					case 2:
						y0--;
						break;
					case 3:
						z0--;
						break;
					case 4:
						z0++;
						break;
					case 5:
						y0++;
						break;
					case 6:
						x0++;
						break;
					default:
						break;
					}
					if (x0 < grid.Nx && x0 >= 0 && y0 < grid.Ny && y0 >= 0 && z0 < grid.Nz && z0 >= 0)
					{
						int ni = grid.cell_ind(x0, y0, z0);
						h_A[i + ni * dof] = h_Ap[i];
					}
				}
			}

			if (closed)
				for (int i = 0; i < dof*dof; i++) h_A[i] += (Scalar)1;

			cudaMemcpy(d_A, h_A, sizeof(Scalar)*dof*dof, cudaMemcpyHostToDevice);
			DnSolve_factorize_lu(cusolverDnHandle, dof, dof, d_A, dof, d_buffer, d_piv, d_info);
			cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			printf("LU factorize info:%d\n", h_info);
		}
	};

	LinearMapping *GenerateDirectSolve(const grid3D & grid, LinearMapping * mapping, bool closed)
	{
		return new DirectSolve3D(grid, mapping, closed);
	}

}
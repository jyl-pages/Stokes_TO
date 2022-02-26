#include "StokeLike.h"
#include "PartialMapping.h"
#include "DescentSmoother.h"
#include "DirectSolveBlockInverse.cuh"
#include "cusolverDn.h"
#include <assert.h>
#include <vector>

namespace StokeLike
{
	class StokeMask3D
	{
	public:
		grid3D grid;
		int cell_flag;
		int face_flag;
		__host__ __device__ bool operator[](int i) const
		{
			if (i < grid.cell_size())
			{
				int x, y, z;
				grid.ind_cell(i, x, y, z);
				return (x + y + z + cell_flag) %2 == 0;
			}
			else
			{
				int d, x, y, z;
				i -= grid.cell_size();
				for (d = 0; d < 3; d++)
				{
					if (i < grid.face_size(d)) break;
					i -= grid.face_size(d);
				}
				grid.ind_face(i, d, x, y, z);
				//printf("%d %d %d %d %d %d\n", i, d, x, y, z, (x^y^z ^ flag) & 0b1);
				return (x + y + z + face_flag) %2 == 0;
			}
		}
	};

	class VelDiagPressureFullMapping3D :public LinearMapping
	{
	public:
		int dof;
		int flag;
		grid3D grid;
		Scalar *d_temp;

		LinearMapping *mapping1;
		LinearMapping *mapping2;

	public:
		VelDiagPressureFullMapping3D(const grid3D& _grid, LinearMapping *_mapping, int _flag, Scalar *_buffer=nullptr, Scalar *_buffer_partial=nullptr)
		{
			assert(_mapping->xDoF() == _mapping->yDoF());
			dof = _mapping->xDoF();
			flag = _flag;
			grid = _grid;

			if (_buffer == nullptr)
				cudaMalloc(&d_temp, sizeof(Scalar)*dof);
			else
				d_temp = _buffer;

			StokeMask3D *wm = new StokeMask3D();
			wm->grid = grid;
			wm->cell_flag = flag;
			wm->face_flag = 0;
			mapping1 = new PartialMapping(wm, _mapping, wm, _buffer_partial);

			StokeMask3D *bm = new StokeMask3D();
			bm->grid = grid;
			bm->cell_flag = flag;
			bm->face_flag = 1;
			mapping2 = new PartialMapping(bm, _mapping, bm, _buffer_partial);
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

	__global__ void GATHER(Scalar *v1, Scalar *v2, grid3D t_grid, int flag, int N)
	{
		const int i = blockDim.x*blockIdx.x + threadIdx.x;
		if (i >= N) return;

		const int t = i / 7, l = i % 7;
		const int b = (i / 7) >> 5;

		int idx = (t & 0b1) << 1, idy = (t & 0b110) >> 1, idz = (t & 0b11000) >> 3;
		idx += (idy + idz + flag) %2;

		const int nbx = t_grid.Nx >> 2, nby = t_grid.Ny >> 2, nbz = t_grid.Nz >> 2;

		int x = ((b % nbx) << 2) + idx;
		int y = (((b / nbx) % nby) << 2) + idy;
		int z = ((b / nbx / nby) << 2) + idz;

		if (l < 6)
		{
			int d = l >> 1;
			x += (d == 0)*(l & 0b1);
			y += (d == 1)*(l & 0b1);
			z += (d == 2)*(l & 0b1);
			int off = t_grid.cell_size();
			for (int p = 0; p < d; p++)
				off += t_grid.face_size(p);
			//printf("%d %d %d %d %d\n", i, x, y, z, d);
			v1[i] = v2[off + t_grid.face_ind(x, y, z, d)];
		}
		else
		{
			v1[i] = v2[t_grid.cell_ind(x, y, z)];
		}
	}

	__global__ void SHATTER(Scalar *v1, Scalar *v2, grid3D t_grid, int flag, int N)
	{
		const int i = blockDim.x*blockIdx.x + threadIdx.x;
		if (i >= N) return;

		const int t = i / 7, l = i % 7;
		const int b = (i / 7) >> 5;

		int idx = (t & 0b1) << 1, idy = (t & 0b110) >> 1, idz = (t & 0b11000) >> 3;
		idx += (idy + idz + flag) %2;

		const int nbx = t_grid.Nx >> 2, nby = t_grid.Ny >> 2, nbz = t_grid.Nz >> 2;

		int x = ((b % nbx) << 2) + idx;
		int y = (((b / nbx) % nby) << 2) + idy;
		int z = ((b / nbx / nby) << 2) + idz;

		if (l < 6)
		{
			int d = l >> 1;
			x += (d == 0)*(l & 0b1);
			y += (d == 1)*(l & 0b1);
			z += (d == 2)*(l & 0b1);
			int off = t_grid.cell_size();
			for (int p = 0; p < d; p++)
				off += t_grid.face_size(p);
			//printf("%d %d %d %d %d\n", i, x, y, z, d);
			v1[off + t_grid.face_ind(x, y, z, d)] = v2[i];
		}
		else
		{
			v1[t_grid.cell_ind(x, y, z)] = v2[i];
		}
	}

	class VelDiagPressureFullBlockMapping :public LinearMapping
	{
	public:
		int dof;
		int flag;
		grid3D grid;
		Scalar *d_temp_p, *d_temp_Ap;

		VelDiagPressureFullMapping3D *mapping;

	public:
		VelDiagPressureFullBlockMapping(VelDiagPressureFullMapping3D *_mapping, Scalar *_buffer_p = nullptr, Scalar *_buffer_Ap = nullptr)
		{
			grid = _mapping->grid;
			dof = grid.cell_size() / 2 * 7;
			flag = _mapping->flag;
			mapping = _mapping;

			if (_buffer_p == nullptr)
				cudaMalloc(&d_temp_p, sizeof(Scalar)*mapping->xDoF());
			else
				d_temp_p = _buffer_p;
			if (_buffer_Ap == nullptr)
				cudaMalloc(&d_temp_Ap, sizeof(Scalar)*mapping->xDoF());
			else
				d_temp_Ap = _buffer_Ap;
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
			cudaMemset(d_temp_p, 0, sizeof(Scalar)*mapping->xDoF());

			SHATTER << <(dof + 63) / 64, 64 >> > (d_temp_p, p, grid, flag, dof);

			mapping->applyMapping(d_temp_Ap, d_temp_p);

			GATHER << <(dof + 63) / 64, 64 >> > (Ap, d_temp_Ap, grid, flag, dof);
		}
	};

	class VelDiagPressureFullSmoothing :public LinearMapping
	{
	public:
		int dof;
		int flag;
		grid3D grid;
		Scalar *d_temp_p, *d_temp_Ap;

		LinearMapping *block_inv;
		int inv_dof;

	public:
		VelDiagPressureFullSmoothing(const grid3D& _grid, LinearMapping *_mapping, int _flag, Scalar *_buffer_p = nullptr, Scalar *_buffer_Ap = nullptr)
		{
			grid = _grid;
			dof = grid.cell_size()+grid.face_size(0)+ grid.face_size(1)+ grid.face_size(2);
			flag = _flag;
			block_inv = _mapping;
			assert(block_inv->xDoF() == block_inv->yDoF());
			inv_dof = block_inv->xDoF();

			if (_buffer_p == nullptr)
				cudaMalloc(&d_temp_p, sizeof(Scalar)*block_inv->xDoF());
			else
				d_temp_p = _buffer_p;
			if (_buffer_Ap == nullptr)
				cudaMalloc(&d_temp_Ap, sizeof(Scalar)*block_inv->xDoF());
			else
				d_temp_Ap = _buffer_Ap;
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
			cudaMemset(Ap, 0, sizeof(Scalar)*dof);

			GATHER << <(inv_dof + 63) / 64, 64 >> > (d_temp_p, p, grid, flag, inv_dof);

			block_inv->applyMapping(d_temp_Ap, d_temp_p);

			SHATTER << <(inv_dof + 63) / 64, 64 >> > (Ap, d_temp_Ap, grid, flag, inv_dof);
		}

		void update() override
		{
			block_inv->update();
		}
	};


	class CellMask3D
	{
	public:
		grid3D grid;
		int flag;
		__host__ __device__ bool operator[](int i) const
		{
			if (i < grid.cell_size())
			{
				int x, y, z;
				grid.ind_cell(i, x, y, z);
				return (x + 2 * y + 3 * z + flag) % 7 == 0;
			}
			else
			{
				int d, x, y, z;
				i -= grid.cell_size();
				for (d = 0; d < 3; d++)
				{
					if (i < grid.face_size(d)) break;
					i -= grid.face_size(d);
				}
				grid.ind_face(i, d, x, y, z);
				int t = (x + 2 * y + 3 * z + flag) % 7;
				return t == 0 || t == 6 - d;
			}
		}
	};

	class VelFullPressureFullMapping3D :public LinearMapping
	{
	public:
		int dof;
		LinearMapping *mapping;

	public:
		VelFullPressureFullMapping3D(const grid3D& _grid, LinearMapping *_mapping, int _flag)
		{
			assert(_mapping->xDoF() == _mapping->yDoF());
			dof = _mapping->xDoF();

			CellMask3D *cm = new CellMask3D();
			cm->grid = _grid;
			cm->flag = _flag;
			mapping = new PartialMapping(cm, _mapping, cm);
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
			mapping->applyMapping(Ap, p);
		}
	};

	std::vector<LinearMapping*> GenerateJacobiPreconditioner(const grid3D & grid, LinearMapping * mapping, int iters)
	{
		Scalar *buffer, *buffer_p, *buffer_Ap, *buffer_partial,*buffer_block;
		cudaMalloc(&buffer, sizeof(Scalar)*mapping->xDoF());
		cudaMalloc(&buffer_p, sizeof(Scalar)*mapping->xDoF());
		cudaMalloc(&buffer_Ap, sizeof(Scalar)*mapping->xDoF());
		cudaMalloc(&buffer_partial, sizeof(Scalar)*mapping->xDoF());
		cudaMalloc(&buffer_block, sizeof(Scalar)*LinearMappingToBlockMatrix_BufferSize<7, 7>(mapping));

		LinearMapping *t0 = new DirectSolveBlockInverse<7>\
			(new VelDiagPressureFullBlockMapping(new VelDiagPressureFullMapping3D(grid, mapping, 0, buffer,buffer_partial), buffer_p, buffer_Ap ),buffer_block);
		LinearMapping *t1 = new DirectSolveBlockInverse<7>\
			(new VelDiagPressureFullBlockMapping(new VelDiagPressureFullMapping3D(grid, mapping, 1, buffer, buffer_partial), buffer_p, buffer_Ap ), buffer_block);
		VelDiagPressureFullSmoothing* t2 = new VelDiagPressureFullSmoothing(grid, t0, 0, buffer_p, buffer_Ap);
		VelDiagPressureFullSmoothing* t3 = new VelDiagPressureFullSmoothing(grid, t1, 1, buffer_p, buffer_Ap);
		//return t2;
		return std::vector<LinearMapping*>({ new DescentSmoother({ t2,t3 }, { iters,iters }, mapping),new DescentSmoother({ t3,t2 }, { iters,iters }, mapping) });
		/*if(order)
			return new DescentSmoother({ t2,t3 }, { iters,iters }, mapping);
		else
			return new DescentSmoother({ t3,t2 }, { iters,iters }, mapping);*/

		//LinearMapping* t[7];
		//for (int i = 0; i < 7; i++)
		//	t[i] = new VelFullPressureFullMapping3D(grid, mapping, i);
		//return new VelFullPressureFullMapping3D(grid, mapping, 0);

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

			// cell to others
			for (int flag = 0; flag < 2; flag++)
			{
				memset(h_p, 0, sizeof(Scalar)*dof);
				for (int i = 0; i < grid.cell_size(); i++)
				{
					int x, y, z;
					grid.ind_cell(i, x, y, z);
					h_p[i] = ((x + y + z + flag) % 2 == 0) ? (Scalar)1 : (Scalar)0;
				}

				cudaMemcpy(d_p, h_p, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
				mapping->applyMapping(d_Ap, d_p);
				cudaMemcpy(h_Ap, d_Ap, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
				cudaDeviceSynchronize();

				for (int i = 0; i < dof; i++)
				{
					if (i < grid.cell_size()) continue;
					int temp = i - grid.cell_size();
					int d, x, y, z;
					for (d = 0; d < 3; d++)
					{
						if (temp < grid.face_size(d)) break;
						temp -= grid.face_size(d);
					}
					grid.ind_face(temp, d, x, y, z);
					int x0 = x, y0 = y, z0 = z;
					int dflag = (x + y + z + flag) % 2;
					if (dflag == 1)
					{
						switch (d)
						{
						case 0:
							x0--;
							break;
						case 1:
							y0--;
							break;
						case 2:
							z0--;
							break;
						default:
							break;
						}
					}
					if (x0 < grid.Nx && x0 >= 0 && y0 < grid.Ny && y0 >= 0 && z0 < grid.Nz && z0 >= 0)
					{
						int ni = grid.cell_ind(x0, y0, z0);
						h_A[i + ni * dof] = h_Ap[i];
					}
				}
			}

			// face to others
			for (int flag = 0; flag < 7; flag++)
			{
				int magic[3] = { 2,4,5 };
				{
					memset(h_p, 0, sizeof(Scalar)*dof);
					int off = grid.cell_size();
					for (int d = 0; d < 3; d++)
					{
						for (int i = 0; i < grid.face_size(d); i++)
						{
							int x, y, z;
							grid.ind_face(i, d, x, y, z);
							h_p[i + off] = ((x + y * 2 + z * 3 + magic[d] + flag) % 7 == 0) ? (Scalar)1 : (Scalar)0;
						}
						off += grid.face_size(d);
					}
				}

				cudaMemcpy(d_p, h_p, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
				mapping->applyMapping(d_Ap, d_p);
				cudaMemcpy(h_Ap, d_Ap, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
				cudaDeviceSynchronize();

				// to cell
				for (int i = 0; i < grid.cell_size(); i++)
				{
					int x, y, z;
					grid.ind_cell(i, x, y, z);
					int d0 = -1, x0 = x, y0 = y, z0 = z;
					int dflag = (x + y * 2 + z * 3 + flag) % 7;
					switch (dflag)
					{
					case 0:
						break;
					case 1:
						d0 = 1;
						y0++;
						break;
					case 2:
						d0 = 2;
						break;
					case 3:
						d0 = 1;
						break;
					case 4:
						d0 = 0;
						x0++;
						break;
					case 5:
						d0 = 0;
						break;
					case 6:
						d0 = 2;
						z0++;
						break;
					default:
						break;
					}

					if (d0 != -1)
						if (x0 < (grid.Nx + 4 * (d0 == 0)) && x0 >= 0 && y0 < (grid.Ny + 4 * (d0 == 1)) && y0 >= 0 && z0 < (grid.Nz + 4 * (d0 == 2)) && z0 >= 0)
					{
						int ni = grid.face_ind(x0, y0, z0, d0);
						ni += grid.cell_size();
						for (int dd = 0; dd < d0; dd++)
							ni += grid.face_size(dd);
						h_A[i + ni * dof] = h_Ap[i];
					}
				}
				

				// to face
				for (int i = grid.cell_size(); i < dof; i++)
				{
					int temp = i - grid.cell_size();
					int d, x, y, z;
					for (d = 0; d < 3; d++)
					{
						if (temp < grid.face_size(d)) break;
						temp -= grid.face_size(d);
					}
					grid.ind_face(temp, d, x, y, z);
					int x0 = x, y0 = y, z0 = z;
					int dflag = (x + y * 2 + z * 3 + magic[d] + flag) % 7;
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
					default:
						break;
					}
					if (x0 < (grid.Nx + 4*(d == 0)) && x0 >= 0 && y0 < (grid.Ny + 4*(d == 1)) && y0 >= 0 && z0 < (grid.Nz + 4*(d == 2)) && z0 >= 0)
					{
						int ni = grid.face_ind(x0, y0, z0, d);
						ni += grid.cell_size();
						for (int dd = 0; dd < d; dd++)
							ni += grid.face_size(dd);
						h_A[i + ni * dof] = h_Ap[i];
					}
				}
			}

			if (closed)
				for (int i = 0; i < grid.cell_size(); i++)
					for (int j = 0; j < grid.cell_size(); j++)
						h_A[i + j * dof] += (Scalar)1;

			//{
			//	std::vector<int> rowInd;
			//	std::vector<int> colInd;
			//	std::vector<Scalar> val;

			//	for (int i = 0; i < dof; i++)
			//		for (int j = 0; j < dof; j++)
			//		{
			//			if (fabs(h_A[j*dof + i]) > 1e-3)
			//			{
			//				rowInd.push_back(j);
			//				colInd.push_back(i);
			//				val.push_back(h_A[j*dof + i]);
			//			}
			//		}

			//	FILE *file = fopen("matrix1.txt", "w");

			//	fprintf(file, "%d %d %d\n", dof, dof, val.size());

			//	for (int i = 0; i < val.size(); i++)
			//		fprintf(file, "%d %d %f\n", rowInd[i], colInd[i], val[i]);

			//	fclose(file);
			//}

			cudaMemcpy(d_A, h_A, sizeof(Scalar)*dof*dof, cudaMemcpyHostToDevice);
			DnSolve_factorize_lu(cusolverDnHandle, dof, dof, d_A, dof, d_buffer, d_piv, d_info);
			cudaMemcpy(&h_info, d_info, sizeof(int), cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();
			printf("LU factorize info:%d\n", h_info);
		}
	};


	LinearMapping * GenerateDirectSolve(const grid3D & grid, LinearMapping * mapping, bool closed)
	{
		return new DirectSolve3D(grid, mapping, closed);
	}
}

#include "cuda_runtime.h"
#include "form01mapping.h"

// 8x8 in coarse grid view
template<typename T, typename F>
__global__ void CellDownSampleKernel(T *c_cell, T *f_cell, grid2D c_grid, grid2D f_grid, F f, T default_value = T(0))
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	const int x = bx * 8 + idx;
	const int y = by * 8 + idy;

	__shared__ T shared_f_v[64];
	__shared__ T shared_c_v[64];

	// left lower
	shared_f_v[idy * 8 + idx] = f_cell[f_grid.cell_ind(bx * 2 * 8 + idx, by * 2 * 8 + idy)];
	__syncthreads();
	if (idx < 4 && idy < 4)
	{
		int tmpx = idx, tmpy = idy;
		shared_c_v[idy * 8 + idx] = f(
			shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
			shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	// right lower
	shared_f_v[idy * 8 + idx] = f_cell[f_grid.cell_ind((bx * 2 + 1) * 8 + idx, by * 2 * 8 + idy)];
	__syncthreads();
	if (idx >= 4 && idy < 4)
	{
		int tmpx = idx - 4, tmpy = idy;
		shared_c_v[idy * 8 + idx] = f(
			shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
			shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	// left upper
	shared_f_v[idy * 8 + idx] = f_cell[f_grid.cell_ind(bx * 2 * 8 + idx, (by * 2 + 1) * 8 + idy)];
	__syncthreads();
	if (idx < 4 && idy >= 4)
	{
		int tmpx = idx, tmpy = idy - 4;
		shared_c_v[idy * 8 + idx] = f(
			shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
			shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	// right upper
	shared_f_v[idy * 8 + idx] = f_cell[f_grid.cell_ind((bx * 2 + 1) * 8 + idx, (by * 2 + 1) * 8 + idy)];
	__syncthreads();
	if (idx >= 4 && idy >= 4)
	{
		int tmpx = idx - 4, tmpy = idy - 4;
		shared_c_v[idy * 8 + idx] = f(
			shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
			shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	c_cell[c_grid.cell_ind(x, y)] = shared_c_v[idy * 8 + idx];
}

// 8x8 in coarse grid view
template<typename T, typename F>
__global__ void FaceXDownSampleKernel(T *c_face_x, T *f_face_x, grid2D c_grid, grid2D f_grid, F f, T default_value = T(0))
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ T shared_f_v[64];
	__shared__ T shared_c_v[64];

	// left lower
	shared_f_v[idy * 8 + idx] = (bx == nbx - 1 && idy) ? default_value : f_face_x[f_grid.face_ind(bx * 2 * 8 + idy, by * 2 * 8 + idx, 0)];
	__syncthreads();
	if (idx < 4 && idy < 4)
	{
		int tmpx = idx, tmpy = idy;
		shared_c_v[idx * 8 + idy] = \
			f(shared_f_v[tmpx * 2 * 8 + tmpy * 2], shared_f_v[(tmpx * 2 + 1) * 8 + tmpy * 2], \
				shared_f_v[tmpx * 2 * 8 + (tmpy * 2 + 1)], shared_f_v[(tmpx * 2 + 1) * 8 + (tmpy * 2 + 1)]);
	}
	__syncthreads();

	// right lower
	shared_f_v[idy * 8 + idx] = (bx == nbx - 1) ? default_value : f_face_x[f_grid.face_ind((bx * 2 + 1) * 8 + idy, by * 2 * 8 + idx, 0)];
	__syncthreads();
	if (idx >= 4 && idy < 4)
	{
		int tmpx = idx - 4, tmpy = idy;
		shared_c_v[idx * 8 + idy] = \
			f(shared_f_v[tmpx * 2 * 8 + tmpy * 2], shared_f_v[(tmpx * 2 + 1) * 8 + tmpy * 2], \
				shared_f_v[tmpx * 2 * 8 + (tmpy * 2 + 1)], shared_f_v[(tmpx * 2 + 1) * 8 + (tmpy * 2 + 1)]);
	}
	__syncthreads();

	// left upper
	shared_f_v[idy * 8 + idx] = (bx == nbx - 1 && idy) ? default_value : f_face_x[f_grid.face_ind(bx * 2 * 8 + idy, (by * 2 + 1) * 8 + idx, 0)];
	__syncthreads();
	if (idx < 4 && idy >= 4)
	{
		int tmpx = idx, tmpy = idy - 4;
		shared_c_v[idx * 8 + idy] = \
			f(shared_f_v[tmpx * 2 * 8 + tmpy * 2], shared_f_v[(tmpx * 2 + 1) * 8 + tmpy * 2], \
				shared_f_v[tmpx * 2 * 8 + (tmpy * 2 + 1)], shared_f_v[(tmpx * 2 + 1) * 8 + (tmpy * 2 + 1)]);
	}
	__syncthreads();

	// right upper
	shared_f_v[idy * 8 + idx] = (bx == nbx - 1) ? default_value : f_face_x[f_grid.face_ind((bx * 2 + 1) * 8 + idy, (by * 2 + 1) * 8 + idx, 0)];
	__syncthreads();
	if (idx >= 4 && idy >= 4)
	{
		int tmpx = idx - 4, tmpy = idy - 4;
		shared_c_v[idx * 8 + idy] = \
			f(shared_f_v[tmpx * 2 * 8 + tmpy * 2], shared_f_v[(tmpx * 2 + 1) * 8 + tmpy * 2], \
				shared_f_v[tmpx * 2 * 8 + (tmpy * 2 + 1)], shared_f_v[(tmpx * 2 + 1) * 8 + (tmpy * 2 + 1)]);
	}
	__syncthreads();

	c_face_x[c_grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)] = shared_c_v[idy * 8 + idx];
}

// 8x8 in coarse grid view
template<typename T, typename F>
__global__ void FaceYDownSampleKernel(T *c_face_y, T *f_face_y, grid2D c_grid, grid2D f_grid, F f, T default_value = T(0))
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ T shared_f_v[64];
	__shared__ T shared_c_v[64];

	// left lower
	shared_f_v[idy * 8 + idx] = (by == nby - 1 && idy) ? default_value : f_face_y[f_grid.face_ind(bx * 2 * 8 + idx, by * 2 * 8 + idy, 1)];
	__syncthreads();
	if (idx < 4 && idy < 4)
	{
		int tmpx = idx, tmpy = idy;
		shared_c_v[idy * 8 + idx] = \
			f(shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
				shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	// right lower
	shared_f_v[idy * 8 + idx] = (by == nby - 1 && idy) ? default_value : f_face_y[f_grid.face_ind((bx * 2 + 1) * 8 + idx, by * 2 * 8 + idy, 1)];
	__syncthreads();
	if (idx >= 4 && idy < 4)
	{
		int tmpx = idx - 4, tmpy = idy;
		shared_c_v[idy * 8 + idx] = \
			f(shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
				shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	// left upper
	shared_f_v[idy * 8 + idx] = (by == nby - 1) ? default_value : f_face_y[f_grid.face_ind(bx * 2 * 8 + idx, (by * 2 + 1) * 8 + idy, 1)];
	__syncthreads();
	if (idx < 4 && idy >= 4)
	{
		int tmpx = idx, tmpy = idy - 4;
		shared_c_v[idy * 8 + idx] = \
			f(shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
				shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	// right upper
	shared_f_v[idy * 8 + idx] = (by == nby - 1) ? default_value : f_face_y[f_grid.face_ind((bx * 2 + 1) * 8 + idx, (by * 2 + 1) * 8 + idy, 1)];
	__syncthreads();
	if (idx >= 4 && idy >= 4)
	{
		int tmpx = idx - 4, tmpy = idy - 4;
		shared_c_v[idy * 8 + idx] = \
			f(shared_f_v[tmpy * 2 * 8 + tmpx * 2], shared_f_v[tmpy * 2 * 8 + (tmpx * 2 + 1)], \
				shared_f_v[(tmpy * 2 + 1) * 8 + tmpx * 2], shared_f_v[(tmpy * 2 + 1) * 8 + (tmpx * 2 + 1)]);
	}
	__syncthreads();

	c_face_y[c_grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)] = shared_c_v[idy * 8 + idx];
}

// 16x16 in fine grid view
template<typename T, typename F>
__global__ void CellUpSampleKernel(T *f_cell, T *c_cell, grid2D f_grid, grid2D c_grid, F f, T default_value = T(0))
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	int tmpx = idx & 0b111, tmpy = 2 * idy + (idx >> 3);
	int tmpbx = (tmpy & 0b01000) >> 3, tmpby = (tmpy & 0b10000) >> 4;
	tmpy = tmpy & 0b111;

	__shared__ T shared_c_v[64];

	if (tmpbx == 0 && tmpby == 0) shared_c_v[tmpy * 8 + tmpx] = c_cell[c_grid.cell_ind(bx * 8 + tmpx, by * 8 + tmpy)];
	__syncthreads();

	f_cell[f_grid.cell_ind((bx * 2 + tmpbx) * 8 + tmpx, (by * 2 + tmpby) * 8 + tmpy)] = \
		f(shared_c_v[(tmpby * 4 + (tmpy >> 1)) * 8 + (tmpbx * 4 + (tmpx >> 1))], tmpx & 1, tmpy & 1);
}

// 16x16 in fine grid view
template<typename T, typename F>
__global__ void FaceXUpSampleKernel(T *f_face_x, T *c_face_x, grid2D f_grid, grid2D c_grid, F f, T default_value = T(0))
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	int tmpx = idx & 0b111, tmpy = 2 * idy + (idx >> 3);
	int tmpbx = (tmpy & 0b01000) >> 3, tmpby = (tmpy & 0b10000) >> 4;
	tmpy = tmpy & 0b111;

	__shared__ T shared_c_v[64];

	if (tmpbx == 0 && tmpby == 0) shared_c_v[tmpy * 8 + tmpx] = c_face_x[c_grid.face_ind(bx * 8 + tmpy, by * 8 + tmpx, 0)];
	__syncthreads();

	if (bx != nbx - 1 || tmpbx != 1)
		f_face_x[f_grid.face_ind((bx * 2 + tmpbx) * 8 + tmpx, (by * 2 + tmpby) * 8 + tmpy, 0)] = \
		f(shared_c_v[(tmpbx * 4 + (tmpx >> 1)) * 8 + (tmpby * 4 + (tmpy >> 1))], tmpx & 1, tmpy & 1);
}

// 16x16 in fine grid view
template<typename T, typename F>
__global__ void FaceYUpSampleKernel(T *f_face_y, T *c_face_y, grid2D f_grid, grid2D c_grid, F f, T default_value = T(0))
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	int tmpx = idx & 0b111, tmpy = 2 * idy + (idx >> 3);
	int tmpbx = (tmpy & 0b01000) >> 3, tmpby = (tmpy & 0b10000) >> 4;
	tmpy = tmpy & 0b111;


	__shared__ T shared_c_v[64];

	if (tmpbx == 0 && tmpby == 0) shared_c_v[tmpy * 8 + tmpx] = c_face_y[c_grid.face_ind(bx * 8 + tmpx, by * 8 + tmpy, 1)];
	__syncthreads();

	if (by != nby - 1 || tmpby != 1)
		f_face_y[f_grid.face_ind((bx * 2 + tmpbx) * 8 + tmpx, (by * 2 + tmpby) * 8 + tmpy, 1)] = \
		f(shared_c_v[(tmpby * 4 + (tmpy >> 1)) * 8 + (tmpbx * 4 + (tmpx >> 1))], tmpx & 1, tmpy & 1);
}

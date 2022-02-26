#include "cuda_runtime.h"
#include "grid3D.h"

namespace grid3DSubsystem
{
	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void CellDownSampleKernel(T *c_cell, T *f_cell, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_cell[f_grid.cell_ind((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz)];
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
				};
				shared_c_v[((idz + tz * 2) * 4 + (idy + ty * 2)) * 4 + (idx + tx * 2)] = f(v);
			}
			__syncthreads();
		}

		c_cell[c_grid.cell_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void CellDownSample(T *c_cell, T *f_cell, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		CellDownSampleKernel <T, F> << <dim3(c_grid.Nx / 4, c_grid.Ny / 4, c_grid.Nz / 4), dim3(4, 4, 4) >> > (c_cell, f_cell, c_grid, f_grid, f, default_value);
	}

	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void FaceXDownSampleKernel(T *c_face_x, T *f_face_x, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.face_valid((bx * 2 + tx) * 4 + idz, (by * 2 + ty) * 4 + idx, (bz * 2 + tz) * 4 + idy, 0) ?
				f_face_x[f_grid.face_ind((bx * 2 + tx) * 4 + idz, (by * 2 + ty) * 4 + idx, (bz * 2 + tz) * 4 + idy, 0)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idx * 2 * 4 + idz * 2) * 4 + idy * 2],
					shared_f_v[((idx * 2 + 1) * 4 + idz * 2) * 4 + idy * 2],
					shared_f_v[(idx * 2 * 4 + idz * 2) * 4 + (idy * 2 + 1)],
					shared_f_v[((idx * 2 + 1) * 4 + idz * 2) * 4 + (idy * 2 + 1)],
					shared_f_v[(idx * 2 * 4 + (idz * 2 + 1)) * 4 + idy * 2],
					shared_f_v[((idx * 2 + 1) * 4 + (idz * 2 + 1)) * 4 + idy * 2],
					shared_f_v[(idx * 2 * 4 + (idz * 2 + 1)) * 4 + (idy * 2 + 1)],
					shared_f_v[((idx * 2 + 1) * 4 + (idz * 2 + 1)) * 4 + (idy * 2 + 1)],
				};
				shared_c_v[((idx + tx * 2) * 4 + (idz + tz * 2)) * 4 + (idy + ty * 2)] = f(v);
			}
			__syncthreads();
		}
		c_face_x[c_grid.face_ind(bx * 4 + idz, by * 4 + idx, bz * 4 + idy, 0)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void FaceXDownSample(T *c_face_x, T *f_face_x, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		FaceXDownSampleKernel<T, F> << <dim3(c_grid.Nx / 4 + 1, c_grid.Ny / 4, c_grid.Nz / 4), dim3(4, 4, 4) >> > (c_face_x, f_face_x, c_grid, f_grid, f, default_value);
	}

	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void FaceYDownSampleKernel(T *c_face_y, T *f_face_y, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.face_valid((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idz, (bz * 2 + tz) * 4 + idy, 1) ?
				f_face_y[f_grid.face_ind((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idz, (bz * 2 + tz) * 4 + idy, 1)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idy * 2 * 4 + idz * 2) * 4 + idx * 2],
					shared_f_v[(idy * 2 * 4 + idz * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[((idy * 2 + 1) * 4 + idz * 2) * 4 + idx * 2],
					shared_f_v[((idy * 2 + 1) * 4 + idz * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[(idy * 2 * 4 + (idz * 2 + 1)) * 4 + idx * 2],
					shared_f_v[(idy * 2 * 4 + (idz * 2 + 1)) * 4 + (idx * 2 + 1)],
					shared_f_v[((idy * 2 + 1) * 4 + (idz * 2 + 1)) * 4 + idx * 2],
					shared_f_v[((idy * 2 + 1) * 4 + (idz * 2 + 1)) * 4 + (idx * 2 + 1)],
				};
				shared_c_v[((idy + ty * 2) * 4 + (idz + tz * 2)) * 4 + (idx + tx * 2)] = f(v);
			}
			__syncthreads();
		}
		c_face_y[c_grid.face_ind(bx * 4 + idx, by * 4 + idz, bz * 4 + idy, 1)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void FaceYDownSample(T *c_face_y, T *f_face_y, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		FaceYDownSampleKernel <T, F> << <dim3(c_grid.Nx / 4, c_grid.Ny / 4 + 1, c_grid.Nz / 4), dim3(4, 4, 4) >> > (c_face_y, f_face_y, c_grid, f_grid, f, default_value);
	}


	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void FaceZDownSampleKernel(T *c_face_y, T *f_face_y, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.face_valid((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz, 2) ?
				f_face_y[f_grid.face_ind((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz, 2)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
				};
				shared_c_v[((idz + tz * 2) * 4 + (idy + ty * 2)) * 4 + (idx + tx * 2)] = f(v);
			}
			__syncthreads();
		}
		c_face_y[c_grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 2)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void FaceZDownSample(T *c_face_z, T *f_face_z, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		FaceZDownSampleKernel<T, F> << <dim3(c_grid.Nx / 4, c_grid.Ny / 4, c_grid.Nz / 4 + 1), dim3(4, 4, 4) >> > (c_face_z, f_face_z, c_grid, f_grid, f, default_value);
	}

	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void EdgeXDownSampleKernel(T *c_edge_x, T *f_edge_x, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.edge_valid((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz, 0) ?
				f_edge_x[f_grid.edge_ind((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz, 0)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
				};
				shared_c_v[((idz + tz * 2) * 4 + (idy + ty * 2)) * 4 + (idx + tx * 2)] = f(v);
			}
			__syncthreads();
		}
		c_edge_x[c_grid.edge_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 0)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void EdgeXDownSample(T *c_edge_x, T *f_edge_x, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		EdgeXDownSampleKernel<T, F> << <dim3(c_grid.Nx / 4, c_grid.Ny / 4 + 1, c_grid.Nz / 4 + 1), dim3(4, 4, 4) >> > (c_edge_x, f_edge_x, c_grid, f_grid, f, default_value);
	}

	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void EdgeYDownSampleKernel(T *c_edge_y, T *f_edge_y, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.edge_valid((bx * 2 + tx) * 4 + idy, (by * 2 + ty) * 4 + idx, (bz * 2 + tz) * 4 + idz, 1) ?
				f_edge_y[f_grid.edge_ind((bx * 2 + tx) * 4 + idy, (by * 2 + ty) * 4 + idx, (bz * 2 + tz) * 4 + idz, 1)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idz * 2 * 4 + idx * 2) * 4 + idy * 2],
					shared_f_v[(idz * 2 * 4 + (idx * 2 + 1)) * 4 + idy * 2],
					shared_f_v[(idz * 2 * 4 + idx * 2) * 4 + (idy * 2 + 1)],
					shared_f_v[(idz * 2 * 4 + (idx * 2 + 1)) * 4 + (idy * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + idx * 2) * 4 + idy * 2],
					shared_f_v[((idz * 2 + 1) * 4 + (idx * 2 + 1)) * 4 + idy * 2],
					shared_f_v[((idz * 2 + 1) * 4 + idx * 2) * 4 + (idy * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + (idx * 2 + 1)) * 4 + (idy * 2 + 1)],
				};
				shared_c_v[((idz + tz * 2) * 4 + (idx + tx * 2)) * 4 + (idy + ty * 2)] = f(v);
			}
			__syncthreads();
		}
		c_edge_y[c_grid.edge_ind(bx * 4 + idy, by * 4 + idx, bz * 4 + idz, 1)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void EdgeYDownSample(T *c_edge_y, T *f_edge_y, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		EdgeYDownSampleKernel<T, F> << <dim3(c_grid.Nx / 4 + 1, c_grid.Ny / 4, c_grid.Nz / 4 + 1), dim3(4, 4, 4) >> > (c_edge_y, f_edge_y, c_grid, f_grid, f, default_value);
	}

	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void EdgeZDownSampleKernel(T *c_edge_z, T *f_edge_z, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.edge_valid((bx * 2 + tx) * 4 + idy, (by * 2 + ty) * 4 + idz, (bz * 2 + tz) * 4 + idx, 2) ?
				f_edge_z[f_grid.edge_ind((bx * 2 + tx) * 4 + idy, (by * 2 + ty) * 4 + idz, (bz * 2 + tz) * 4 + idx, 2)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idy * 2 * 4 + idx * 2) * 4 + idz * 2],
					shared_f_v[(idy * 2 * 4 + (idx * 2 + 1)) * 4 + idz * 2],
					shared_f_v[((idy * 2 + 1) * 4 + idx * 2) * 4 + idz * 2],
					shared_f_v[((idy * 2 + 1) * 4 + (idx * 2 + 1)) * 4 + idz * 2],
					shared_f_v[(idy * 2 * 4 + idx * 2) * 4 + (idz * 2 + 1)],
					shared_f_v[(idy * 2 * 4 + (idx * 2 + 1)) * 4 + (idz * 2 + 1)],
					shared_f_v[((idy * 2 + 1) * 4 + idx * 2) * 4 + (idz * 2 + 1)],
					shared_f_v[((idy * 2 + 1) * 4 + (idx * 2 + 1)) * 4 + (idz * 2 + 1)],
				};
				shared_c_v[((idy + ty * 2) * 4 + (idx + tx * 2)) * 4 + (idz + tz * 2)] = f(v);
			}
			__syncthreads();
		}
		c_edge_z[c_grid.edge_ind(bx * 4 + idy, by * 4 + idz, bz * 4 + idx, 2)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void EdgeZDownSample(T *c_edge_z, T *f_edge_z, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		EdgeZDownSampleKernel<T, F> << <dim3(c_grid.Nx / 4 + 1, c_grid.Ny / 4 + 1, c_grid.Nz / 4), dim3(4, 4, 4) >> > (c_edge_z, f_edge_z, c_grid, f_grid, f, default_value);
	}

	// 4x4x4 in coarse grid view
	template<typename T, typename F>
	__global__ void NodeDownSampleKernel(T *c_node, T *f_node, grid3D c_grid, grid3D f_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		__shared__ T shared_f_v[64];
		__shared__ T shared_c_v[64];

		for (int tz = 0; tz < 2; tz++) for (int ty = 0; ty < 2; ty++) for (int tx = 0; tx < 2; tx++)
		{
			shared_f_v[(idz * 4 + idy) * 4 + idx] = f_grid.node_valid((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz) ?
				f_node[f_grid.node_ind((bx * 2 + tx) * 4 + idx, (by * 2 + ty) * 4 + idy, (bz * 2 + tz) * 4 + idz)] : default_value;
			__syncthreads();
			if (idx < 2 && idy < 2 && idz < 2)
			{
				T v[8] = {
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[(idz * 2 * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + idy * 2) * 4 + (idx * 2 + 1)],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + idx * 2],
					shared_f_v[((idz * 2 + 1) * 4 + (idy * 2 + 1)) * 4 + (idx * 2 + 1)],
				};
				shared_c_v[((idz + tz * 2) * 4 + (idy + ty * 2)) * 4 + (idx + tx * 2)] = f(v);
			}
			__syncthreads();
		}

		c_node[c_grid.node_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz)] = shared_c_v[(idz * 4 + idy) * 4 + idx];
	}

	template<typename T, typename F>
	void NodeDownSample(T *c_node, T *f_node, const grid3D& c_grid, const grid3D& f_grid, F f, T default_value = T(0))
	{
		NodeDownSampleKernel <T, F> << <dim3(c_grid.Nx / 4+1, c_grid.Ny / 4+1, c_grid.Nz / 4+1), dim3(4, 4, 4) >> > (c_node, f_node, c_grid, f_grid, f, default_value);
	}


	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void CellUpSampleKernel(T *f_cell, T *c_cell, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_cell[c_grid.cell_ind(bx * 4 + tmpx, by * 4 + tmpy, bz * 4 + tmpz)];
		__syncthreads();

		f_cell[f_grid.cell_ind((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpy, (bz * 2 + tmpbz) * 4 + tmpz)] = \
			f(shared_c_v[((tmpbz * 2 + (tmpz >> 1)) * 4 + (tmpby * 2 + (tmpy >> 1))) * 4 + (tmpbx * 2 + (tmpx >> 1))], tmpx & 0b1, tmpy & 0b1, tmpz & 0b1);
	}

	template<typename T, typename F>
	void CellUpSample(T *f_cell, T *c_cell, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		CellUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8, f_grid.Ny / 8, f_grid.Nz / 8), dim3(8, 8, 8) >> > (f_cell, c_cell, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void FaceXUpSampleKernel(T *f_face_x, T *c_face_x, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_face_x[c_grid.face_ind(bx * 4 + tmpz, by * 4 + tmpx, bz * 4 + tmpy, 0)];
		__syncthreads();

		if (f_grid.face_valid((bx * 2 + tmpbx) * 4 + tmpz, (by * 2 + tmpby) * 4 + tmpx, (bz * 2 + tmpbz) * 4 + tmpy, 0))
			f_face_x[f_grid.face_ind((bx * 2 + tmpbx) * 4 + tmpz, (by * 2 + tmpby) * 4 + tmpx, (bz * 2 + tmpbz) * 4 + tmpy, 0)] = \
			f(shared_c_v[((tmpbx * 2 + (tmpz >> 1)) * 4 + (tmpbz * 2 + (tmpy >> 1))) * 4 + (tmpby * 2 + (tmpx >> 1))], tmpz & 0b1, tmpx & 0b1, tmpy & 0b1);
	}

	template<typename T, typename F>
	void FaceXUpSample(T *f_face_x, T *c_face_x, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		FaceXUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8 + 1, f_grid.Ny / 8, f_grid.Nz / 8), dim3(8, 8, 8) >> > (f_face_x, c_face_x, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void FaceYUpSampleKernel(T *f_face_y, T *c_face_y, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_face_y[c_grid.face_ind(bx * 4 + tmpx, by * 4 + tmpz, bz * 4 + tmpy, 1)];
		__syncthreads();

		if (f_grid.face_valid((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpz, (bz * 2 + tmpbz) * 4 + tmpy, 1))
			f_face_y[f_grid.face_ind((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpz, (bz * 2 + tmpbz) * 4 + tmpy, 1)] = \
			f(shared_c_v[((tmpby * 2 + (tmpz >> 1)) * 4 + (tmpbz * 2 + (tmpy >> 1))) * 4 + (tmpbx * 2 + (tmpx >> 1))], tmpx & 0b1, tmpz & 0b1, tmpy & 0b1);
	}

	template<typename T, typename F>
	void FaceYUpSample(T *f_face_y, T *c_face_y, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		FaceYUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8, f_grid.Ny / 8 + 1, f_grid.Nz / 8), dim3(8, 8, 8) >> > (f_face_y, c_face_y, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void FaceZUpSampleKernel(T *f_face_z, T *c_face_z, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_face_z[c_grid.face_ind(bx * 4 + tmpx, by * 4 + tmpy, bz * 4 + tmpz, 2)];
		__syncthreads();

		if (f_grid.face_valid((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpy, (bz * 2 + tmpbz) * 4 + tmpz, 2))
			f_face_z[f_grid.face_ind((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpy, (bz * 2 + tmpbz) * 4 + tmpz, 2)] = \
			f(shared_c_v[((tmpbz * 2 + (tmpz >> 1)) * 4 + (tmpby * 2 + (tmpy >> 1))) * 4 + (tmpbx * 2 + (tmpx >> 1))], tmpx & 0b1, tmpy & 0b1, tmpz & 0b1);
	}

	template<typename T, typename F>
	void FaceZUpSample(T *f_face_z, T *c_face_z, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		FaceZUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8, f_grid.Ny / 8, f_grid.Nz / 8 + 1), dim3(8, 8, 8) >> > (f_face_z, c_face_z, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void EdgeXUpSampleKernel(T *f_edge_x, T *c_edge_x, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_edge_x[c_grid.edge_ind(bx * 4 + tmpx, by * 4 + tmpy, bz * 4 + tmpz, 0)];
		__syncthreads();

		if (f_grid.edge_valid((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpy, (bz * 2 + tmpbz) * 4 + tmpz, 0))
			f_edge_x[f_grid.edge_ind((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpy, (bz * 2 + tmpbz) * 4 + tmpz, 0)] = \
			f(shared_c_v[((tmpbz * 2 + (tmpz >> 1)) * 4 + (tmpby * 2 + (tmpy >> 1))) * 4 + (tmpbx * 2 + (tmpx >> 1))], tmpx & 0b1, tmpy & 0b1, tmpz & 0b1);
	}

	template<typename T, typename F>
	void EdgeXUpSample(T *f_edge_x, T *c_edge_x, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		EdgeXUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8, f_grid.Ny / 8 + 1, f_grid.Nz / 8 + 1), dim3(8, 8, 8) >> > (f_edge_x, c_edge_x, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void EdgeYUpSampleKernel(T *f_edge_y, T *c_edge_y, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_edge_y[c_grid.edge_ind(bx * 4 + tmpy, by * 4 + tmpx, bz * 4 + tmpz, 1)];
		__syncthreads();

		if (f_grid.edge_valid((bx * 2 + tmpbx) * 4 + tmpy, (by * 2 + tmpby) * 4 + tmpx, (bz * 2 + tmpbz) * 4 + tmpz, 1))
			f_edge_y[f_grid.edge_ind((bx * 2 + tmpbx) * 4 + tmpy, (by * 2 + tmpby) * 4 + tmpx, (bz * 2 + tmpbz) * 4 + tmpz, 1)] = \
			f(shared_c_v[((tmpbz * 2 + (tmpz >> 1)) * 4 + (tmpbx * 2 + (tmpy >> 1))) * 4 + (tmpby * 2 + (tmpx >> 1))], tmpy & 0b1, tmpx & 0b1, tmpz & 0b1);
	}

	template<typename T, typename F>
	void EdgeYUpSample(T *f_edge_y, T *c_edge_y, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		EdgeYUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8 + 1, f_grid.Ny / 8, f_grid.Nz / 8 + 1), dim3(8, 8, 8) >> > (f_edge_y, c_edge_y, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void EdgeZUpSampleKernel(T *f_edge_z, T *c_edge_z, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_edge_z[c_grid.edge_ind(bx * 4 + tmpy, by * 4 + tmpz, bz * 4 + tmpx, 2)];
		__syncthreads();

		if (f_grid.edge_valid((bx * 2 + tmpbx) * 4 + tmpy, (by * 2 + tmpby) * 4 + tmpz, (bz * 2 + tmpbz) * 4 + tmpx, 2))
			f_edge_z[f_grid.edge_ind((bx * 2 + tmpbx) * 4 + tmpy, (by * 2 + tmpby) * 4 + tmpz, (bz * 2 + tmpbz) * 4 + tmpx, 2)] = \
			f(shared_c_v[((tmpby * 2 + (tmpz >> 1)) * 4 + (tmpbx * 2 + (tmpy >> 1))) * 4 + (tmpbz * 2 + (tmpx >> 1))], tmpy & 0b1, tmpz & 0b1, tmpx & 0b1);
	}

	template<typename T, typename F>
	void EdgeZUpSample(T *f_edge_z, T *c_edge_z, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		EdgeZUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8 + 1, f_grid.Ny / 8 + 1, f_grid.Nz / 8), dim3(8, 8, 8) >> > (f_edge_z, c_edge_z, f_grid, c_grid, f, default_value);
	}

	// 8x8x8 in fine grid view
	template<typename T, typename F>
	__global__ void NodeUpSampleKernel(T *f_node, T *c_node, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		const int bx = blockIdx.x;
		const int by = blockIdx.y;
		const int bz = blockIdx.z;

		const int idx = threadIdx.x;
		const int idy = threadIdx.y;
		const int idz = threadIdx.z;

		const int tid = (idz * 8 + idy) * 8 + idx;
		const int tmpx = tid & 0b11, tmpy = (tid & 0b1100) >> 2, tmpz = (tid & 0b110000) >> 4;
		const int tmpbx = (tid & 0b1000000) >> 6, tmpby = (tid & 0b10000000) >> 7, tmpbz = (tid & 0b100000000) >> 8;

		__shared__ T shared_c_v[64];

		if (tmpbx == 0 && tmpby == 0 && tmpbz == 0) shared_c_v[(tmpz * 4 + tmpy) * 4 + tmpx] = c_node[c_grid.node_ind(bx * 4 + tmpx, by * 4 + tmpy, bz * 4 + tmpz)];
		__syncthreads();

		f_node[f_grid.node_ind((bx * 2 + tmpbx) * 4 + tmpx, (by * 2 + tmpby) * 4 + tmpy, (bz * 2 + tmpbz) * 4 + tmpz)] = \
			f(shared_c_v[((tmpbz * 2 + (tmpz >> 1)) * 4 + (tmpby * 2 + (tmpy >> 1))) * 4 + (tmpbx * 2 + (tmpx >> 1))], tmpx & 0b1, tmpy & 0b1, tmpz & 0b1);
	}

	template<typename T, typename F>
	void NodeUpSample(T *f_node, T *c_node, grid3D f_grid, grid3D c_grid, F f, T default_value = T(0))
	{
		NodeUpSampleKernel<T, F> << <dim3(f_grid.Nx / 8+1, f_grid.Ny / 8+1, f_grid.Nz / 8+1), dim3(8, 8, 8) >> > (f_node, c_node, f_grid, c_grid, f, default_value);
	}
}
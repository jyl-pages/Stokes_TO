#include "grid3D.h"

grid3D::grid3D()
{
	Nx = 0;
	Ny = 0;
	Nz = 0;
}

grid3D::grid3D(int _Nx, int _Ny, int _Nz)
{
	Nx = _Nx;
	Ny = _Ny;
	Nz = _Nz;
}

void grid3D::init(int _Nx, int _Ny, int _Nz)
{
	Nx = _Nx;
	Ny = _Ny;
	Nz = _Nz;
}

__host__ __device__ int grid3D::cell_size() const
{
	return Nx * Ny*Nz;
}

__host__ __device__ int grid3D::face_size(int d) const
{
	return (Nx + 4 * (d == 0))*(Ny + 4 * (d == 1))*(Nz + 4 * (d == 2));
}

__host__ __device__ int grid3D::face_size_sum() const
{
	return face_size(0) + face_size(1) + face_size(2);
}

__host__ __device__ int grid3D::edge_size(int d) const
{
	return (Nx + 4 * (d != 0))*(Ny + 4 * (d != 1))*(Nz + 4 * (d != 2));
}

__host__ __device__ int grid3D::edge_size_sum() const
{
	return edge_size(0) + edge_size(1) + edge_size(2);
}

__host__ __device__ int grid3D::node_size() const
{
	return (Nx + 4)*(Ny + 4)*(Nz + 4);
}

__host__ __device__ int grid3D::cell_ind(int x, int y, int z) const
{
	int nbx = Nx >> 2, nby = Ny >> 2, nbz = Nz >> 2;
	int bx = x >> 2, by = y >> 2, bz = z >> 2;
	int idx = x & 0b11, idy = y & 0b11, idz = z & 0b11;
	return ((bz*nby + by)*nbx + bx) * 64 + ((idz * 4 + idy) * 4 + idx);
}

__host__ __device__ int grid3D::face_ind(int x, int y, int z, int d) const
{
	int nbx = (Nx >> 2) + (d == 0), nby = (Ny >> 2) + (d == 1), nbz = (Nz >> 2) + (d == 2);
	int bx = x >> 2, by = y >> 2, bz = z >> 2;
	int idx = x & 0b11, idy = y & 0b11, idz = z & 0b11;
	return ((bz*nby + by)*nbx + bx) * 64 + (d == 0)*((idx * 4 + idz) * 4 + idy) + (d == 1)*((idy * 4 + idz) * 4 + idx) + (d == 2)*((idz * 4 + idy) * 4 + idx);
}

__host__ __device__ int grid3D::edge_ind(int x, int y, int z, int d) const
{
	int nbx = (Nx >> 2) + (d != 0), nby = (Ny >> 2) + (d != 1), nbz = (Nz >> 2) + (d != 2);
	int bx = x >> 2, by = y >> 2, bz = z >> 2;
	int idx = x & 0b11, idy = y & 0b11, idz = z & 0b11;
	return ((bz*nby + by)*nbx + bx) * 64 + (d == 0)*((idz * 4 + idy) * 4 + idx) + (d == 1)*((idz * 4 + idx) * 4 + idy) + (d == 2)*((idy * 4 + idx) * 4 + idz);
}

__host__ __device__ int grid3D::node_ind(int x, int y, int z) const
{
	int nbx = (Nx >> 2) + 1, nby = (Ny >> 2) + 1, nbz = (Nz >> 2) + 1;
	int bx = x >> 2, by = y >> 2, bz = z >> 2;
	int idx = x & 0b11, idy = y & 0b11, idz = z & 0b11;
	return ((bz*nby + by)*nbx + bx) * 64 + ((idz * 4 + idy) * 4 + idx);
}

__host__ __device__ void grid3D::ind_cell(int i, int & x, int & y, int & z) const
{
	int nbx = Nx >> 2, nby = Ny >> 2, nbz = Nz >> 2;
	int idx = i & 0b11, idy = (i & 0b1100) >> 2, idz = (i & 0b110000) >> 4;
	int b = i >> 6;
	x = ((b % nbx) << 2) + idx;
	y = (((b / nbx) % nby) << 2) + idy;
	z = ((b / nbx / nby) << 2) + idz;
}

__host__ __device__ void grid3D::ind_face(int i, int d, int & x, int & y, int & z) const
{
	int nbx = (Nx >> 2) + (d == 0), nby = (Ny >> 2) + (d == 1), nbz = (Nz >> 2) + (d == 2);
	int idx = i & 0b11, idy = (i & 0b1100) >> 2, idz = (i & 0b110000) >> 4;
	int b = i >> 6;

	x = ((b % nbx) << 2) + (d == 0)*idz + (d == 1)*idx + (d == 2)*idx;
	y = (((b / nbx) % nby) << 2) + (d == 0)*idx + (d == 1)*idz + (d == 2)*idy;
	z = ((b / nbx / nby) << 2) + (d == 0)*idy + (d == 1)*idy + (d == 2)*idz;
}

__host__ __device__ void grid3D::ind_edge(int i, int d, int & x, int & y, int & z) const
{
	int nbx = (Nx >> 2) + (d != 0), nby = (Ny >> 2) + (d != 1), nbz = (Nz >> 2) + (d != 2);
	int idx = i & 0b11, idy = (i & 0b1100) >> 2, idz = (i & 0b110000) >> 4;
	int b = i >> 6;

	x = ((b % nbx) << 2) + (d == 0)*idx + (d == 1)*idy + (d == 2)*idy;
	y = (((b / nbx) % nby) << 2) + (d == 0)*idy + (d == 1)*idx + (d == 2)*idz;
	z = ((b / nbx / nby) << 2) + (d == 0)*idz + (d == 1)*idz + (d == 2)*idx;
}

__host__ __device__ void grid3D::ind_node(int i, int & x, int & y, int & z) const
{
	int nbx = (Nx >> 2) + 1, nby = (Ny >> 2) + 1, nbz = (Nz >> 2) + 1;
	int idx = i & 0b11, idy = (i & 0b1100) >> 2, idz = (i & 0b110000) >> 4;
	int b = i >> 6;
	x = ((b % nbx) << 2) + idx;
	y = (((b / nbx) % nby) << 2) + idy;
	z = ((b / nbx / nby) << 2) + idz;
}

__host__ __device__ bool grid3D::cell_valid(int x, int y, int z) const
{
	return x >= 0 && x < Nx && y >= 0 && y < Ny && z >= 0 && z < Nz;
}

__host__ __device__ bool grid3D::face_valid(int x, int y, int z, int d) const
{
	return x >= 0 && x < (Nx + (d == 0)) && y >= 0 && y < (Ny + (d == 1)) && z >= 0 && z < (Nz + (d == 2));
}

__host__ __device__ bool grid3D::edge_valid(int x, int y, int z, int d) const
{
	return x >= 0 && x < (Nx + (d != 0)) && y >= 0 && y < (Ny + (d != 1)) && z >= 0 && z < (Nz + (d != 2));
}

__host__ __device__ bool grid3D::node_valid(int x, int y, int z) const
{
	return x >= 0 && x <= Nx && y >= 0 && y <= Ny && z >= 0 && z <= Nz;
}

void grid3DOperator::ApplyMask(const grid3D & grid, Scalar * v, bool * fixed, int N)
{
}

void grid3DOperator::ApplyFaceRedudantFixed(const grid3D & grid, bool * fixed, int d)
{
}

void grid3DOperator::ApplyEdgeRedudantFixed(const grid3D & grid, bool * fixed, int d)
{
}

void grid3DOperator::ApplyNodeRedudantFixed(const grid3D & grid, bool * fixed)
{
}

// for blockDim = (4, 4, 4)
// iterate through node
__global__ void D0MappingKernel(grid3D grid, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z, Scalar * node)
{
	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	const int x = bx * 4 + idx;
	const int y = by * 4 + idy;
	const int z = bz * 4 + idz;

	__shared__ Scalar shared_node_data[64];

	Scalar grad_x(0), grad_y(0), grad_z(0);

	// 4x4x4 node data load
	shared_node_data[(idz * 4 + idy) * 4 + idx] = grid.node_valid(x, y, z) ? node[grid.node_ind(x, y, z)] : (Scalar)0;
	__syncthreads();

	grad_x -= shared_node_data[(idz * 4 + idy) * 4 + idx];
	if (idx != 3) grad_x += shared_node_data[(idz * 4 + idy) * 4 + (idx + 1)];

	grad_y -= shared_node_data[(idz * 4 + idy) * 4 + idx];
	if (idy != 3) grad_y += shared_node_data[(idz * 4 + (idy + 1)) * 4 + idx];

	grad_z -= shared_node_data[(idz * 4 + idy) * 4 + idx];
	if (idz != 3) grad_z += shared_node_data[((idz + 1) * 4 + idy) * 4 + idx];

	if (idx == 3)
		if (grid.node_valid(x + 1, y, z)) grad_x += node[grid.node_ind(x + 1, y, z)];
	if (idy == 3)
		if (grid.node_valid(x, y + 1, z)) grad_y += node[grid.node_ind(x, y + 1, z)];
	if (idz == 3)
		if (grid.node_valid(x, y, z + 1)) grad_z += node[grid.node_ind(x, y, z + 1)];

	if (grid.edge_valid(x, y, z, 0)) edge_x[grid.edge_ind(x, y, z, 0)] = grad_x;
	if (grid.edge_valid(x, y, z, 1)) edge_y[grid.edge_ind(x, y, z, 1)] = grad_y;
	if (grid.edge_valid(x, y, z, 2)) edge_z[grid.edge_ind(x, y, z, 2)] = grad_z;
}

void grid3DOperator::D0Mapping(const grid3D & grid, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z, Scalar * node)
{
	cudaMemset(edge_x, 0, sizeof(Scalar)*grid.edge_size(0));
	cudaMemset(edge_y, 0, sizeof(Scalar)*grid.edge_size(1));
	cudaMemset(edge_z, 0, sizeof(Scalar)*grid.edge_size(2));
	D0MappingKernel << <dim3((grid.Nx >> 2) + 1, (grid.Ny >> 2) + 1, (grid.Nz >> 2) + 1), dim3(4, 4, 4) >> > (grid, edge_x, edge_y, edge_z, node);
}

// for blockDim = (4, 4, 4)
// iterate through node
__global__ void D1MappingKernel(grid3D grid, Scalar * face_x, Scalar * face_y, Scalar * face_z, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z)
{
	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	const int x = bx * 4 + idx;
	const int y = by * 4 + idy;
	const int z = bz * 4 + idz;

	__shared__ Scalar shared_edge_data[64];

	Scalar curl_x(0), curl_y(0), curl_z(0);

	{
		// 4x4x4 edge_x data load
		const int ex = bx * 4 + idx, ey = by * 4 + idy, ez = bz * 4 + idz;
		shared_edge_data[(idz * 4 + idy) * 4 + idx] = grid.edge_valid(ex, ey, ez, 0) ? edge_x[grid.edge_ind(ex, ey, ez, 0)] : (Scalar)0;
		__syncthreads();

		// for curl_y
		// front edge_x
		curl_y += shared_edge_data[(idz * 4 + idy) * 4 + idx];
		// back edge_x
		curl_y += idz == 3 ? (Scalar)0 : -shared_edge_data[((idz + 1) * 4 + idy) * 4 + idx];

		// for curl_z
		// down edge_x
		curl_z += -shared_edge_data[(idz * 4 + idy) * 4 + idx];
		// up edge_x
		curl_z += idy == 3 ? (Scalar)0 : shared_edge_data[(idz * 4 + (idy + 1)) * 4 + idx];
	}
	__syncthreads();

	{
		// 4x4x4 edge_y data load
		const int ex = bx * 4 + idy, ey = by * 4 + idx, ez = bz * 4 + idz;
		shared_edge_data[(idz * 4 + idy) * 4 + idx] = grid.edge_valid(ex, ey, ez, 1) ? edge_y[grid.edge_ind(ex, ey, ez, 1)] : (Scalar)0;
		__syncthreads();

		// for curl_x
		// front edge_y
		curl_x += -shared_edge_data[(idz * 4 + idx) * 4 + idy];
		// back edge_y
		curl_x += idz == 3 ? (Scalar)0 : shared_edge_data[((idz + 1) * 4 + idx) * 4 + idy];

		// for curl_z
		// left edge_y
		curl_z +=  shared_edge_data[(idz * 4 + idx) * 4 + idy];
		// right edge_y
		curl_z += idx == 3 ? (Scalar)0 : -shared_edge_data[(idz * 4 + (idx + 1)) * 4 + idy];
	}
	__syncthreads();

	{
		// 4x4x4 edge_z data load
		const int ex = bx * 4 + idy, ey = by * 4 + idz, ez = bz * 4 + idx;
		shared_edge_data[(idz * 4 + idy) * 4 + idx] = grid.edge_valid(ex, ey, ez, 2) ? edge_z[grid.edge_ind(ex, ey, ez, 2)] : (Scalar)0;
		__syncthreads();

		// for curl_x
		// down edge_z
		curl_x += shared_edge_data[(idy * 4 + idx) * 4 + idz];
		// up edge_z
		curl_x += idy == 3 ? (Scalar)0 : -shared_edge_data[((idy + 1) * 4 + idx) * 4 + idz];

		// for curl_y
		// left edge_z
		curl_y += -shared_edge_data[(idy * 4 + idx) * 4 + idz];
		// right edge_z
		curl_y += idx == 3 ? (Scalar)0 : shared_edge_data[(idy * 4 + (idx + 1)) * 4 + idz];
	}
	__syncthreads();

	if (idx == 3)
	{
		if (grid.edge_valid(x + 1, y, z, 2)) curl_y += edge_z[grid.edge_ind(x + 1, y, z, 2)];
		if (grid.edge_valid(x + 1, y, z, 1)) curl_z += -edge_y[grid.edge_ind(x + 1, y, z, 1)];
	}

	if (idy == 3)
	{
		if (grid.edge_valid(x, y + 1, z, 2)) curl_x += -edge_z[grid.edge_ind(x, y + 1, z, 2)];
		if (grid.edge_valid(x, y + 1, z, 0)) curl_z += +edge_x[grid.edge_ind(x, y + 1, z, 0)];
	}

	if (idz == 3)
	{
		if (grid.edge_valid(x, y, z + 1, 1)) curl_x += edge_y[grid.edge_ind(x, y, z + 1, 1)];
		if (grid.edge_valid(x, y, z + 1, 0)) curl_y += -edge_x[grid.edge_ind(x, y, z + 1, 0)];
	}

	if (grid.face_valid(x, y, z, 0)) face_x[grid.face_ind(x, y, z, 0)] = curl_x;
	if (grid.face_valid(x, y, z, 1)) face_y[grid.face_ind(x, y, z, 1)] = curl_y;
	if (grid.face_valid(x, y, z, 2)) face_z[grid.face_ind(x, y, z, 2)] = curl_z;
}

void grid3DOperator::D1Mapping(const grid3D & grid, Scalar * face_x, Scalar * face_y, Scalar * face_z, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z)
{
	cudaMemset(face_x, 0, sizeof(Scalar)*grid.face_size(0));
	cudaMemset(face_y, 0, sizeof(Scalar)*grid.face_size(1));
	cudaMemset(face_z, 0, sizeof(Scalar)*grid.face_size(2));
	D1MappingKernel << <dim3((grid.Nx >> 2) + 1, (grid.Ny >> 2) + 1, (grid.Nz >> 2) + 1), dim3(4, 4, 4) >> > (grid, face_x, face_y, face_z, edge_x, edge_y, edge_z);
}

// for blockDim = (4, 4, 4)
// iterate through cell
__global__ void D2MappingKernel(grid3D grid, Scalar *cell, Scalar * face_x, Scalar * face_y, Scalar * face_z)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;
	const int nbz = gridDim.z;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	__shared__ Scalar shared_face_data[80];

	Scalar div;
	{
		// 5x4x4 face_x data load
		shared_face_data[(idz * 4 + idy) * 4 + idx] = face_x[grid.face_ind(bx * 4 + idz, by * 4 + idx, bz * 4 + idy, 0)];
		if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_x[grid.face_ind(bx * 4 + 4, by * 4 + idx, bz * 4 + idy, 0)];
		__syncthreads();

		// left x-axis faces
		div = -shared_face_data[(idx * 4 + idz) * 4 + idy];

		// right x-axis faces
		div += shared_face_data[((idx + 1) * 4 + idz) * 4 + idy];
	}
	__syncthreads();

	{
		// 4x5x4 face_y data load
		shared_face_data[(idz * 4 + idy) * 4 + idx] = face_y[grid.face_ind(bx * 4 + idx, by * 4 + idz, bz * 4 + idy, 1)];
		if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_y[grid.face_ind(bx * 4 + idx, by * 4 + 4, bz * 4 + idy, 1)];
		__syncthreads();

		// down y-axis faces
		div += -shared_face_data[(idy * 4 + idz) * 4 + idx];

		// up y-axis faces
		div += shared_face_data[((idy + 1) * 4 + idz) * 4 + idx];
	}
	__syncthreads();

	{
		// 4x4x5 face_y data load
		shared_face_data[(idz * 4 + idy) * 4 + idx] = face_z[grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 2)];
		if (idz == 0) shared_face_data[(4 * 4 + idy) * 4 + idx] = face_z[grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + 4, 2)];
		__syncthreads();

		// front z-axis faces
		div += -shared_face_data[(idz * 4 + idy) * 4 + idx];

		// back y-axis faces
		div += shared_face_data[((idz + 1) * 4 + idy) * 4 + idx];
	}
	__syncthreads();

	cell[grid.cell_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz)] = div;
}


void grid3DOperator::D2Mapping(const grid3D & grid, Scalar * cell, Scalar * face_x, Scalar * face_y, Scalar * face_z)
{
	cudaMemset(cell, 0, sizeof(Scalar)*grid.cell_size());
	D2MappingKernel << <dim3((grid.Nx >> 2), (grid.Ny >> 2), (grid.Nz >> 2)), dim3(4, 4, 4) >> > (grid, cell, face_x, face_y, face_z);
}

// for blockDim = (8, 8)
// iterate through cell
__global__ void Cod0MappingKernel(grid3D grid, Scalar *face_x, Scalar *face_y, Scalar *face_z, Scalar *cell)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;
	const int nbz = gridDim.z;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	__shared__ Scalar shared_cell_data[64];

	// 8x8 cell data load
	shared_cell_data[(idz * 4 + idy) * 4 + idx] = cell[grid.cell_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz)];
	__syncthreads();

	const Scalar cell_data = shared_cell_data[(idz * 4 + idy) * 4 + idx];
	int face_ind;

	// x-axis faces
	face_ind = grid.face_ind(bx * 4 + (idx + 1), by * 4 + idy, bz * 4 + idz, 0);
	if (idx != 3) face_x[face_ind] = shared_cell_data[(idz * 4 + idy) * 4 + (idx + 1)] - cell_data;

	// y-axis faces
	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + (idy + 1), bz * 4 + idz, 1);
	if (idy != 3) face_y[face_ind] = shared_cell_data[(idz * 4 + (idy + 1)) * 4 + idx] - cell_data;

	// z-axis faces
	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + (idz + 1), 2);
	if (idz != 3) face_z[face_ind] = shared_cell_data[((idz + 1) * 4 + idy) * 4 + idx] - cell_data;

	// x-axis shared faces
	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 0);
	if (idx == 0) atomicAdd(face_x + face_ind, cell_data);

	face_ind = grid.face_ind(bx * 4 + (idx + 1), by * 4 + idy, bz * 4 + idz, 0);
	if (idx == 3) atomicAdd(face_x + face_ind, -cell_data);

	// y-axis shared faces
	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 1);
	if (idy == 0) atomicAdd(face_y + face_ind, cell_data);

	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + (idy + 1), bz * 4 + idz, 1);
	if (idy == 3) atomicAdd(face_y + face_ind, -cell_data);

	// x-axis shared faces
	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + idz, 2);
	if (idz == 0) atomicAdd(face_z + face_ind, cell_data);

	face_ind = grid.face_ind(bx * 4 + idx, by * 4 + idy, bz * 4 + (idz + 1), 2);
	if (idz == 3) atomicAdd(face_z + face_ind, -cell_data);
}


void grid3DOperator::Cod0Mapping(const grid3D & grid, Scalar * face_x, Scalar * face_y, Scalar * face_z, Scalar * cell)
{
	cudaMemset(face_x, 0, sizeof(Scalar)*grid.face_size(0));
	cudaMemset(face_y, 0, sizeof(Scalar)*grid.face_size(1));
	cudaMemset(face_z, 0, sizeof(Scalar)*grid.face_size(2));
	Cod0MappingKernel << <dim3((grid.Nx >> 2), (grid.Ny >> 2), (grid.Nz >> 2)), dim3(4, 4, 4) >> > (grid, face_x, face_y, face_z, cell);
}

// for blockDim = (4, 4, 4)
// iterate through nodes
__global__ void Cod1MappingKernel(grid3D grid, Scalar *edge_x, Scalar* edge_y, Scalar *edge_z, Scalar *face_x, Scalar *face_y, Scalar *face_z)
{
	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	const int x = bx * 4 + idx;
	const int y = by * 4 + idy;
	const int z = bz * 4 + idz;

	__shared__ Scalar shared_face_data[64];

	Scalar curl_x(0), curl_y(0), curl_z(0);

	{
		// 4x4x4 face_x data load
		const int fx = bx * 4 + idz, fy = by * 4 + idx, fz = bz * 4 + idy;
		shared_face_data[(idz * 4 + idy) * 4 + idx] = grid.face_valid(fx, fy, fz, 0) ? face_x[grid.face_ind(fx, fy, fz, 0)] : (Scalar)0;
		__syncthreads();

		// for curl_y
		// front face_x
		curl_y += idz == 0 ? (Scalar)0 : shared_face_data[(idx * 4 + idz - 1) * 4 + idy];
		// back faces_x
		curl_y += -shared_face_data[(idx * 4 + idz) * 4 + idy];

		// for curl_z
		// down face_x
		curl_z += idy == 0 ? (Scalar)0 : -shared_face_data[(idx * 4 + idz) * 4 + idy - 1];
		// up face_x
		curl_z += shared_face_data[(idx * 4 + idz) * 4 + idy];
	}
	__syncthreads();

	{
		// 4x4x4 face_y data load
		const int fx = bx * 4 + idx, fy = by * 4 + idz, fz = bz * 4 + idy;
		shared_face_data[(idz * 4 + idy) * 4 + idx] = grid.face_valid(fx, fy, fz, 1) ? face_y[grid.face_ind(fx, fy, fz, 1)] : (Scalar)0;
		__syncthreads();

		// for curl_x
		// front face_y
		curl_x += idz == 0 ? (Scalar)0 : -shared_face_data[(idy * 4 + idz - 1) * 4 + idx];
		// back faces_y
		curl_x += shared_face_data[(idy * 4 + idz) * 4 + idx];

		// for curl_z
		// left face_y
		curl_z += idx == 0 ? (Scalar)0 : shared_face_data[(idy * 4 + idz) * 4 + idx - 1];
		// right face_y
		curl_z += -shared_face_data[(idy * 4 + idz) * 4 + idx];
	}
	__syncthreads();

	{
		// 4x4x4 face_z data load
		const int fx = bx * 4 + idx, fy = by * 4 + idy, fz = bz * 4 + idz;
		shared_face_data[(idz * 4 + idy) * 4 + idx] = grid.face_valid(fx, fy, fz, 2) ? face_z[grid.face_ind(fx, fy, fz, 2)] : (Scalar)0;
		__syncthreads();

		// for curl_x
		// down face_z
		curl_x += idy == 0 ? (Scalar)0 : shared_face_data[(idz * 4 + idy - 1) * 4 + idx];
		// up faces_z
		curl_x += -shared_face_data[(idz * 4 + idy) * 4 + idx];

		// for curl_y
		// left face_z
		curl_y += idx == 0 ? (Scalar)0 : -shared_face_data[(idz * 4 + idy) * 4 + idx - 1];
		// right face_z
		curl_y += shared_face_data[(idz * 4 + idy) * 4 + idx];
	}
	__syncthreads();

	if (idx == 0)
	{
		if (grid.face_valid(x - 1, y, z, 2)) curl_y += -face_z[grid.face_ind(x - 1, y, z, 2)];
		if (grid.face_valid(x - 1, y, z, 1)) curl_z += face_y[grid.face_ind(x - 1, y, z, 1)];
	}

	if (idy == 0)
	{
		if (grid.face_valid(x, y - 1, z, 2)) curl_x += face_z[grid.face_ind(x, y - 1, z, 2)];
		if (grid.face_valid(x, y - 1, z, 0)) curl_z += -face_x[grid.face_ind(x, y - 1, z, 0)];
	}

	if (idz == 0)
	{
		if (grid.face_valid(x, y, z - 1, 1)) curl_x += -face_y[grid.face_ind(x, y, z - 1, 1)];
		if (grid.face_valid(x, y, z - 1, 0)) curl_y += face_x[grid.face_ind(x, y, z - 1, 0)];
	}

	if (grid.edge_valid(x, y, z, 0)) edge_x[grid.edge_ind(x, y, z, 0)] = curl_x;
	if (grid.edge_valid(x, y, z, 1)) edge_y[grid.edge_ind(x, y, z, 1)] = curl_y;
	if (grid.edge_valid(x, y, z, 2)) edge_z[grid.edge_ind(x, y, z, 2)] = curl_z;
}


void grid3DOperator::Cod1Mapping(const grid3D & grid, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z, Scalar * face_x, Scalar * face_y, Scalar * face_z)
{
	cudaMemset(edge_x, 0, sizeof(Scalar)*grid.edge_size(0));
	cudaMemset(edge_y, 0, sizeof(Scalar)*grid.edge_size(1));
	cudaMemset(edge_z, 0, sizeof(Scalar)*grid.edge_size(2));

	Cod1MappingKernel << <dim3((grid.Nx >> 2) + 1, (grid.Ny >> 2) + 1, (grid.Nz >> 2) + 1), dim3(4, 4, 4) >> > (grid, edge_x, edge_y, edge_z, face_x, face_y, face_z);
}

__global__ void Cod2MappingKernel(grid3D grid, Scalar * node, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z)
{
	const int bx = blockIdx.x;
	const int by = blockIdx.y;
	const int bz = blockIdx.z;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;
	const int idz = threadIdx.z;

	const int x = bx * 4 + idx;
	const int y = by * 4 + idy;
	const int z = bz * 4 + idz;

	__shared__ Scalar shared_edge_data[64];

	Scalar div(0);

	{
		// 4x4x4 edge_x data load
		int ex = bx * 4 + idx, ey = by * 4 + idy, ez = bz * 4 + idz;
		shared_edge_data[(idz * 4 + idy) * 4 + idx] = grid.edge_valid(ex, ey, ez, 0) ? edge_x[grid.edge_ind(ex, ey, ez, 0)] : (Scalar)0;
		__syncthreads();
		
		div += shared_edge_data[(idz * 4 + idy) * 4 + idx];
		if (idx != 0) div -= shared_edge_data[(idz * 4 + idy) * 4 + (idx - 1)];

		if (idx == 0)
			if (grid.edge_valid(x - 1, y, z, 0)) 
				div -= edge_x[grid.edge_ind(x - 1, y, z, 0)];
	}
	__syncthreads();

	{
		// 4x4x4 edge_y data load
		int ex = bx * 4 + idy, ey = by * 4 + idx, ez = bz * 4 + idz;
		shared_edge_data[(idz * 4 + idy) * 4 + idx] = grid.edge_valid(ex, ey, ez, 1) ? edge_y[grid.edge_ind(ex, ey, ez, 1)] : (Scalar)0;
		__syncthreads();

		div += shared_edge_data[(idz * 4 + idx) * 4 + idy];
		if (idy != 0) div -= shared_edge_data[(idz * 4 + idx) * 4 + (idy - 1)];

		if (idy == 0)
			if (grid.edge_valid(x, y - 1, z, 1))
				div -= edge_y[grid.edge_ind(x, y - 1, z, 1)];
	}
	__syncthreads();

	{
		// 4x4x4 edge_z data load
		int ex = bx * 4 + idy, ey = by * 4 + idz, ez = bz * 4 + idx;
		shared_edge_data[(idz * 4 + idy) * 4 + idx] = grid.edge_valid(ex, ey, ez, 2) ? edge_z[grid.edge_ind(ex, ey, ez, 2)] : (Scalar)0;
		__syncthreads();

		div += shared_edge_data[(idy * 4 + idx) * 4 + idz];
		if (idz != 0) div -= shared_edge_data[(idy * 4 + idx) * 4 + (idz - 1)];

		if (idz == 0)
			if (grid.edge_valid(x, y, z - 1, 2))
				div -= edge_z[grid.edge_ind(x, y, z - 1, 2)];
	}
	__syncthreads();

	if (grid.node_valid(x, y, z)) node[grid.node_ind(x, y, z)] = div;
}

void grid3DOperator::Cod2Mapping(const grid3D & grid, Scalar * node, Scalar * edge_x, Scalar * edge_y, Scalar * edge_z)
{
	cudaMemset(node, 0, sizeof(Scalar)*grid.node_size());
	Cod2MappingKernel << <dim3((grid.Nx >> 2) + 1, (grid.Ny >> 2) + 1, (grid.Nz >> 2) + 1), dim3(4, 4, 4) >> > (grid, node, edge_x, edge_y, edge_z);
}

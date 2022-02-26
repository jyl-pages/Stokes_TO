#include "form01mapping.h"
#include "gpuUtils.h"

void grid2D::init(int _Nx,int _Ny)
{
	Nx = _Nx;
	Ny = _Ny;
}

__host__ __device__ int grid2D::cell_size() const
{
	return Nx * Ny;
}

__host__ __device__ int grid2D::face_size(int d) const
{
	return (Nx + 8 * (d == 0))*(Ny + 8 * (d == 1));
}

__host__ __device__ int grid2D::node_size() const
{
	return (Nx + 8)*(Ny + 8);
}

__host__ __device__ int grid2D::cell_ind(int x, int y) const
{
	return ((y >> 3)*(Nx >> 3) + (x >> 3)) * 64 + ((y & 7) * 8 + (x & 7));
}

__host__ __device__ int grid2D::face_ind(int x, int y, int d) const
{
	int b_ind = (y >> 3)*((Nx>>3) + (d == 0)) + (x >> 3);
	int idx = x & 7, idy = y & 7;
	return b_ind * 64 + ((d == 0)*(idx * 8 + idy) + (d == 1)*(idy * 8 + idx));
}

__host__ __device__ int grid2D::node_ind(int x, int y) const
{
	return ((y >> 3)*((Nx >> 3) + 1) + (x >> 3)) * 64 + ((y & 7) * 8 + (x & 7));
}

__host__ __device__ void grid2D::ind_cell(int i, int & x, int & y) const
{
	int idx = i & 0b111;
	int idy = (i & 0b111000) >> 3;
	int b = i >> 6;
	x = ((b % (Nx >> 3))<<3) + idx;
	y = ((b / (Nx >> 3))<<3) + idy;
}

__host__ __device__ void grid2D::ind_face(int i, int d, int & x, int & y) const
{
}

__host__ __device__ void grid2D::ind_node(int i, int & x, int & y) const
{
}

template<int k>
__global__ void apply_mask(Scalar *v, bool *mask, int N)
{
	const int ind = blockIdx.x*blockDim.x + threadIdx.x;
	if (ind >= N) return;
	const bool b = mask[ind];

	for (int i = 0; i < k; i++)
		if (b) v[ind + i * N] = Scalar(0);
}

void grid2DOperator::ApplyMask(const grid2D & grid, Scalar * v, bool * fixed, int N)
{
	cwise_mapping_wrapper(v, fixed, [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; }, N);
}

// for blockDim = (8, 8)
__global__ void ApplyXFaceRedudantFixedKernel(grid2D grid, bool *fixed, int xMax, int yMax)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	const int x = bx * 8 + idy;
	const int y = by * 8 + idx;

	const int face_ind = grid.face_ind(x, y, 0);

	fixed[face_ind] |= (x >= xMax || y >= yMax);
	//if (x >= xMax || y >= yMax) fixed[face_ind] = true;
}

// for blockDim = (8, 8)
__global__ void ApplyYFaceRedudantFixedKernel(grid2D grid, bool *fixed, int xMax, int yMax)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	const int x = bx * 8 + idx;
	const int y = by * 8 + idy;

	const int face_ind = grid.face_ind(x, y, 1);

	fixed[face_ind] |= (x >= xMax || y >= yMax);
	//if (x >= xMax || y >= yMax) fixed[face_ind] = true;
}

void grid2DOperator::ApplyFaceRedudantFixed(const grid2D& grid, bool * fixed, int d)
{
	if (d == 0)
		ApplyXFaceRedudantFixedKernel << <dim3((grid.Nx>>3) + 1, (grid.Ny>>3)), dim3(8, 8) >> > (grid, fixed, grid.Nx + 1, grid.Ny);
	else
		ApplyYFaceRedudantFixedKernel << <dim3((grid.Nx>>3), (grid.Ny>>3)+1), dim3(8, 8) >> > (grid, fixed, grid.Nx, grid.Ny + 1);
}

__global__ void ApplyNodeRedudantFixedKernel(grid2D grid, bool *fixed, int xMax, int yMax)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	const int x = bx * 8 + idx;
	const int y = by * 8 + idy;

	const int node_ind = grid.node_ind(x, y);

	fixed[node_ind] |= (x >= xMax || y >= yMax);
}

void grid2DOperator::ApplyNodeRedudantFixed(const grid2D& grid, bool *fixed)
{
	ApplyNodeRedudantFixedKernel << <dim3((grid.Nx >> 3) + 1, (grid.Ny >> 3) + 1), dim3(8, 8) >> > (grid, fixed, grid.Nx + 1, grid.Ny + 1);
}

// for blockDim = (8, 8)
// iterate through cell
__global__ void Cod0MappingKernel(grid2D grid, Scalar *face_x, Scalar *face_y, Scalar *cell)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ Scalar shared_cell_data[64];

	// 8x8 cell data load
	shared_cell_data[idy * 8 + idx] = cell[grid.cell_ind(bx * 8 + idx, by * 8 + idy)];
	__syncthreads();

	const Scalar cell_data = shared_cell_data[idy * 8 + idx];
	int face_ind;

	// x-axis faces
	face_ind = grid.face_ind(bx * 8 + idx + 1, by * 8 + idy, 0);
	if (idx != 7) face_x[face_ind] = shared_cell_data[idy * 8 + (idx + 1)] - cell_data;

	// y-axis faces
	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + idy + 1, 1);
	if (idy != 7) face_y[face_ind] = shared_cell_data[(idy + 1) * 8 + idx] - cell_data;

	// x-axis shared faces
	face_ind = grid.face_ind(bx * 8 + 0, by * 8 + idy, 0);
	if (idx == 0) atomicAdd(face_x + face_ind, cell_data);

	face_ind = grid.face_ind((bx + 1) * 8 + 0, by * 8 + idy, 0);
	if (idx == 7) atomicAdd(face_x + face_ind, -cell_data);

	// y-axis shared faces
	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + 0, 1);
	if (idy == 0) atomicAdd(face_y + face_ind, cell_data);

	face_ind = grid.face_ind(bx * 8 + idx, (by + 1) * 8 + 0, 1);
	if (idy == 7) atomicAdd(face_y + face_ind, -cell_data);
}

void grid2DOperator::Cod0Mapping(const grid2D& grid, Scalar *face_x, Scalar *face_y, Scalar *cell)
{
	cudaMemset(face_x, 0, sizeof(Scalar)*grid.face_size(0));
	cudaMemset(face_y, 0, sizeof(Scalar)*grid.face_size(1));
	Cod0MappingKernel << <dim3((grid.Nx >> 3), (grid.Ny >> 3)), dim3(8, 8) >> > (grid, face_x, face_y, cell);
}

// for blockDim = (8, 8)
// iterate through cell
__global__ void D1MappingKernel(grid2D grid, Scalar *cell, Scalar *face_x, Scalar *face_y)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ Scalar shared_face_data[72];
	
	Scalar div;
	{
		// 9x8 face_x data load
		shared_face_data[idy * 8 + idx] = face_x[grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)];
		if (idy == 0) shared_face_data[64 + idx] = face_x[grid.face_ind((bx + 1) * 8 + idy, by * 8 + idx, 0)];
		__syncthreads();

		// left x-axis faces
		div = -shared_face_data[idx * 8 + idy];

		// right x-axis faces
		div += shared_face_data[(idx + 1) * 8 + idy];
	}
	__syncthreads();

	{
		// 8x9 face_y data load
		shared_face_data[idy * 8 + idx] = face_y[grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)];
		if (idy == 0) shared_face_data[64 + idx] = face_y[grid.face_ind(bx * 8 + idx, (by + 1) * 8 + idy, 1)];
		__syncthreads();

		// down y-axis faces
		div += shared_face_data[idy * 8 + idx];

		// up y-axis faces
		div += -shared_face_data[(idy + 1) * 8 + idx];
	}

	cell[grid.cell_ind(bx * 8 + idx, by * 8 + idy)] = div;
}


void grid2DOperator::D1Mapping(const grid2D& grid, Scalar *cell, Scalar *face_x, Scalar *face_y)
{
	cudaMemset(cell, 0, sizeof(Scalar)*grid.cell_size());
	D1MappingKernel << <dim3((grid.Nx>>3), (grid.Ny>>3)), dim3(8, 8) >> > (grid, cell, face_x, face_y);
}

// for blockDim = (8, 8)
// iterate through nodes
__global__ void Cod1MappingKernel(grid2D grid, Scalar *node, Scalar *face_x, Scalar *face_y)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ Scalar shared_face_data[72];

	Scalar curl;

	{
		// 9x8 face_x data load
		shared_face_data[idy * 9 + (idx + 1)] = (by == nby - 1) ? (Scalar)0 : face_x[grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)];
		if (idx == 7) shared_face_data[idy * 9 + 0] = (by == 0 || (idy != 0 && bx == nbx - 1)) ? (Scalar)0 : face_x[grid.face_ind(bx * 8 + idy, (by - 1) * 8 + idx, 0)];
		__syncthreads();

		// down x-axis faces
		curl = shared_face_data[idx * 9 + idy];

		// up x-axis faces
		curl += -shared_face_data[idx * 9 + (idy + 1)];

	}
	__syncthreads();

	{
		// 9x8 face_y data load
		shared_face_data[idy * 9 + (idx + 1)] = (bx == nbx - 1) ? (Scalar)0 : face_y[grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)];
		if (idx == 7) shared_face_data[idy * 9 + 0] = (bx == 0 || (idy != 0 && by == nby - 1)) ? (Scalar)0 : face_y[grid.face_ind((bx - 1) * 8 + idx, by * 8 + idy, 1)];
		__syncthreads();

		// left y-axis faces
		curl += -shared_face_data[idy * 9 + idx];

		// right y-axis faces
		curl += shared_face_data[idy * 9 + (idx + 1)];
	}

	node[grid.node_ind(bx * 8 + idx, by * 8 + idy)] = curl;
}

//// for blockDim = (8, 8)
//// iterate through cells
//__global__ void Cod1MappingKernel(grid2D grid, Scalar *node, Scalar *face_x, Scalar *face_y)
//{
//	const int nbx = gridDim.x;
//	const int nby = gridDim.y;
//
//	const int bx = blockIdx.x;
//	const int by = blockIdx.y;
//
//	const int idx = threadIdx.x;
//	const int idy = threadIdx.y;
//
//	__shared__ Scalar shared_face_data[72];
//	__shared__ Scalar shared_node_data[81];
//
//	{
//		shared_node_data[idy * 9 + idx] = (Scalar)0;
//		if (idx == 0) shared_node_data[idy * 9 + 8] = (Scalar)0;
//		if (idy == 0) shared_node_data[8 * 9 + idx] = (Scalar)0;
//		if (idx == 0 && idy == 0) shared_node_data[8 * 9 + 8] = (Scalar)0;
//		__syncthreads();
//	}
//
//	{
//		// 9x8 face_x data load
//		shared_face_data[idy * 8 + idx] = face_x[grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)];
//		if (idy == 0) shared_face_data[64 + idx] = face_x[grid.face_ind((bx + 1) * 8 + idy, by * 8 + idx, 0)];
//		__syncthreads();
//
//		// down x-axis faces
//		shared_node_data[(idy + 1) * 9 + idx] += shared_face_data[idx * 8 + idy];
//		if (idx == 0) shared_node_data[(idy + 1) * 9 + 8] += shared_face_data[8 * 8 + idy];
//		__syncthreads();
//
//		// up x-axis faces
//		shared_node_data[idy * 9 + idx] -= shared_face_data[idx * 8 + idy];
//		if (idx == 0) shared_node_data[idy * 9 + 8] -= shared_face_data[8 * 8 + idy];
//		__syncthreads();
//	}
//
//	{
//		// 9x8 face_y data load
//		shared_face_data[idy * 8 + idx] = face_y[grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)];
//		if (idy == 0) shared_face_data[64 + idx] = face_y[grid.face_ind(bx * 8 + idx, (by + 1) * 8 + idy, 1)];
//		__syncthreads();
//
//		// left y-axis faces
//		shared_node_data[idy * 9 + (idx + 1)] -= shared_face_data[idy * 8 + idx];
//		if(idy==0) shared_node_data[8 * 9 + (idx + 1)] -= shared_face_data[8 * 8 + idx];
//		__syncthreads();
//
//		// right y-axis faces
//		shared_node_data[idy * 9 + idx] += shared_face_data[idy * 8 + idx];
//		if (idy == 0) shared_node_data[8 * 9 + idx] += shared_face_data[8 * 8 + idx];
//		__syncthreads();
//	}
//
//	{
//		atomicAdd(node + grid.node_ind(bx * 8 + idx, by * 8 + idy), shared_node_data[idy * 9 + idx]);
//		if (idx == 0) atomicAdd(node + grid.node_ind(bx * 8 + 8, by * 8 + idy), shared_node_data[idy * 9 + 8]);
//		if (idy == 0) atomicAdd(node + grid.node_ind(bx * 8 + idx, by * 8 + 8), shared_node_data[8 * 9 + idx]);
//		if (idx == 0 && idy == 0) atomicAdd(node + grid.node_ind(bx * 8 + 8, by * 8 + 8), shared_node_data[8 * 9 + 8]);
//	}
//}

void grid2DOperator::Cod1Mapping(const grid2D& grid, Scalar *node, Scalar *face_x, Scalar *face_y)
{
	cudaMemset(node, 0, sizeof(Scalar)*grid.node_size());
	Cod1MappingKernel << <dim3((grid.Nx >> 3) + 1, (grid.Ny >> 3) + 1), dim3(8, 8) >> > (grid, node, face_x, face_y);
}

// for blockDim = (8, 8)
// iterate through nodes
__global__ void D0MappingKernel(grid2D grid, Scalar *face_x, Scalar *face_y, Scalar *node)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ Scalar shared_node_data[64];

	// 8x8 node data load
	shared_node_data[idy * 8 + idx] = (idx == 0 || bx != nbx - 1) && (idy == 0 || by != nby - 1) ? node[grid.node_ind(bx * 8 + idx, by * 8 + idy)] : (Scalar)0;
	__syncthreads();

	const Scalar node_data = shared_node_data[idy * 8 + idx];
	int face_ind;

	// x-axis face
	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + idy, 0);
	if (idy != 7 && by != nby - 1 && (idx == 0 || bx != nbx - 1)) face_x[face_ind] = shared_node_data[(idy + 1) * 8 + idx] - node_data;

	// y-axis face
	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + idy, 1);
	if (idx != 7 && bx != nbx - 1 && (idy == 0 || by != nby - 1)) face_y[face_ind] = shared_node_data[idy * 8 + (idx + 1)] - node_data;

	// x-axis shared face
	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + idy - 1, 0);
	if (idy == 0 && by != 0 && (idx == 0 || bx != nbx - 1)) atomicAdd(face_x + face_ind, node_data);

	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + idy, 0);
	if (idy == 7 && by != nby - 1 && (idx == 0 || bx != nbx - 1)) atomicAdd(face_x + face_ind, -node_data);

	// y-axis shared face
	face_ind = grid.face_ind(bx * 8 + idx - 1, by * 8 + idy, 1);
	if (idx == 0 && bx != 0 && (idy == 0 || by != nby - 1)) atomicAdd(face_y + face_ind, node_data);

	face_ind = grid.face_ind(bx * 8 + idx, by * 8 + idy, 1);
	if (idx == 7 && bx != nbx - 1 && (idy == 0 || by != nby - 1)) atomicAdd(face_y + face_ind, -node_data);
}

void grid2DOperator::D0Mapping(const grid2D& grid, Scalar *face_x, Scalar *face_y, Scalar *node)
{
	cudaMemset(face_x, 0, sizeof(Scalar)*grid.face_size(0));
	cudaMemset(face_y, 0, sizeof(Scalar)*grid.face_size(1));
	D0MappingKernel << <dim3((grid.Nx >> 3) + 1, (grid.Ny >> 3) + 1), dim3(8, 8) >> > (grid, face_x, face_y, node);
}

// for blockDim = (8, 8) 
// no boundary
// ugly implementation
// iterate through face
__global__ void FaceXLapMappingKernel(grid2D grid, Scalar *lap_face_x, Scalar *face_x)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ Scalar shared_face_data[64];

	// 8x8 face_x data load
	shared_face_data[idy * 8 + idx] = face_x[grid.face_ind(bx * 8 + idy, by * 8 + idx, 0)];
	__syncthreads();

	const Scalar v = shared_face_data[idx * 8 + idy];
	Scalar lap = (Scalar)4 * v;

	if (idx != 0) lap -= shared_face_data[(idx - 1) * 8 + idy];

	if (idx != 7) lap -= shared_face_data[(idx + 1) * 8 + idy];

	if (idy != 0) lap -= shared_face_data[idx * 8 + (idy - 1)];

	if (idy != 7) lap -= shared_face_data[idx * 8 + (idy + 1)];

	atomicAdd(lap_face_x + grid.face_ind(bx * 8 + idx, by * 8 + idy, 0), lap);

	if (bx != 0 && idx == 0) atomicAdd(lap_face_x + grid.face_ind((bx - 1) * 8 + 7, by * 8 + idy, 0), -v);

	if (bx != nbx - 1 && idx == 7) atomicAdd(lap_face_x + grid.face_ind((bx + 1) * 8 + 0, by * 8 + idy, 0), -v);

	if (by != 0 && idy == 0) atomicAdd(lap_face_x + grid.face_ind(bx * 8 + idx, (by - 1) * 8 + 7, 0), -v);

	if (by != nby - 1 && idy == 7) atomicAdd(lap_face_x + grid.face_ind(bx * 8 + idx, (by + 1) * 8 + 0, 0), -v);
}

// for blockDim = (8, 8) 
// no boundary
// ugly implementation
// iterate through face
__global__ void FaceYLapMappingKernel(grid2D grid, Scalar *lap_face_y, Scalar *face_y)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	__shared__ Scalar shared_face_data[64];

	// 8x8 face_y data load
	shared_face_data[idy * 8 + idx] = face_y[grid.face_ind(bx * 8 + idx, by * 8 + idy, 1)];
	__syncthreads();

	const Scalar v = shared_face_data[idy * 8 + idx];
	Scalar lap = (Scalar)4 * v;

	if (idx != 0) lap -= shared_face_data[idy * 8 + (idx - 1)];

	if (idx != 7) lap -= shared_face_data[idy * 8 + (idx + 1)];

	if (idy != 0) lap -= shared_face_data[(idy - 1) * 8 + idx];

	if (idy != 7) lap -= shared_face_data[(idy + 1) * 8 + idx];

	atomicAdd(lap_face_y + grid.face_ind(bx * 8 + idx, by * 8 + idy, 1), lap);

	if (bx != 0 && idx == 0) atomicAdd(lap_face_y + grid.face_ind((bx - 1) * 8 + 7, by * 8 + idy, 1), -v);

	if (bx != nbx - 1 && idx == 7) atomicAdd(lap_face_y + grid.face_ind((bx + 1) * 8 + 0, by * 8 + idy, 1), -v);

	if (by != 0 && idy == 0) atomicAdd(lap_face_y + grid.face_ind(bx * 8 + idx, (by - 1) * 8 + 7, 1), -v);

	if (by != nby - 1 && idy == 7) atomicAdd(lap_face_y + grid.face_ind(bx * 8 + idx, (by + 1) * 8 + 0, 1), -v);
}


void grid2DOperator::cellLapMapping(const grid2D& grid, Scalar * lap_cell, Scalar * cell)
{
	cudaMemset(lap_cell, 0, sizeof(Scalar)*grid.cell_size());
	FaceYLapMappingKernel << <dim3((grid.Nx >> 3), (grid.Ny >> 3)), dim3(8, 8) >> > (grid, lap_cell, cell);
}

void grid2DOperator::faceLapMapping(const grid2D& grid, Scalar * lap_face_x, Scalar * lap_face_y, Scalar * face_x, Scalar * face_y)
{
	cudaMemset(lap_face_x, 0, sizeof(Scalar)*grid.face_size(0));
	cudaMemset(lap_face_y, 0, sizeof(Scalar)*grid.face_size(1));
	FaceXLapMappingKernel << <dim3((grid.Nx >> 3) + 1, (grid.Ny >> 3)), dim3(8, 8) >> > (grid, lap_face_x, face_x);
	FaceYLapMappingKernel << <dim3((grid.Nx >> 3), (grid.Ny >> 3) + 1), dim3(8, 8) >> > (grid, lap_face_y, face_y);
}
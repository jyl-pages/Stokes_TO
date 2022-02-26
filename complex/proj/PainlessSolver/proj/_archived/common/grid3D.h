#pragma once

#include "para.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


class grid3D
{
public:
	int Nx, Ny, Nz;

	grid3D();

	grid3D(int _Nx, int _Ny, int _Nz);

	void init(int _Nx, int _Ny, int _Nz);

	__host__ __device__ int cell_size() const;

	__host__ __device__ int face_size(int d) const;

	__host__ __device__ int face_size_sum() const;

	__host__ __device__ int edge_size(int d) const;

	__host__ __device__ int edge_size_sum() const;

	__host__ __device__ int node_size() const;

	__host__ __device__ int cell_ind(int x, int y,int z) const;

	__host__ __device__ int face_ind(int x, int y, int z, int d) const;

	__host__ __device__ int edge_ind(int x, int y, int z, int d) const;

	__host__ __device__ int node_ind(int x, int y, int z) const;

	__host__ __device__ void ind_cell(int i, int& x, int& y, int &z) const;

	__host__ __device__ void ind_face(int i, int d, int& x, int& y, int& z) const;

	__host__ __device__ void ind_edge(int i, int d, int& x, int& y, int& z) const;

	__host__ __device__ void ind_node(int i, int& x, int& y, int& z) const;

	__host__ __device__ bool cell_valid(int x, int y, int z) const;

	__host__ __device__ bool face_valid(int x, int y, int z, int d) const;

	__host__ __device__ bool edge_valid(int x, int y, int z, int d) const;

	__host__ __device__ bool node_valid(int x, int y, int z) const;

};

namespace grid3DOperator
{
	void ApplyMask(const grid3D& grid, Scalar *v, bool *fixed, int N);

	void ApplyFaceRedudantFixed(const grid3D& grid, bool *fixed, int d);

	void ApplyEdgeRedudantFixed(const grid3D& grid, bool *fixed, int d);

	void ApplyNodeRedudantFixed(const grid3D& grid, bool *fixed);

	void D0Mapping(const grid3D& grid, Scalar *edge_x, Scalar *edge_y, Scalar *edge_z, Scalar *node);

	void D1Mapping(const grid3D& grid, Scalar *face_x, Scalar *face_y, Scalar *face_z, Scalar *edge_x, Scalar *edge_y, Scalar *edge_z);

	void D2Mapping(const grid3D& grid, Scalar *cell, Scalar *face_x, Scalar *face_y, Scalar *face_z);

	void Cod0Mapping(const grid3D& grid, Scalar *face_x, Scalar *face_y, Scalar *face_z, Scalar *cell);

	void Cod1Mapping(const grid3D& grid, Scalar *edge_x, Scalar *edge_y, Scalar *edge_z, Scalar *face_x, Scalar *face_y, Scalar *face_z);

	void Cod2Mapping(const grid3D& grid, Scalar *node, Scalar *edge_x, Scalar *edge_y, Scalar *edge_z);
};
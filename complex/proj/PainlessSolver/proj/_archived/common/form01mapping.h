#pragma once

#include "para.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"


class grid2D
{
public:
	int Nx, Ny;

	void init(int _Nx, int _Ny);

	__host__ __device__ int cell_size() const;

	__host__ __device__ int face_size(int d) const;

	__host__ __device__ int node_size() const;

	__host__ __device__ int cell_ind(int x, int y) const;

	__host__ __device__ int face_ind(int x, int y, int d) const;

	__host__ __device__ int node_ind(int x, int y) const;

	__host__ __device__ void ind_cell(int i, int& x, int& y) const;

	__host__ __device__ void ind_face(int i, int d, int& x, int& y) const;

	__host__ __device__ void ind_node(int i, int& x, int& y) const;

};

namespace grid2DOperator
{
	void ApplyMask(const grid2D& grid, Scalar *v, bool *fixed,int N);

	void ApplyFaceRedudantFixed(const grid2D& grid, bool *fixed, int d);

	void ApplyNodeRedudantFixed(const grid2D& grid, bool *fixed);

	void Cod0Mapping(const grid2D& grid, Scalar *face_x, Scalar *face_y, Scalar *cell);

	void D1Mapping(const grid2D& grid, Scalar *cell, Scalar *face_x, Scalar *face_y);

	void Cod1Mapping(const grid2D& grid, Scalar *node, Scalar *face_x, Scalar *face_y);

	void D0Mapping(const grid2D& grid, Scalar *face_x, Scalar *face_y, Scalar *node);

	void cellLapMapping(const grid2D& grid, Scalar *lap_cell, Scalar *cell);

	void faceLapMapping(const grid2D& grid, Scalar *lap_face_x, Scalar *lap_face_y, Scalar *face_x, Scalar *face_y);

};
#pragma once

#include "form01mapping.h"
#include "grid3D.h"
#include <nvfunctional>

class GridMask2D
{
public:
	grid2D grid;
	nvstd::function<bool(int, int)> f;

	GridMask2D(const grid2D& _grid, const nvstd::function<bool(int, int)>& _f)
	{
		grid = _grid;
		f = _f;
	}

	__host__ __device__ bool operator[](int i) const
	{
		int x, y;
		grid.ind_cell(i, x, y);
		return f(x, y);
	}
};
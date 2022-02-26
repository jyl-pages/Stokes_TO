#pragma once

#include "cpuUtils.h"
#include "form01mapping.h"
#include "grid3D.h"
#include "TypeFunc.h"

template<int d>
class PoissonDescriptor
{
	Typedef_VectorDii(d);
public:
	VectorDi N;
	//int Nx, Ny;
	int size;
	int fsize;
	using Grid = typename If<d == 2, grid2D, grid3D >::Type;
	Grid grid;
	Scalar *h_vol, *d_vol;//h on CPU and d on GPU
	bool *h_fixed, *d_fixed;

public:
	void init(const VectorDi& _n);
	//void init(int _Nx, int _Ny);

	void setFixed(bool *fixed);

	void setVol(Scalar *vol);

	void finish();

	void toHost();

	void toDevice();

};
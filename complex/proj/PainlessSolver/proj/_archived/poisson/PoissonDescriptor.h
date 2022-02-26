#pragma once

#include "cpuUtils.h"
#include "form01mapping.h"

class PoissonDescriptor
{
public:
	int Nx, Ny;
	int size;
	grid2D grid;
	Scalar *h_vol, *d_vol;
	bool *h_fixed, *d_fixed;

public:
	void init(int _Nx, int _Ny);

	void setFixed(bool *fixed);

	void setVol(Scalar *vol);

	void finish();

	void toHost();

	void toDevice();

};
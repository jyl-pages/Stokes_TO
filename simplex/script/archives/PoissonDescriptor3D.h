#pragma once

#include "cpuUtils.h"
#include "grid3D.h"

class PoissonDescriptor3D
{
public:
	int Nx, Ny, Nz;
	int size;
	int fsize;
	grid3D grid;
	Scalar *h_vol, *d_vol;
	bool *h_fixed, *d_fixed;

public:
	void init(int _Nx, int _Ny, int _Nz);

	void setFixed(bool *fixed);

	void setVol(Scalar *vol);

	void finish();

	void toHost();

	void toDevice();

};
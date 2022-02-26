#pragma once
#include "cpuUtils.h"
#include "form01mapping.h"

class StokeFlowDescriptor
{
public:
	int Nx, Ny;
	int size;
	Scalar lap_coeff = (Scalar)1;
	grid2D grid;
	Scalar *h_vol, *d_vol;
	Scalar *h_penalty, *d_penalty;
	bool *h_fixed, *d_fixed;

public:
	void init(int _Nx, int _Ny);

	void setFixed(bool *fixed);

	void setVol(Scalar *vol);

	void setPenalty(Scalar *penalty);

	void finish();

	void toHost();

	void toDevice();
};

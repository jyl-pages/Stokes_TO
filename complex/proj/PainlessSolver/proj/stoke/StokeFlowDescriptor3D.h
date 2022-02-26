#pragma once
#include "para.h"
#include "grid3D.h"

class StokeFlowDescriptor3D
{
public:
	int Nx, Ny, Nz;
	int size;
	int tsize;
	grid3D grid;
	Scalar *h_vol, *d_vol;
	Scalar *h_penalty, *d_penalty;
	bool *h_fixed, *d_fixed;

public:
	void init(int _Nx, int _Ny, int _Nz);

	void setFixed(bool *fixed);

	void setVol(Scalar *vol);

	void setPenalty(Scalar *penalty);

	void finish();

	void toHost();

	void toDevice();
};

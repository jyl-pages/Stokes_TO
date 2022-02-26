#pragma once
#include "LinearMapping.h"
#include "form01mapping.h"
#include "PoissonDescriptor.h"

class PoissonMapping :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	grid2D grid;
	PoissonDescriptor *descr;

	Scalar *d_temp;

	void init(PoissonDescriptor *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class PoissonMappingFixed :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	grid2D grid;
	PoissonDescriptor *descr;

	Scalar *d_temp;

	void init(PoissonDescriptor *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};
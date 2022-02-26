#pragma once
#include "LinearMapping.h"
#include "grid3D.h"
#include "PoissonDescriptor.h"

class PoissonMapping3D :public LinearMapping
{
public:
	int Nx, Ny, Nz;
	int dof;
	grid3D grid;
	PoissonDescriptor<3> *descr;

	Scalar *d_temp;

	void init(PoissonDescriptor<3> *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class PoissonMapping3DFixed :public LinearMapping
{
public:
	int Nx, Ny, Nz;
	int dof;
	grid3D grid;
	PoissonDescriptor<3> *descr;

	Scalar *d_temp;

	void init(PoissonDescriptor<3> *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

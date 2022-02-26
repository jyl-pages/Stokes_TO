#pragma once
#include "para.h"
#include "LinearMapping.h"
#include "PoissonDescriptor.h"

class PoissonDownSample :public LinearMapping
{
public:
	int f_Nx, f_Ny;
	int c_Nx, c_Ny;
	int xdof, ydof;
	PoissonDescriptor *f_descr;
	PoissonDescriptor *c_descr;

	PoissonDownSample(PoissonDescriptor *_f_descr, PoissonDescriptor *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class PoissonUpSample :public LinearMapping
{
public:
	int f_Nx, f_Ny;
	int c_Nx, c_Ny;
	int xdof, ydof;
	PoissonDescriptor *f_descr;
	PoissonDescriptor *c_descr;

	PoissonUpSample(PoissonDescriptor *_f_descr, PoissonDescriptor *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

void updateSubSystem(PoissonDescriptor& c_descr, const PoissonDescriptor& f_descr);

PoissonDescriptor createSubSystem(const PoissonDescriptor& f_descr);
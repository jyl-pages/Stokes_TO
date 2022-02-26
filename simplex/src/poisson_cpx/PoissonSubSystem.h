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
	PoissonDescriptor<2> *f_descr;
	PoissonDescriptor<2> *c_descr;

	PoissonDownSample(PoissonDescriptor<2> *_f_descr, PoissonDescriptor<2> *_c_descr);

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
	PoissonDescriptor<2>*f_descr;
	PoissonDescriptor<2>*c_descr;

	PoissonUpSample(PoissonDescriptor<2>*_f_descr, PoissonDescriptor<2>*_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

void updateSubSystem(PoissonDescriptor<2>& c_descr, const PoissonDescriptor<2>& f_descr);

PoissonDescriptor<2> createSubSystem(const PoissonDescriptor<2>& f_descr);
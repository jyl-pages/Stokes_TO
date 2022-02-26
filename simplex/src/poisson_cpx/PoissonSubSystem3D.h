#pragma once
#include "para.h"
#include "LinearMapping.h"
#include "PoissonDescriptor.h"

class PoissonDownSample3D :public LinearMapping
{
public:
	int xdof, ydof;
	PoissonDescriptor<3> *f_descr;
	PoissonDescriptor<3>*c_descr;

	PoissonDownSample3D(PoissonDescriptor<3>*_f_descr, PoissonDescriptor<3>*_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class PoissonUpSample3D :public LinearMapping
{
public:
	int xdof, ydof;
	PoissonDescriptor<3> *f_descr;
	PoissonDescriptor<3> *c_descr;

	PoissonUpSample3D(PoissonDescriptor<3> *_f_descr, PoissonDescriptor<3> *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

void updateSubSystem(PoissonDescriptor<3>& c_descr, const PoissonDescriptor<3>& f_descr);

PoissonDescriptor<3> createSubSystem(const PoissonDescriptor<3>& f_descr);
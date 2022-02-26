#pragma once
#include "para.h"
#include "LinearMapping.h"
#include "PoissonDescriptor3D.h"

class PoissonDownSample3D :public LinearMapping
{
public:
	int xdof, ydof;
	PoissonDescriptor3D *f_descr;
	PoissonDescriptor3D *c_descr;

	PoissonDownSample3D(PoissonDescriptor3D *_f_descr, PoissonDescriptor3D *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class PoissonUpSample3D :public LinearMapping
{
public:
	int xdof, ydof;
	PoissonDescriptor3D *f_descr;
	PoissonDescriptor3D *c_descr;

	PoissonUpSample3D(PoissonDescriptor3D *_f_descr, PoissonDescriptor3D *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

void updateSubSystem(PoissonDescriptor3D& c_descr, const PoissonDescriptor3D& f_descr);

PoissonDescriptor3D createSubSystem(const PoissonDescriptor3D& f_descr);
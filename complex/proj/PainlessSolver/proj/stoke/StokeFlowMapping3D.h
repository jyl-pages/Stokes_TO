#pragma once
#include "LinearMapping.h"
#include "grid3D.h"
#include "StokeFlowDescriptor3D.h"

class StokeFlowMapping3D :public LinearMapping
{
public:
	int Nx, Ny,Nz;
	int dof;
	grid3D grid;
	StokeFlowDescriptor3D *descr;

	Scalar *d_temp;

	void init(StokeFlowDescriptor3D *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class StokeFlowMapping3DFixed :public LinearMapping
{
public:
	int Nx, Ny, Nz;
	int dof;
	grid3D grid;
	StokeFlowDescriptor3D *descr;

	Scalar *d_temp;

	void init(StokeFlowDescriptor3D *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

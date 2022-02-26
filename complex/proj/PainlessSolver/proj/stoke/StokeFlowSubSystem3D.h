#pragma once
#include "para.h"
#include "LinearMapping.h"
#include "StokeFlowDescriptor3D.h"

class StokeFlowDownSample3D :public LinearMapping
{
public:
	int f_Nx, f_Ny, f_Nz;
	int c_Nx, c_Ny, c_Nz;
	int xdof, ydof;
	StokeFlowDescriptor3D *f_descr;
	StokeFlowDescriptor3D *c_descr;

	void init(StokeFlowDescriptor3D *_f_descr, StokeFlowDescriptor3D *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class StokeFlowUpSample3D :public LinearMapping
{
public:
	int f_Nx, f_Ny, f_Nz;
	int c_Nx, c_Ny, c_Nz;
	int xdof, ydof;
	StokeFlowDescriptor3D *f_descr;
	StokeFlowDescriptor3D *c_descr;

	void init(StokeFlowDescriptor3D *_f_descr, StokeFlowDescriptor3D *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

void updateSubSystem(StokeFlowDescriptor3D& c_descr, const StokeFlowDescriptor3D& f_descr);

StokeFlowDescriptor3D createSubSystem(const StokeFlowDescriptor3D& f_descr);
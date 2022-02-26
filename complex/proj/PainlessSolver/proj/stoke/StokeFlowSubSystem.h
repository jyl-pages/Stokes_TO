#pragma once
#include "para.h"
#include "LinearMapping.h"
#include "StokeFlowDescriptor.h"

class StokeFlowDownSample :public LinearMapping
{
public:
	int f_Nx, f_Ny;
	int c_Nx, c_Ny;
	int xdof, ydof;
	StokeFlowDescriptor *f_descr;
	StokeFlowDescriptor *c_descr;

	void init(StokeFlowDescriptor *_f_descr, StokeFlowDescriptor *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

class StokeFlowUpSample :public LinearMapping
{
public:
	int f_Nx, f_Ny;
	int c_Nx, c_Ny;
	int xdof, ydof;
	StokeFlowDescriptor *f_descr;
	StokeFlowDescriptor *c_descr;

	void init(StokeFlowDescriptor *_f_descr, StokeFlowDescriptor *_c_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

void updateSubSystem(StokeFlowDescriptor& c_descr, const StokeFlowDescriptor& f_descr);

StokeFlowDescriptor createSubSystem(const StokeFlowDescriptor& f_descr);
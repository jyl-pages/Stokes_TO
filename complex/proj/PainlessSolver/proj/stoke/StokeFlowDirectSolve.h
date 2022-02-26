#pragma once
#include "LinearMapping.h"
#include "StokeFlowDescriptor.h"
#include "cusolverDn.h"

class StokeFlowDirectSolve :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	StokeFlowDescriptor *descr;
	Scalar *h_A, *d_A;
	int *d_piv;
	Scalar *d_buffer;
	int h_info, *d_info;
	cusolverDnHandle_t cusolverDnHandle;

	bool closed = false;

	void init(StokeFlowDescriptor* _descr);

	void finish();

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;

};
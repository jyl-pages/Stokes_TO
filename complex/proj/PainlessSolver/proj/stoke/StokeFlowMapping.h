#pragma once
#include "LinearMapping.h"
#include "form01mapping.h"
#include "StokeFlowDescriptor.h"
#include "cublas_v2.h"

class StokeFlowMapping :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	grid2D grid;
	StokeFlowDescriptor *descr;

	Scalar *d_temp;

	cublasHandle_t cublasHandle;

	void init(StokeFlowDescriptor *_descr);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};
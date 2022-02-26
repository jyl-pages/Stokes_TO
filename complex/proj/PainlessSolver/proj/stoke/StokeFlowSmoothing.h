#pragma once
#include "LinearMapping.h"
#include "StokeFlowDescriptor.h"
#include "StokeFlowMapping.h"
#include "cublas_v2.h"

class StokeFlowSmoothing :public LinearMapping
{
public:
	int Nx, Ny;
	int dof;
	StokeFlowDescriptor *descr;
	StokeFlowMapping *stoke_flow_mapping;
	Scalar *d_p, *d_Ap;
	Scalar *d_r;
	int order = 0;
	int iter_num = 1;
	Scalar alpha = (Scalar)1;

	cublasHandle_t cublasHandle;

	void init(StokeFlowMapping *_stoke_flow_mapping);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;
};

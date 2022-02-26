#pragma once

#include "LinearMapping.h"
#include "cublas_v2.h"

class MultigridSolve :public LinearMapping
{
public:
	int l;
	int *dof;
	LinearMapping **mapping;
	LinearMapping **downSample;
	LinearMapping **upSample;
	LinearMapping **preSmoother;
	LinearMapping **postSmoother;

	Scalar **r;
	Scalar **x;
	Scalar **p;
	Scalar **Ap;

	cublasHandle_t cublasHandle;

public:
	void init(int _l);

	void finish();

	int xDoF() override;

	int yDoF() override;

	void Vcycle();

	void applyMapping(Scalar *Ap, Scalar *p) override;

	void update() override;
};

#pragma once
#include "LinearMapping.h"
#include "Eigen/Dense"
#include "cublas_v2.h"

class ConjugatedGradient
{
public:
	LinearMapping *linear_mapping;
	LinearMapping *preconditioner = nullptr;
	int dof;
	Scalar *d_p, *d_Ap, *d_z, *d_r, *d_x;

	Scalar CG_absolute_residual_thres = (Scalar)1e-3;
	Scalar CG_relative_residual_thres = (Scalar)1e-4;

	bool verbose = true;
	int max_iter;

	cublasHandle_t cublasHandle;

	void Init(int _max_iter);

	void Solve(Scalar *x, Scalar const *b);

};
#pragma once
#include "LinearMapping.h"
#include "Eigen/Dense"
#include "cpuUtils.h"

class GMRES
{
public:
	LinearMapping *linear_mapping;
	LinearMapping *preconditioner = nullptr;
	int dof;
	VECTOR x, b;
	Scalar *d_q_buf, *d_Aq_buf;

	Scalar GMRES_norm_thres = (Scalar)1e-3;
	Scalar GMRES_absolute_residual_thres = (Scalar)1e-3;
	Scalar GMRES_relative_residual_thres = (Scalar)1e-4;

	bool debug = false;
	bool verbose = true;
	int max_iter;
	void Init(int _max_iter);

	void Solve();
};
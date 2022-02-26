#pragma once

#include "LinearMapping.h"
#include <vector>

class DescentSmoother :public LinearMapping
{
public:
	int dof;
	std::vector<LinearMapping*> P;
	LinearMapping *A;
	std::vector<int> iters;

	Scalar *d_r, *d_p, *d_Ap;

public:
	DescentSmoother(std::vector<LinearMapping*> _P, std::vector<int> _iters, LinearMapping *_A);

	int xDoF() override;

	int yDoF() override;

	void applyMapping(Scalar *Ap, Scalar *p) override;

	void update() override;
};
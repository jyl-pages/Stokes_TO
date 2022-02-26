#include "DescentSmoother.h"
#include <assert.h>
#include "gpuUtils.h"
#include "cuda_runtime_api.h"

DescentSmoother::DescentSmoother(std::vector<LinearMapping*> _P, std::vector<int> _iters, LinearMapping * _A)
{
	assert(_A.xDoF() == _A.yDoF());
	assert(_P.xDoF() == _P.yDoF());
	assert(_A.xDoF() == _P.xDoF());
	dof = _A->xDoF();
	P = _P;
	A = _A;
	iters = _iters;

	cudaMalloc(&d_r, sizeof(Scalar)*dof);
	cudaMalloc(&d_p, sizeof(Scalar)*dof);
	cudaMalloc(&d_Ap, sizeof(Scalar)*dof);
}

int DescentSmoother::xDoF()
{
	return dof;
}

int DescentSmoother::yDoF()
{
	return dof;
}

void DescentSmoother::applyMapping(Scalar *Ap, Scalar *p)
{
	auto add = [=] __device__(Scalar& tv1, Scalar tv2) { tv1 += tv2; };
	auto sub = [=] __device__(Scalar& tv1, Scalar tv2) { tv1 -= tv2; };

	cudaMemcpy(d_r, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
	cudaMemset(Ap, 0, sizeof(Scalar)*dof);
	for (int i = 0; i < P.size(); i++)
		for (int j = 0; j < iters[i]; j++)
		{
			P[i]->applyMapping(d_p, d_r);
			cwise_mapping_wrapper(Ap, d_p, add, dof);
			if (i == P.size() - 1 && j == iters[i] - 1) break;
			A->applyMapping(d_Ap, d_p);
			cwise_mapping_wrapper(d_r, d_Ap, sub, dof);
		}
}

void DescentSmoother::update()
{
	for (LinearMapping *p : P)
		p->update();
	A->update();
}
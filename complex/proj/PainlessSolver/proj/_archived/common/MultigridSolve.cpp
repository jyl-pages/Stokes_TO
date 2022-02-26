#include "MultigridSolve.h"
#include "cuda_runtime_api.h"
#include "gpuUtils.h"
#include <assert.h>

void MultigridSolve::init(int _l)
{
	l = _l;

	mapping = (LinearMapping**)malloc(sizeof(LinearMapping*)*l);
	downSample = (LinearMapping**)malloc(sizeof(LinearMapping*)*(l - 1));
	upSample = (LinearMapping**)malloc(sizeof(LinearMapping*)*(l - 1));
	preSmoother = (LinearMapping**)malloc(sizeof(LinearMapping*)*l);
	postSmoother= (LinearMapping**)malloc(sizeof(LinearMapping*)*l);

	r = (Scalar**)malloc(sizeof(Scalar*)*l);
	x = (Scalar**)malloc(sizeof(Scalar*)*l);
	p = (Scalar**)malloc(sizeof(Scalar*)*l);
	Ap = (Scalar**)malloc(sizeof(Scalar*)*l);

	cublasCreate(&cublasHandle);
}

void MultigridSolve::finish()
{
	dof = (int*)malloc(sizeof(int)*l);
	for (int i = 0; i < l; i++)
	{
		assert(mapping[i]->xDoF() == mapping[i]->yDoF());
		dof[i] = mapping[i]->xDoF();
	}
	for (int i = 0; i < l - 1; i++)
	{
		assert(downSample[i]->xDoF() == dof[i + 1]);
		assert(downSample[i]->yDoF() == dof[i]);
		assert(upSample[i]->xDoF() == dof[i]);
		assert(upSample[i]->yDoF() == dof[i + 1]);
	}
	for (int i = 0; i < l; i++)
	{
		if (preSmoother[i] != nullptr)
		{
			assert(preSmoother[i]->xDoF() == dof[i]);
			assert(preSmoother[i]->yDoF() == dof[i]);
		}
		if (postSmoother[i] != nullptr)
		{
			assert(postSmoother[i]->xDoF() == dof[i]);
			assert(postSmoother[i]->yDoF() == dof[i]);
		}
	}

	for (int i = 0; i < l; i++)
	{
		cudaMalloc((void**)(r + i), sizeof(Scalar)*dof[i]);
		cudaMalloc((void**)(x + i), sizeof(Scalar)*dof[i]);
		cudaMalloc((void**)(p + i), sizeof(Scalar)*dof[i]);
		cudaMalloc((void**)(Ap + i), sizeof(Scalar)*dof[i]);
	}
}


int MultigridSolve::xDoF()
{
	return dof[l - 1];
}

int MultigridSolve::yDoF()
{
	return dof[l - 1];
}

void MultigridSolve::Vcycle()
{
	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;

	for (int i = l - 1; i >= 0; i--)
	{
		cudaMemset(p[i], 0, sizeof(Scalar)*dof[i]);
		if (preSmoother[i]) preSmoother[i]->applyMapping(p[i], r[i]);
		mapping[i]->applyMapping(Ap[i], p[i]);
		Axpy(cublasHandle, dof[i], &one, p[i], 1, x[i], 1);
		Axpy(cublasHandle, dof[i], &neg_one, Ap[i], 1, r[i], 1);
		if (i != 0)
		{
			downSample[i - 1]->applyMapping(r[i - 1], r[i]);
			cudaMemset(x[i - 1], 0, sizeof(Scalar)*dof[i - 1]);
		}
	}

	for (int i = 0; i < l; i++)
	{
		if (i != 0)
		{
			upSample[i - 1]->applyMapping(p[i], x[i - 1]);
			mapping[i]->applyMapping(Ap[i], p[i]);
			Axpy(cublasHandle, dof[i], &one, p[i], 1, x[i], 1);
			Axpy(cublasHandle, dof[i], &neg_one, Ap[i], 1, r[i], 1);
		}
		cudaMemset(p[i], 0, sizeof(Scalar)*dof[i]);
		if (postSmoother[i]) postSmoother[i]->applyMapping(p[i], r[i]);
		mapping[i]->applyMapping(Ap[i], p[i]);
		Axpy(cublasHandle, dof[i], &one, p[i], 1, x[i], 1);
		Axpy(cublasHandle, dof[i], &neg_one, Ap[i], 1, r[i], 1);
	}
}

void MultigridSolve::applyMapping(Scalar *Ap, Scalar *p)
{
	cudaMemcpy(r[l - 1], p, sizeof(Scalar)*dof[l - 1], cudaMemcpyDeviceToDevice);
	cudaMemset(x[l - 1], 0, sizeof(Scalar)*dof[l - 1]);
	Vcycle();
	cudaMemcpy(Ap, x[l - 1], sizeof(Scalar)*dof[l - 1], cudaMemcpyDeviceToDevice);
}

void MultigridSolve::update()
{
	for (int i = 0; i < l; i++)
	{
		if(mapping[i]) mapping[i]->update();
		if(preSmoother[i]) preSmoother[i]->update();
		if(postSmoother[i]) postSmoother[i]->update();
	}
}

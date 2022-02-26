#include "ConjugatedGradient.h"
#include <memory>
#include <iostream>
#include "cuda_runtime_api.h"
#include "gpuUtils.h"

void ConjugatedGradient::Init(int _max_iter)
{
	max_iter = _max_iter;
	assert(linear_mapping->xDoF() == linear_mapping->yDoF());
	if (preconditioner)
	{
		assert(preconditioner->xDoF() == preconditioner->yDoF());
		assert(linear_mapping->xDoF() == linear_mapping->xDoF());
	}
	dof = linear_mapping->xDoF();
	cudaMalloc((void**)&d_p, sizeof(Scalar)*dof);
	cudaMalloc((void**)&d_Ap, sizeof(Scalar)*dof);
	cudaMalloc((void**)&d_z, sizeof(Scalar)*dof);
	cudaMalloc((void**)&d_r, sizeof(Scalar)*dof);
	cudaMalloc((void**)&d_x, sizeof(Scalar)*dof);

	cublasCreate(&cublasHandle);
}

void ConjugatedGradient::Solve(Scalar *x, Scalar const *b)
{
	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;
	Scalar norm0, old_norm, new_norm, dot;

	cudaMemcpy(d_r, b, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
	cudaMemset(d_x, 0, sizeof(Scalar)*dof);

	if (preconditioner) preconditioner->applyMapping(d_z, d_r);
	else cudaMemcpy(d_z, d_r, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

	cudaMemcpy(d_p, d_z, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

	Dot(cublasHandle, dof, d_z, 1, d_r, 1, &old_norm);
	cudaDeviceSynchronize();
	norm0 = old_norm;

	for (int i = 0; i < max_iter; i++)
	{
		if (verbose) printf("iter: %d res:%lf\n", i, old_norm);

		linear_mapping->applyMapping(d_Ap, d_p);

		Dot(cublasHandle, dof, d_p, 1, d_Ap, 1, &dot);
		cudaDeviceSynchronize();

		Scalar alpha = old_norm / dot;
		Scalar neg_alpha = -alpha;

		Axpy(cublasHandle, dof, &alpha, d_p, 1, d_x, 1);

		Axpy(cublasHandle, dof, &neg_alpha, d_Ap, 1, d_r, 1);

		Dot(cublasHandle, dof, d_r, 1, d_r, 1, &new_norm);
		cudaDeviceSynchronize();
		if (new_norm < CG_absolute_residual_thres)
			break;

		if (preconditioner) preconditioner->applyMapping(d_z, d_r);
		else cudaMemcpy(d_z, d_r, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

		Dot(cublasHandle, dof, d_z, 1, d_r, 1, &new_norm);
		cudaDeviceSynchronize();


		if (new_norm < CG_relative_residual_thres*norm0)
			break;

		Scalar beta = new_norm / old_norm;

		Scal(cublasHandle, dof, &beta, d_p, 1);
		Axpy(cublasHandle, dof, &one, d_z, 1, d_p, 1);

		old_norm = new_norm;
	}

	cudaMemcpy(x, d_x, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();
}

#include "GMRES.h"
#include <memory>
#include <iostream>
#include "cuda_runtime.h"

void GMRES::Init(int _max_iter)
{
	max_iter = _max_iter;
	assert(linear_mapping->xDoF() == linear_mapping->yDoF());
	if (preconditioner)
	{
		assert(preconditioner->xDoF() == preconditioner->yDoF());
		assert(linear_mapping->xDoF() == linear_mapping->xDoF());
	}
	dof = linear_mapping->xDoF();
	cudaMalloc(&d_q_buf, sizeof(Scalar)*dof);
	cudaMalloc(&d_Aq_buf, sizeof(Scalar)*dof);
}

void GMRES::Solve()
{
	VECTOR Aq(dof);
	MATRIX Q(dof, 1);
	MATRIX H(1, 0);
	VECTOR e0(1);

	//Q.col(0) = b.normalized();
	//e0(0) = b.norm();
	cudaMemcpy(d_Aq_buf, b.data(), sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
	if(preconditioner)
	{
		cudaMemcpy(d_q_buf, d_Aq_buf, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
		preconditioner->applyMapping(d_Aq_buf, d_q_buf);
	}
	cudaMemcpy(Q.col(0).data(), d_Aq_buf, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
	cudaDeviceSynchronize();

	//std::cout << Q.col(0) << std::endl;

	e0(0) = Q.col(0).norm();
	Q.col(0).normalize();

	if(verbose) std::cout << "e0:\n" << e0 << std::endl;

	bool stop_flag = false;
	for (int i = 1; i < max_iter; i++)
	{
		Q.conservativeResize(Eigen::NoChange, i + 1);
		H.conservativeResize(i + 1, i);
		for (int j = 0; j < i - 1; j++) for (int k = j + 2; k < i + 1; k++) H(k, j) = (Scalar)0;
		e0.conservativeResize(i + 1);
		e0(i) = (Scalar)0;

		cudaMemcpy(d_q_buf, Q.col(i - 1).data(), sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
		linear_mapping->applyMapping(d_Aq_buf, d_q_buf);
		if (preconditioner)
		{
			cudaMemcpy(d_q_buf, d_Aq_buf, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);
			preconditioner->applyMapping(d_Aq_buf, d_q_buf);
		}
		cudaMemcpy(Aq.data(), d_Aq_buf, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();

		for (int j = 0; j < i; j++)
		{
			Scalar dot = Aq.dot(Q.col(j));
			H(j, i - 1) = dot;
		}

		{
			for (int j = 0; j < i; j++)
				Aq -= H(j, i - 1) * Q.col(j);
			if (Aq.norm() < GMRES_norm_thres) stop_flag = true;
			
			H(i, i - 1) = Aq.norm();
			Q.col(i) = Aq.normalized();
		}

		VECTOR y = H.colPivHouseholderQr().solve(e0);

		//std::cout << "H:\n" << H << std::endl;
		//std::cout << "y:\n" << y << std::endl;

		x = Q.leftCols(i) *y;

		Scalar residual = (e0 - H * y).norm();

		if(verbose) printf("iter: %d Aq norm: %lf\n", i, H(i, i - 1));
		if(verbose) printf("iter: %d res: %lf\n", i, residual);

		if (residual < GMRES_absolute_residual_thres) stop_flag = true;
		if (residual < GMRES_relative_residual_thres*e0(0)) stop_flag = true;

		if (stop_flag) break;
	}

	if (verbose)
	{
		cudaMemcpy(d_q_buf, x.data(), sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
		linear_mapping->applyMapping(d_Aq_buf, d_q_buf);
		cudaMemcpy(Aq.data(), d_Aq_buf, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
		cudaDeviceSynchronize();
		printf("final absolute res:%lf\n", (b - Aq).norm());
		printf("final relative res:%lf\n", (b - Aq).norm() / b.norm());
	}
	//x = b - Aq;
	//x = Q.col(0)*e0(0);
	//x = Q.col(1);
	//x = Aq;
}

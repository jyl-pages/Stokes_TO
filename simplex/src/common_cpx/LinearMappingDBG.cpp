#include "LinearMappingDBG.h"
#include <cstdio>
#include <memory>
#include <assert.h>
#include "cuda_runtime.h"
#include <vector>

namespace LinearMappingDBG
{
	void LinearMapping2DenseMatrix(char *path, LinearMapping *mapping)
	{
		Scalar *h_p, *d_p;
		Scalar *h_Ap, *d_Ap;
		//assert(mapping->xDoF() == mapping->yDoF());
		int dof = mapping->xDoF();
		h_p = (Scalar*)malloc(sizeof(Scalar)*dof);
		h_Ap = (Scalar*)malloc(sizeof(Scalar)*dof);
		cudaMalloc(&d_p, sizeof(Scalar)*dof);
		cudaMalloc(&d_Ap, sizeof(Scalar)*dof);

		Scalar *A;
		A = (Scalar*)malloc(sizeof(Scalar)*dof*dof);

		for (int i = 0; i < dof; i++)
		{
			memset(h_p, 0, sizeof(Scalar)*dof);
			h_p[i] = (Scalar)1;
			cudaMemcpy(d_p, h_p, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
			mapping->applyMapping(d_Ap, d_p);
			cudaMemcpy(h_Ap, d_Ap, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();

			for (int j = 0; j < dof; j++)
				A[i*dof + j] = h_Ap[j];
		}

		FILE *file = fopen(path, "w");

		for (int i = 0; i < dof; i++)
			for (int j = 0; j < dof; j++)
			{
				fprintf(file, "%f", A[j*dof + i]);
				if (j == dof - 1)
					fprintf(file, "\n");
				else fprintf(file, ", ");
			}

		fclose(file);

		free(h_p);
		free(h_Ap);
		cudaFree(d_p);
		cudaFree(d_Ap);
	}

	void LinearMapping2SparseMatrix(char *path, LinearMapping *mapping, Scalar thres)
	{
		Scalar *h_p, *d_p;
		Scalar *h_Ap, *d_Ap;
		//assert(mapping->xDoF() == mapping->yDoF());
		int dof = mapping->xDoF();
		h_p = (Scalar*)malloc(sizeof(Scalar)*dof);
		h_Ap = (Scalar*)malloc(sizeof(Scalar)*dof);
		cudaMalloc(&d_p, sizeof(Scalar)*dof);
		cudaMalloc(&d_Ap, sizeof(Scalar)*dof);

		
		std::vector<int> rowInd;
		std::vector<int> colInd;
		std::vector<Scalar> val;
		for (int i = 0; i < dof; i++)
		{
			memset(h_p, 0, sizeof(Scalar)*dof);
			h_p[i] = (Scalar)1;
			cudaMemcpy(d_p, h_p, sizeof(Scalar)*dof, cudaMemcpyHostToDevice);
			mapping->applyMapping(d_Ap, d_p);
			cudaMemcpy(h_Ap, d_Ap, sizeof(Scalar)*dof, cudaMemcpyDeviceToHost);
			cudaDeviceSynchronize();

			for (int j = 0; j < dof; j++)
			{
				if (fabs(h_Ap[j]) > thres)
				{
					rowInd.push_back(j);
					colInd.push_back(i);
					val.push_back(h_Ap[j]);
				}
			}
		}

		FILE *file = fopen(path, "w");

		fprintf(file, "%d %d %d\n", dof, dof, val.size());

		for (int i = 0; i < val.size(); i++)
			fprintf(file, "%d %d %f\n", rowInd[i], colInd[i], val[i]);

		fclose(file);

		free(h_p);
		free(h_Ap);
		cudaFree(d_p);
		cudaFree(d_Ap);
	}
}
#include "StokeFlowMapping.h"
#include "gpuUtils.h"
#include <cstdlib>
#include <memory.h>

void StokeFlowMapping::init(StokeFlowDescriptor *_descr)
{
	descr = _descr;
	Nx = descr->Nx;
	Ny = descr->Ny;

	grid.init(Nx, Ny);

	dof = grid.cell_size() + grid.face_size(0) + grid.face_size(1);
	cudaMalloc(&d_temp, sizeof(Scalar)*(dof + grid.node_size()));

	cublasCreate(&cublasHandle);
}

int StokeFlowMapping::xDoF()
{
	return dof;
}

int StokeFlowMapping::yDoF()
{
	return dof;
}

void StokeFlowMapping::applyMapping(Scalar *Ap, Scalar *p)
{
	int cell_off = 0;
	int face_x_off = cell_off + grid.cell_size();
	int face_y_off = face_x_off + grid.face_size(0);
	int node_off = face_y_off + grid.face_size(1);

	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;
	Scalar mu = descr->lap_coeff*(Scalar)1;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [=] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [=] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };
	auto inv = [=]__device__(Scalar& tv, Scalar tvol) { tv /= tvol; };
	auto neg = [=]__device__(Scalar& tv) { tv = -tv; };

	//// lap v
	//{
	//	grid2DOperator::faceLapMapping(grid, d_temp + face_x_off, d_temp + face_y_off, p + face_x_off, p + face_y_off);

	//	cwise_mapping_wrapper(d_temp + face_x_off, descr->d_fixed + face_x_off, fix, grid.face_size(0));
	//	cwise_mapping_wrapper(d_temp + face_y_off, descr->d_fixed + face_y_off, fix, grid.face_size(1));

	//	Axpy(cublasHandle, grid.face_size(0), &mu, d_temp + face_x_off, 1, Ap + face_x_off, 1);
	//	Axpy(cublasHandle, grid.face_size(1), &mu, d_temp + face_y_off, 1, Ap + face_y_off, 1);
	//}

	// grad div v
	{
		cudaMemcpy(d_temp + face_x_off, p + face_x_off, sizeof(Scalar)*grid.face_size(0), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_temp + face_y_off, p + face_y_off, sizeof(Scalar)*grid.face_size(1), cudaMemcpyDeviceToDevice);

		cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

		grid2DOperator::D1Mapping(grid, d_temp + cell_off, d_temp + face_x_off, d_temp + face_y_off);

		//cwise_mapping_wrapper(d_temp + cell_off, descr->d_vol + cell_off, inv, grid.cell_size());
		cwise_mapping_wrapper(d_temp + cell_off, descr->d_fixed + cell_off, fix, grid.cell_size());

		cwise_mapping_wrapper(d_temp + cell_off, neg, grid.cell_size());


		grid2DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + cell_off);

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_fixed + face_x_off, fix, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_fixed + face_y_off, fix, grid.face_size(1));

		Axpy(cublasHandle, grid.face_size(0), &one, d_temp + face_x_off, 1, Ap + face_x_off, 1);
		Axpy(cublasHandle, grid.face_size(1), &one, d_temp + face_y_off, 1, Ap + face_y_off, 1);
	}

	// - curl curl v
	{
		grid2DOperator::Cod1Mapping(grid, d_temp + node_off, p + face_x_off, p + face_y_off);

		//cwise_mapping_wrapper(d_temp + node_off, descr->d_vol + cell_off, multi, grid.cell_size());

		grid2DOperator::D0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + node_off);

		//cwise_mapping_wrapper(d_temp + face_x_off, neg, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

		Axpy(cublasHandle, grid.face_size(0), &one, d_temp + face_x_off, 1, Ap + face_x_off, 1);
		Axpy(cublasHandle, grid.face_size(1), &one, d_temp + face_y_off, 1, Ap + face_y_off, 1);
	}

	// -grad p
	{
		grid2DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, p + cell_off);

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

		Axpy(cublasHandle, grid.face_size(0), &neg_one, d_temp + face_x_off, 1, Ap + face_x_off, 1);
		Axpy(cublasHandle, grid.face_size(1), &neg_one, d_temp + face_y_off, 1, Ap + face_y_off, 1);
	}

	// div v
	{
		cudaMemcpy(d_temp + face_x_off, p + face_x_off, sizeof(Scalar)*grid.face_size(0), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_temp + face_y_off, p + face_y_off, sizeof(Scalar)*grid.face_size(1), cudaMemcpyDeviceToDevice);

		cwise_mapping_wrapper(d_temp + face_y_off, neg, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));

		grid2DOperator::D1Mapping(grid, d_temp + cell_off, d_temp + face_x_off, d_temp + face_y_off);

		cwise_mapping_wrapper(d_temp + cell_off, descr->d_fixed + cell_off, fix, grid.cell_size());

		//cwise_mapping_wrapper(d_temp + cell_off, neg, grid.cell_size());

		Axpy(cublasHandle, grid.cell_size(), &one, d_temp + cell_off, 1, Ap + cell_off, 1);
	}

	// penalty
	{
		cudaMemcpy(d_temp, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

		cwise_mapping_wrapper(d_temp, descr->d_penalty, multi, dof);

		Axpy(cublasHandle, dof, &one, d_temp, 1, Ap, 1);
		//Axpy(cublasHandle, dof, &neg_one, d_temp, 1, Ap, 1);
	}
}

#include "FluidEnergyGrad3D.h"
#include "gpuUtils.h"

void FluidEnergyGrad3D::init(StokeFlowDescriptor3D *_descr)
{
	descr = _descr;
	Nx = descr->Nx;
	Ny = descr->Ny;
	Nz = descr->Nz;

	grid.init(Nx, Ny, Nz);

	dof = grid.cell_size() + grid.face_size(0) + grid.face_size(1) + grid.face_size(2);
	cudaMalloc(&d_temp, sizeof(Scalar)*(dof + grid.edge_size_sum()));
}

int FluidEnergyGrad3D::xDoF()
{
	return dof;
}

int FluidEnergyGrad3D::yDoF()
{
	return dof;
}

void FluidEnergyGrad3D::applyMapping(Scalar *Ap, Scalar *p)
{
	int cell_off = 0;

	int face_x_off = cell_off + grid.cell_size();
	int face_y_off = face_x_off + grid.face_size(0);
	int face_z_off = face_y_off + grid.face_size(1);

	int edge_x_off = face_z_off + grid.face_size(2);
	int edge_y_off = edge_x_off + grid.edge_size(0);
	int edge_z_off = edge_y_off + grid.edge_size(1);

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);

	auto fix = [] __device__(Scalar& tv, bool tfixed) { if (tfixed) tv = (Scalar)0; };
	auto multi = [] __device__(Scalar& tv, Scalar tvol) { tv *= tvol; };
	auto add = []__device__(Scalar& tv1, Scalar tv2) { tv1 += tv2; };
	auto sub = []__device__(Scalar& tv1, Scalar tv2) { tv1 -= tv2; };

	// grad div v
	{
		cudaMemcpy(d_temp + face_x_off, p + face_x_off, sizeof(Scalar)*grid.face_size(0), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_temp + face_y_off, p + face_y_off, sizeof(Scalar)*grid.face_size(1), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_temp + face_z_off, p + face_z_off, sizeof(Scalar)*grid.face_size(2), cudaMemcpyDeviceToDevice);

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_z_off, descr->d_vol + face_z_off, multi, grid.face_size(2));

		grid3DOperator::D2Mapping(grid, d_temp + cell_off, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off);

		cwise_mapping_wrapper(d_temp + cell_off, descr->d_fixed + cell_off, fix, grid.cell_size());

		grid3DOperator::Cod0Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off, d_temp + cell_off);

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_fixed + face_x_off, fix, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_fixed + face_y_off, fix, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_z_off, descr->d_fixed + face_z_off, fix, grid.face_size(2));

		cwise_mapping_wrapper(Ap + face_x_off, d_temp + face_x_off, sub, grid.face_size(0));
		cwise_mapping_wrapper(Ap + face_y_off, d_temp + face_y_off, sub, grid.face_size(1));
		cwise_mapping_wrapper(Ap + face_z_off, d_temp + face_z_off, sub, grid.face_size(2));

	}

	// - curl curl v
	{
		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());

		cudaMemcpy(d_temp + face_x_off, p + face_x_off, sizeof(Scalar)*grid.face_size(0), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_temp + face_y_off, p + face_y_off, sizeof(Scalar)*grid.face_size(1), cudaMemcpyDeviceToDevice);
		cudaMemcpy(d_temp + face_z_off, p + face_z_off, sizeof(Scalar)*grid.face_size(2), cudaMemcpyDeviceToDevice);

		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_fixed + face_x_off, fix, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_fixed + face_y_off, fix, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_z_off, descr->d_fixed + face_z_off, fix, grid.face_size(2));

		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());

		grid3DOperator::Cod1Mapping(grid, d_temp + edge_x_off, d_temp + edge_y_off, d_temp + edge_z_off, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off);

		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());

		grid3DOperator::D1Mapping(grid, d_temp + face_x_off, d_temp + face_y_off, d_temp + face_z_off, d_temp + edge_x_off, d_temp + edge_y_off, d_temp + edge_z_off);

		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());

		cwise_mapping_wrapper(d_temp + face_x_off, descr->d_vol + face_x_off, multi, grid.face_size(0));
		cwise_mapping_wrapper(d_temp + face_y_off, descr->d_vol + face_y_off, multi, grid.face_size(1));
		cwise_mapping_wrapper(d_temp + face_z_off, descr->d_vol + face_z_off, multi, grid.face_size(2));

		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());

		cwise_mapping_wrapper(Ap + face_x_off, d_temp + face_x_off, add, grid.face_size(0));
		cwise_mapping_wrapper(Ap + face_y_off, d_temp + face_y_off, add, grid.face_size(1));
		cwise_mapping_wrapper(Ap + face_z_off, d_temp + face_z_off, add, grid.face_size(2));

		cudaDeviceSynchronize();
		checkCudaErrors(cudaGetLastError());
	}


	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());


	// penalty
	{
		cudaMemcpy(d_temp, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

		cwise_mapping_wrapper(d_temp, descr->d_penalty, multi, dof);

		cwise_mapping_wrapper(Ap, d_temp, add, dof);
	}


	cudaDeviceSynchronize();
	checkCudaErrors(cudaGetLastError());

}

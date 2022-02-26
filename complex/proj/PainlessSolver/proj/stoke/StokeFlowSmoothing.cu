#include "StokeFlowSmoothing.h"
#include "gpuUtils.h"
#include "cuda_runtime.h"
#include "device_launch_parameters.h"

//#define GS

__forceinline__ __device__ void lap2x2inv(Scalar *y, Scalar *x, Scalar a, Scalar b, Scalar f)
{
	Scalar s = (Scalar)1 / (a*b - f*f);
	y[0] = (b*x[0] + f*x[1])*s;
	y[1] = (f*x[0] + a*x[1])*s;
}

__forceinline__ __device__ void lap2x2inv_inplace(Scalar *y, Scalar a, Scalar b, Scalar f)
{
	Scalar x[2] = { y[0],y[1] };
	Scalar s = (Scalar)1 / (a*b - f*f);
	y[0] = (b*x[0] + f*x[1])*s;
	y[1] = (f*x[0] + a*x[1])*s;
}

__forceinline__ __device__ void left_apply_diag4x4inv(Scalar *y, Scalar *x, Scalar *diag, Scalar *f)
{
	lap2x2inv(y, x, diag[0], diag[1], f[0]);
	lap2x2inv(y + 2, x + 2, diag[2], diag[3], f[1]);
}

__forceinline__ __device__ void right_apply_diag4x4inv(Scalar *y, Scalar *x, Scalar *diag, Scalar *f)
{
	lap2x2inv(y, x, diag[0], diag[1], f[0]);
	lap2x2inv(y + 2, x + 2, diag[2], diag[3], f[1]);
}

__forceinline__ __device__ void left_apply_diag4x4inv_inplace(Scalar *y, Scalar *diag, Scalar *f)
{
	lap2x2inv_inplace(y, diag[0], diag[1], f[0]);
	lap2x2inv_inplace(y + 2, diag[2], diag[3], f[1]);
}

__forceinline__ __device__ void right_apply_diag4x4inv_inplace(Scalar *y, Scalar *diag, Scalar *f)
{
	lap2x2inv_inplace(y, diag[0], diag[1], f[0]);
	lap2x2inv_inplace(y + 2, diag[2], diag[3], f[1]);
}

// local layout
// -------
//|   3   |
//|0  4  1|
//|   2   |
// -------
__forceinline__ __device__ void BlockInverse(Scalar *x, Scalar *r, Scalar *div, Scalar *penalty, bool *fixed, Scalar lap_coef)
{
	Scalar f[2], v_diag[4];
	Scalar temp[4];
	Scalar p_diag = (Scalar)0;

	f[0] = (fixed[0] || fixed[1]) ? (Scalar)0 : lap_coef;
	f[1] = (fixed[2] || fixed[3]) ? (Scalar)0 : lap_coef;

	v_diag[0] = v_diag[1] = v_diag[2] = v_diag[3] = (Scalar)4 * lap_coef;
	for (int i = 0; i < 4; i++) v_diag[i] += penalty[i];
	//for (int i = 0; i < 4; i++) v_diag[i] -= penalty[i];

	right_apply_diag4x4inv(temp, div, v_diag, f);
	for (int i = 0; i < 4; i++) x[i] = r[i];
	x[4] = r[4]; for (int i = 0; i < 4; i++) x[4] -= temp[i] * r[i];

	left_apply_diag4x4inv_inplace(x, v_diag, f);
	for (int i = 0; i < 4; i++) p_diag += temp[i] * div[i];
	p_diag = -p_diag + penalty[4];
	//p_diag = -p_diag - penalty[4];
	if (fixed[4]) x[4] = (Scalar)0; else x[4] /= p_diag;

	left_apply_diag4x4inv(temp, div, v_diag, f);
	for (int i = 0; i < 4; i++) x[i] -= x[4] * temp[i];
}

__global__ void SmoothKernel(Scalar *p, Scalar *vx, Scalar *vy, Scalar *r_p, Scalar *r_vx, Scalar *r_vy, \
	Scalar *area_vx, Scalar *area_vy, Scalar *penalty_p, Scalar *penalty_vx, Scalar *penalty_vy, \
	bool *fixed_p, bool *fixed_vx, bool *fixed_vy, grid2D grid, int flag, Scalar lap_coef = (Scalar)1)
{
	const int nbx = gridDim.x;
	const int nby = gridDim.y;

	const int bx = blockIdx.x;
	const int by = blockIdx.y;

	const int idx = threadIdx.x;
	const int idy = threadIdx.y;

	const int x = bx * 8 + idx;
	const int y = by * 8 + idy;

#ifdef GS
	if ((idx % 5 + (idy * 2) % 5 + flag) % 5)
		return;
#else
	if ((idx & 1) ^ (idy & 1) ^ flag)
		return;
#endif

	int ind[5] = { grid.face_ind(x, y, 0) ,grid.face_ind(x + 1, y, 0),grid.face_ind(x, y, 1),grid.face_ind(x, y + 1, 1),grid.cell_ind(x , y) };
	Scalar local_x[5], local_r[5];
	Scalar div[4];
	Scalar penalty[5];
	bool fixed[5];

	local_r[0] = r_vx[ind[0]];
	local_r[1] = r_vx[ind[1]];
	local_r[2] = r_vy[ind[2]];
	local_r[3] = r_vy[ind[3]];
	local_r[4] = r_p[ind[4]];

	div[0] = -area_vx[ind[0]];
	div[1] = area_vx[ind[1]];
	div[2] = -area_vy[ind[2]];
	div[3] = area_vy[ind[3]];

	fixed[0] = fixed_vx[ind[0]];
	fixed[1] = fixed_vx[ind[1]];
	fixed[2] = fixed_vy[ind[2]];
	fixed[3] = fixed_vy[ind[3]];
	fixed[4] = fixed_p[ind[4]];

	penalty[0] = penalty_vx[ind[0]];
	penalty[1] = penalty_vx[ind[1]];
	penalty[2] = penalty_vy[ind[2]];
	penalty[3] = penalty_vy[ind[3]];
	penalty[4] = penalty_p[ind[4]];

	BlockInverse(local_x, local_r, div, penalty, fixed, lap_coef);

	vx[ind[0]] = local_x[0];
	vx[ind[1]] = local_x[1];
	vy[ind[2]] = local_x[2];
	vy[ind[3]] = local_x[3];
	p[ind[4]] = local_x[4];
}

void StokeFlowSmoothing::init(StokeFlowMapping *_stoke_flow_mapping)
{
	descr = _stoke_flow_mapping->descr;
	stoke_flow_mapping = _stoke_flow_mapping;

	Nx = descr->Nx;
	Ny = descr->Ny;

	dof = stoke_flow_mapping->dof;

	cudaMalloc(&d_p, sizeof(Scalar)*dof);
	cudaMalloc(&d_Ap, sizeof(Scalar)*dof);
	cudaMalloc(&d_r, sizeof(Scalar)*dof);

	cublasCreate(&cublasHandle);
}

int StokeFlowSmoothing::xDoF()
{
	return dof;
}

int StokeFlowSmoothing::yDoF()
{
	return dof;
}

// d_r: residual
// d_p: smoothed direction
// d_Ap: mapping buffer
void StokeFlowSmoothing::applyMapping(Scalar *Ap, Scalar *p)
{
	int cell_off = 0;
	int face_x_off = cell_off + stoke_flow_mapping->grid.cell_size();
	int face_y_off = face_x_off + stoke_flow_mapping->grid.face_size(0);

	Scalar one = (Scalar)1, neg_one = (Scalar)-1, zero = (Scalar)0;

	cudaMemset(Ap, 0, sizeof(Scalar)*dof);
	cudaMemcpy(d_r, p, sizeof(Scalar)*dof, cudaMemcpyDeviceToDevice);

	for (int iter = 0; iter < iter_num; iter++)
	{
#ifdef GS
		for (int i = 0; i < 5; i++)
		{
			int flag = order ? 4 - i : i;
#else
		for (int i = 0; i < 2; i++)
		{
			int flag = order ? 1 - i : i;
#endif
			cudaMemset(d_p, 0, sizeof(Scalar)*dof);
			SmoothKernel << <dim3(Nx / 8, Ny / 8), dim3(8, 8) >> > \
				(d_p + cell_off, d_p + face_x_off, d_p + face_y_off, \
					d_r + cell_off, d_r + face_x_off, d_r + face_y_off, \
					descr->d_vol + face_x_off, descr->d_vol + face_y_off, \
					descr->d_penalty + cell_off, descr->d_penalty + face_x_off, descr->d_penalty + face_y_off, \
					descr->d_fixed + cell_off, descr->d_fixed + face_x_off, descr->d_fixed + face_y_off, \
					descr->grid, flag, descr->lap_coeff);
			//cwise_mapping_wrapper(d_p, descr->d_fixed, [=] __device__(Scalar& tv, bool tfixed) { if(tfixed) tv=(Scalar)0; }, dof);
			Scal(cublasHandle, dof, &alpha, d_p, 1);
			stoke_flow_mapping->applyMapping(d_Ap, d_p);
			Axpy(cublasHandle, dof, &one, d_p, 1, Ap, 1);
			Axpy(cublasHandle, dof, &neg_one, d_Ap, 1, d_r, 1);
		}
	}
}

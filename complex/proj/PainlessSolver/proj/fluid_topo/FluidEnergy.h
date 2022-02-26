#pragma once
#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

#include "StokeFlowDescriptor.h"
#include "FluidEnergyGrad.h"

class FluidEnergy
{
public:
	int grid_size = 128;

	MacGrid<2> mac_grid;
	Grid<2> grid;

	StokeFlowDescriptor descr;
	FluidEnergyGrad grad;

	Scalar *h_x, *d_x;
	Scalar *h_Ax, *d_Ax;
	Scalar *h_penalty, *d_penalty;

	FaceField<Scalar, 2> vel_grad;
	Field<Scalar, 2> p_grad;

	Field<Scalar, 2> cell_penalty_grad;
	FaceField<Scalar, 2> face_penalty_grad;

public:
	void init(int _grid_size);

	void init_boundary(const Field<Scalar, 2>& cell_vol, const Field<int, 2>& cell_fixed, \
		const FaceField<Scalar, 2>& face_vol, const FaceField<int, 2>& face_fixed);

	void update_penalty(const Field<Scalar, 2>& cell_penalty, const FaceField<Scalar, 2>& face_penalty);

	void compute_gradient_x(const FaceField<Scalar, 2>& vel, const Field<Scalar, 2>& p);

	void compute_gradient_q(const FaceField<Scalar, 2>& vel, const Field<Scalar, 2>& p);

	Scalar compute_energy(const FaceField<Scalar, 2>& vel, const Field<Scalar, 2>& p);
};
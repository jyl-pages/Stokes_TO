#pragma once
#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

#include "StokeFlowDescriptor3D.h"
#include "FluidEnergyGrad3D.h"

class FluidEnergy3D
{
public:
	int grid_size;

	MacGrid<3> mac_grid;
	Grid<3> grid;

	StokeFlowDescriptor3D descr;
	FluidEnergyGrad3D grad;

	Scalar *h_x, *d_x;
	Scalar *h_Ax, *d_Ax;
	Scalar *h_penalty, *d_penalty;

	FaceField<Scalar, 3> vel_grad;
	Field<Scalar, 3> p_grad;

	Field<Scalar, 3> cell_penalty_grad;
	FaceField<Scalar, 3> face_penalty_grad;

public:
	void init(int _grid_size);

	void init_boundary(const Field<Scalar, 3>& cell_vol, const Field<int, 3>& cell_fixed, \
		const FaceField<Scalar, 3>& face_vol, const FaceField<int, 3>& face_fixed);

	void update_penalty(const Field<Scalar, 3>& cell_penalty, const FaceField<Scalar, 3>& face_penalty);

	void compute_gradient_x(const FaceField<Scalar, 3>& vel, const Field<Scalar, 3>& p);

	void compute_gradient_q(const FaceField<Scalar, 3>& vel, const Field<Scalar, 3>& p);

	Scalar compute_energy(const FaceField<Scalar, 3>& vel, const Field<Scalar, 3>& p);
};
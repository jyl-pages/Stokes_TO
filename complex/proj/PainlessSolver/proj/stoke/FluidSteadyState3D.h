#pragma once

#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

#include "GMRES.h"
#include "StokeFlowDescriptor3D.h"

class FluidSteadyState3D
{
public:
	int grid_size = 32;
	int l = 4;

	MacGrid<3> mac_grid;
	Grid<3> grid;

	StokeFlowDescriptor3D *mg_descr;
	//StokeFlowDirectSolve *direct;
	bool closed;
	GMRES gmres;

	FaceField<Scalar, 3> vel;
	Field<Scalar, 3> p;

public:
	void init(int _grid_size);

	void init_boundary(const Field<Scalar, 3>& cell_vol, const Field<int, 3>& cell_fixed, \
		const FaceField<Scalar, 3>& face_vol, const FaceField<int, 3>& face_fixed, bool _closed);

	void update_b(const Field<Scalar, 3>& cell_b, const FaceField<Scalar, 3>& face_b);

	void update_penalty(const Field<Scalar, 3>& cell_penalty, const FaceField<Scalar, 3>& face_penalty);

	void solve();
};
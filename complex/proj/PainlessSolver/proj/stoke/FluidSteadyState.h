#pragma once

#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

#include "GMRES.h"
#include "StokeFlowDescriptor.h"
#include "StokeFlowSubSystem.h"
#include "StokeFlowDirectSolve.h"

class FluidSteadyState
{
public:
	int grid_size = 128;
	int l = 5;

	MacGrid<2> mac_grid;
	Grid<2> grid;

	StokeFlowDescriptor *mg_descr;
	StokeFlowDirectSolve *direct;
	bool closed;
	GMRES gmres;

	FaceField<Scalar, 2> vel;
	Field<Scalar, 2> p;

public:
	void init(int _grid_size);

	void init_boundary(const Field<Scalar,2>& cell_vol, const Field<int,2>& cell_fixed, \
		const FaceField<Scalar, 2>& face_vol, const FaceField<int, 2>& face_fixed, bool _closed);

	void update_b(const Field<Scalar, 2>& cell_b, const FaceField<Scalar, 2>& face_b);

	void update_penalty(const Field<Scalar, 2>& cell_penalty, const FaceField<Scalar, 2>& face_penalty);

	void solve();
};
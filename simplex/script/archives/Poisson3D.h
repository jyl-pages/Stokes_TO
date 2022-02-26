#pragma once
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
#include "PoissonDescriptor.h"
#include "ConjugatedGradient.h"

class Poisson3D
{
public:
	Vector3i grid_size = Vector3i(32, 32, 32);
	Grid<3> grid;
	int l = 4;

	PoissonDescriptor<3> *mg_descr;
	ConjugatedGradient cg;
	bool closed;

	Field<Scalar, 3> x, b;
	Scalar *temp_x, *temp_b;

public:
	void init(Vector3i _grid_size, Scalar _dx);

	void solve();

	void solve_fast();

	void init_boundary(const FaceField<Scalar, 3>& face_vol, const Field<int, 3>& cell_fixed, bool _closed);

	void update_b(const Field<Scalar, 3>& cell_b);
};

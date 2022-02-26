#pragma once
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
#include "PoissonDescriptor.h"
#include "ConjugatedGradient.h"

class Poisson
{
public:
	int grid_size = 128;
	Grid<2> grid;
	int l=5;

	PoissonDescriptor *mg_descr;
	ConjugatedGradient cg;
	bool closed;

	Field<Scalar, 2> x, b;
	Scalar *temp_x, *temp_b;

public:
	void init(int _grid_size, Scalar _dx);

	void solve();

	void solve_fast();

	void init_boundary(const FaceField<Scalar, 2>& face_vol, const Field<int, 2>& cell_fixed, bool _closed);

	void update_b(const Field<Scalar, 2>& cell_b);
};

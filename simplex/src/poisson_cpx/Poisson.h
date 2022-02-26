//////////////////////////////////////////////////////////////////////////
// Poisson solver interface aligning with simplex for CPX solver
// Copyright (c) (2018-), Yueyang Xianzang, Jinyuan Liu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
#include "PoissonDescriptor.h"
#include "ConjugatedGradient.h"

template<int d>
class Poisson
{
	Typedef_VectorDii(d);
public:
	VectorDi grid_size;
	Grid<d> grid;
	int l=5;

	//mg_descr describes A for equation Ax=b
	PoissonDescriptor<d> *mg_descr;
	ConjugatedGradient cg;
	bool closed;

	Field<Scalar, d> x, b;
	Scalar *temp_x, *temp_b;

public:
	Poisson() {
		if (d == 2) {
			grid_size = VectorDi::Ones() * 128;
			l = 5;
		}
		else {
			grid_size = VectorDi::Ones() * 32;
			l = 4;
		}
	}

	void init(VectorDi _grid_size, Scalar _dx);

	void solve();

	void solve_fast();

	void init_boundary(const FaceField<Scalar, d>& face_vol, const Field<int, d>& cell_fixed, bool _closed);

	void update_b(const Field<Scalar, d>& cell_b);
};

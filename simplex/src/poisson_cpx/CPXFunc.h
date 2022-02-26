//////////////////////////////////////////////////////////////////////////
// Interfaces of basic functions in CPX solver, with Field representation
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "Poisson.h"
#include "BoundaryCondition.h"

namespace CPXFunc {
	


	//make velocity a divergence-free field. alpha is alpha in FluidEuler, which means the proportion of fluid on each face. it's equal to use_alpha_for_correction=true and use_alpha_for_update_b=false in Projection.
	//template<int d> void Eliminate_Divergence(const MacGrid<d>& mac_grid, FaceField<real, d>& velocity, const FaceField<real, d>& alpha, const Field<ushort, d>& type, const BoundaryConditionMacGrid<d>& bc);
}
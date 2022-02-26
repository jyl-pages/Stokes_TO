#pragma once
#include "LinearMapping.h"
#include "form01mapping.h"
#include "grid3D.h"

namespace PoissonLike
{
	LinearMapping *GenerateJacobiPreconditioner(const grid2D& grid, LinearMapping *mapping);

	LinearMapping *GenerateDirectSolve(const grid2D& grid, LinearMapping *mapping, bool closed);

	LinearMapping *GenerateJacobiPreconditioner(const grid3D& grid, LinearMapping *mapping);

	LinearMapping *GenerateDirectSolve(const grid3D& grid, LinearMapping *mapping, bool closed);
}
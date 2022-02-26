#pragma once
#include "LinearMapping.h"
#include "form01mapping.h"
#include "grid3D.h"
#include <vector>

namespace StokeLike
{
	// Local subsystem on one cell
	// Jacobi for velocity

	LinearMapping *GenerateJacobiPreconditioner(const grid2D& grid, LinearMapping *mapping);

	LinearMapping *GenerateDirectSolve(const grid2D& grid, LinearMapping *mapping);

	std::vector<LinearMapping*> GenerateJacobiPreconditioner(const grid3D& grid, LinearMapping *mapping, int iters=5);

	LinearMapping *GenerateDirectSolve(const grid3D& grid, LinearMapping *mapping, bool closed);
}
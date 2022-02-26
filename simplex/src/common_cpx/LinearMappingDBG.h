#pragma once
#include "LinearMapping.h"

namespace LinearMappingDBG
{
	// only for small matrices
	void LinearMapping2DenseMatrix(char *path, LinearMapping *mapping);

	void LinearMapping2SparseMatrix(char *path, LinearMapping *mapping, Scalar thres);
}
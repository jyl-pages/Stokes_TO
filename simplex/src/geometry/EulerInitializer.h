//////////////////////////////////////////////////////////////////////////
// grid initializer
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __EulerInitializer_h__
#define __EulerInitializer_h__
#include "MacGrid.h"
#include "BoundaryCondition.h"

//It can assist in initializing MacGrid and BoundaryCondition of some Eularian solver structure
//For example, it can be used in FluidEulerDriver to initialize the grid with certain boundary conditions
template<int d>
class EulerInitializer {
	Typedef_VectorDii(d);
public:
	////user-specified values:
	real length;//it corresponds with 1 in prop
	int scale;
	VectorDi prop;//the main computation domain is prop*scale. Like, prop=(1,1,1) means a cube.
	Eigen::Matrix<int, 3, 2> bc_width;
	Eigen::Matrix<real, 3, 2> bc_val;

	//Generated values, must be accessed after calling Generate_Parameters()
	VectorDi cell_counts;
	real dx;
	VectorD domain_min;
	MacGrid<d> mac_grid;

	bool domain_set = false, bw_set = false, bv_set = false, all_set = false;

	//These three functions must be called before Generate_Parameters()
	//Note: boundary with 0 will still generate Neumann boundary conditions. If you want no boundary conditions, set -1
	//TODO: check if the boundary condition is legal
	void Set_Domain(real _length, int _scale, const VectorDi& _prop);
	void Set_Boundary_Width(int x0 = 1, int x1 = 1, int y0 = 1, int y1 = 1, int z0 = 1, int z1 = 1);
	void Set_Boundary_Value(real vx0 = 0.0, real vx1 = 0.0, real vy0 = 0.0, real vy1 = 0.0, real vz0 = 0.0, real vz1 = 0.0);

	//Calculate all parameters
	void Generate_Parameters(void);

	//The Following functions must be called after Generate_Parameters()
	bool In_Thick_Boundary_Grid(const VectorDi& cell, int axis, int side, int width);
	bool In_Thick_Boundary_MacGrid(int axis, const VectorDi& face, int chk_axis, int side, int width);
	//It can be called any time after Generate_Parameters(). It won't delete anything in bc. It's like an "append" operation.
	void Fill_Boundary_Condition(BoundaryConditionMacGrid<d>& bc);
};

#endif


#pragma once
#include "Driver.h"
#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

#include "FluidTopoOpt3D.h"
class FluidTopoOptDriver3D : public Driver
{
public:
	int grid_size = 32;
	int hole_size = 4;

	real frac;

	MacGrid<3> mac_grid;
	Grid<3> grid;
	Field<Scalar, 3> cell_vol, cell_b, cell_penalty;
	Field<int, 3> cell_fixed;
	FaceField<Scalar, 3> face_vol, face_b, face_penalty;
	FaceField<int, 3> face_fixed;

	FluidTopoOpt3D topo;

	//for output (as double)
	Field<double, 3> field;
	FaceField<double, 3> facefield;

	virtual void Initialize() override;

	virtual void Write_Output_Files(int frame) override;

	virtual void Run() override;

	void init_boundary();

	void init_penalty();
};
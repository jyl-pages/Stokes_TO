#pragma once
#include "Driver.h"
#include "MacGrid.h"
#include "Grid.h"
#include "Field.h"
#include "FaceField.h"

#include "FluidTopoOpt.h"
class FluidTopoOptDriver : public Driver
{
public:
	int grid_size = 128;
	int hole_size = 16;

	real frac;

	MacGrid<2> mac_grid;
	Grid<2> grid;
	Field<Scalar, 2> cell_vol, cell_b, cell_penalty;
	Field<int, 2> cell_fixed;
	FaceField<Scalar, 2> face_vol, face_b, face_penalty;
	FaceField<int, 2> face_fixed;

	FluidTopoOpt topo;

	//for output (as double)
	Field<double, 2> field;
	FaceField<double, 2> facefield;
	Field<double, 2> p_field;

	virtual void Initialize() override;

	virtual void Write_Output_Files(int frame) override;

	virtual void Run() override;

	void init_boundary();

	void init_penalty();
};
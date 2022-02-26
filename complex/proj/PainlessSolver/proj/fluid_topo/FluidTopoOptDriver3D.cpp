#include "FluidTopoOptDriver3D.h"
#include "Particles.h"
#include "Timer.h"
#include <iostream>
#include <fstream>

void FluidTopoOptDriver3D::Initialize()
{
	Vector3i sz = Vector3i(grid_size, grid_size, grid_size);
	grid.Initialize(sz, 1.);
	mac_grid.Initialize(grid);

	cell_vol.Resize(sz);
	cell_fixed.Resize(sz);
	cell_penalty.Resize(sz);
	cell_b.Resize(sz);

	face_vol.Resize(sz);
	face_fixed.Resize(sz);
	face_penalty.Resize(sz);
	face_b.Resize(sz);

	field.Resize(sz);
	facefield.Resize(sz);

	topo.Init(grid_size, (Scalar)frac, last_frame);
	init_boundary();
	init_penalty();

	topo.Init_Boundary(cell_vol, cell_fixed, cell_b, face_vol, face_fixed, face_b, cell_penalty);

	topo.Write_Output_Files_Callback = std::bind(&FluidTopoOptDriver3D::Write_Output_Files, this, std::placeholders::_1);
}

void FluidTopoOptDriver3D::Write_Output_Files(int frame)
{
	Timer<real> timer;
	timer.Reset();
	Driver::Write_Output_Files(frame);
	if (frame == 0)
	{
		std::string file_name = frame_dir + "/grid";
		mac_grid.grid.Write_To_File_3d(file_name);
		std::cout << "Write to file " << file_name << std::endl;
	}

	for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++) for (int k = 0; k < grid_size; k++)
		field(Vector3i(i, j, k)) = (double)topo.design(Vector3i(i, j, k));

	for (int d = 0; d < 3; d++)
		for (int i = 0; i < grid_size + (d == 0); i++) for (int j = 0; j < grid_size + (d == 1); j++) for (int k = 0; k < grid_size + (d == 2); k++)
			facefield(d, Vector3i(i, j, k)) = (double)topo.vel(d, Vector3i(i, j, k));

	{
		std::string file_name = frame_dir + "/velocity";
		facefield.Write_To_File_3d(file_name);
	}

	{
		std::string file_name = frame_dir + "/x";
		field.Write_To_File_3d(file_name);
	}


	// mimic boundary output
	// boundary is not in common
	{
		std::string file_name = frame_dir + "/psi_D";
		Particles<3> particles;
		iterate_cell_d(iter, grid, 3)
		{
			Vector3i cell = iter.Coord();
			if (cell_fixed(cell)) particles.X(particles.Add_Element()) = grid.Center(cell);
		}
		particles.Write_To_File_3d(file_name);
	}

	// no iterate_face_d
	// expand manually
	{
		std::string file_name = frame_dir + "/psi_N";
		Particles<3> particles;
		for (int d = 0; d < 3; d++)
			for (int i = 0; i < grid_size + (d == 0); i++) for (int j = 0; j < grid_size + (d == 1); j++) for (int k = 0; k < grid_size + (d == 2); k++)
			{
				Vector3i face = Vector3i(i, j, k);
				if (face_fixed(d, face)) particles.X(particles.Add_Element()) = mac_grid.Face_Center(d, face);
			}
		particles.Write_To_File_3d(file_name);
	}
	timer.Elapse_And_Output("IO");
}

void FluidTopoOptDriver3D::Run()
{
	topo.Optimize();
}

void FluidTopoOptDriver3D::init_boundary()
{
	cell_vol.Fill((Scalar)0); cell_fixed.Fill(0); cell_b.Fill(0);
	face_vol.Fill((Scalar)0); face_fixed.Fill(0); face_b.Fill(0);
	/// <summary>
	/// add new cases here
	/// </summary>
	if (test == 0)
	{	//64 x 64 x 64 
		for (int j = 0; j < grid_size; j++)	for (int k = 0; k < grid_size; k++)
		{
			face_fixed(0, Vector3i(0, j, k)) = 1;
			face_fixed(0, Vector3i(grid_size, j, k)) = 1;
		}

		for (int i = 0; i < grid_size; i++) for (int k = 0; k < grid_size; k++)
		{
			face_fixed(1, Vector3i(i, 0, k)) = 1;
			face_fixed(1, Vector3i(i, grid_size, k)) = 1;
		}

		for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++)
		{
			face_fixed(2, Vector3i(i, j, 0)) = 1;
			face_fixed(2, Vector3i(i, j, grid_size)) = 1;
		}

		cell_vol.Fill((Scalar)1);
		face_vol.Fill((Scalar)1);

		{
			for (int j = 28; j < 36; j++) for (int k = 28; k < 36; k++)
			{
				cell_b(Vector3i(0, j, k)) = 600;
				face_b(0, Vector3i(1, j, k)) = 600;
			}

			for (int j = 10; j < 18; j++) for (int k = 34; k < 42; k++)
			{
				cell_b(Vector3i(grid_size - 1, j, k)) = -200;
				face_b(0, Vector3i(grid_size - 1, j, k)) = 200;
			}
		}

		{
			int axis = 2;
			for (int i = 28; i < 36; i++) for (int j = 22; j < 30; j++)
			{
				cell_b(Vector3i(i, j, grid_size - 1)) = -200;
				face_b(axis, Vector3i(i, j, grid_size - 1)) = 200;
			}
			for (int i = 22; i < 30; i++) for (int j = 28; j < 36; j++)
			{
				cell_b(Vector3i(i, j, 0)) = -200;
				face_b(axis, Vector3i(i, j, 1)) = -200;
			}
		}
	}
}

void FluidTopoOptDriver3D::init_penalty()
{
	cell_penalty.Fill((Scalar)0);
	face_penalty.Fill((Scalar)0);
	if (test == 0)
	{
		topo.closed = true;
	}
}
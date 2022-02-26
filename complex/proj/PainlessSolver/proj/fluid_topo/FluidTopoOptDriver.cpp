#include "FluidTopoOptDriver.h"
#include "Particles.h"
#include "Timer.h"

void FluidTopoOptDriver::Initialize()
{
	grid.Initialize(Vector2i(grid_size, grid_size), 1.);
	mac_grid.Initialize(grid);

	cell_vol.Resize(Vector2i(grid_size, grid_size));
	cell_fixed.Resize(Vector2i(grid_size, grid_size));
	cell_penalty.Resize(Vector2i(grid_size, grid_size));
	cell_b.Resize(Vector2i(grid_size, grid_size));

	face_vol.Resize(Vector2i(grid_size, grid_size));
	face_fixed.Resize(Vector2i(grid_size, grid_size));
	face_penalty.Resize(Vector2i(grid_size, grid_size));
	face_b.Resize(Vector2i(grid_size, grid_size));

	field.Resize(Vector2i(grid_size, grid_size));
	facefield.Resize(Vector2i(grid_size, grid_size));
	p_field.Resize(Vector2i(grid_size, grid_size));

	topo.Init(grid_size, (Scalar)frac, last_frame);
	init_boundary();
	init_penalty();

	topo.Init_Boundary(cell_vol, cell_fixed, cell_b, face_vol, face_fixed, face_b, cell_penalty);

	topo.Write_Output_Files_Callback = std::bind(&FluidTopoOptDriver::Write_Output_Files, this, std::placeholders::_1);
}

void FluidTopoOptDriver::Write_Output_Files(int frame)
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

	for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++) field(Vector2i(i, j)) = (double)topo.design(Vector2i(i, j));
	for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++) p_field(Vector2i(i, j)) = (double)topo.p(Vector2i(i, j));
	for (int i = 0; i <= grid_size; i++) for (int j = 0; j < grid_size; j++) facefield(0, Vector2i(i, j)) = (double)topo.vel(0, Vector2i(i, j));
	for (int i = 0; i < grid_size; i++) for (int j = 0; j <= grid_size; j++) facefield(1, Vector2i(i, j)) = (double)topo.vel(1, Vector2i(i, j));

	{
		std::string file_name = frame_dir + "/velocity";
		facefield.Write_To_File_3d(file_name);
	}

	{
		Field<double, 2> field_wall; field_wall.Resize(Vector2i(grid_size, grid_size));
		for (int i = 0; i < grid_size; i++) for (int j = 0; j < grid_size; j++)
		{
			field_wall(Vector2i(i, j)) = field(Vector2i(i, j));
		}
		std::string file_name = frame_dir + "/x";
		field_wall.Write_To_File_3d(file_name);
	}


	// mimic boundary output
	// boundary is not in common
	{
		std::string file_name = frame_dir + "/psi_D";
		Particles<2> particles;
		iterate_cell_d(iter, grid, 2)
		{
			Vector2i cell = iter.Coord();
			if (cell_fixed(cell)) particles.X(particles.Add_Element()) = grid.Center(cell);
		}
		particles.Write_To_File_3d(file_name);
	}

	// no iterate_face_d
	// expand manually
	{
		std::string file_name = frame_dir + "/psi_N";
		Particles<2> particles;
		for (int i = 0; i <= grid_size; i++) for (int j = 0; j < grid_size; j++)
		{
			Vector2i face = Vector2i(i, j);
			if (face_fixed(0, face)) particles.X(particles.Add_Element()) = mac_grid.Face_Center(0, face);
		}
		for (int j = 0; j <= grid_size; j++) for (int i = 0; i < grid_size; i++)
		{
			Vector2i face = Vector2i(i, j);
			if (face_fixed(1, face)) particles.X(particles.Add_Element()) = mac_grid.Face_Center(1, face);
		}
		particles.Write_To_File_3d(file_name);
	}
	timer.Elapse_And_Output("IO");
}

void FluidTopoOptDriver::Run()
{
	topo.Optimize();
}

void FluidTopoOptDriver::init_boundary()
{
	cell_vol.Fill((Scalar)0); cell_fixed.Fill(0); cell_b.Fill(0);
	face_vol.Fill((Scalar)0); face_fixed.Fill(0); face_b.Fill(0);
	/// <summary>
	/// add new cases here
	/// </summary>
	if (test == 0)
	{// ./fluid_topo.exe -driver 0 -test 0 -s 256 -o output -frac 0.15 -lf 40
		for (int i = 0; i < grid_size; i++)	face_fixed(0, Vector2i(0, i)) = 1;
		for (int i = 0; i < grid_size; i++) face_fixed(0, Vector2i(grid_size, i)) = 1;

		for (int i = 0; i < grid_size; i++) face_fixed(1, Vector2i(i, 0)) = 1;
		for (int i = 0; i < grid_size; i++) face_fixed(1, Vector2i(i, grid_size)) = 1;

		cell_vol.Fill((Scalar)1);
		face_vol.Fill((Scalar)1);

		int O = 16;
		
		int in_half_width = grid_size / 32;
		int out_half_width = grid_size / (4 * O);

		Scalar e = 60;
		Scalar r = e * (2 * in_half_width + 1) / (2 * out_half_width + 1) / O;
		//e -= 0.4 * r * (2 * out_half_width + 1) / (2 * in_half_width + 1);
		{
			for (int i = 0; i < grid_size; i++)
			{
				if (abs(i - grid_size / 2) > in_half_width) continue;
				cell_b(Vector2i(0, i)) = e;
				face_b(0, Vector2i(1, i)) = e;
			}
		}

		for (int o = 0; o < O; o++)
		{
			for (int i = 0; i < grid_size; i++)
			{
				if (abs(i - grid_size / (2 * O)*(2 * o + 1)) > out_half_width) continue;
				if (o == 0 || o == O - 1) {
					cell_b(Vector2i(grid_size - 1, i)) = -r;
					face_b(0, Vector2i(grid_size - 1, i)) = r;
				}
				else {
					cell_b(Vector2i(grid_size - 1, i)) = -r;
					face_b(0, Vector2i(grid_size - 1, i)) = r;
				}
			}
		}

	}
}

void FluidTopoOptDriver::init_penalty()
{
	cell_penalty.Fill((Scalar)0);
	face_penalty.Fill((Scalar)0);
	if (test == 0)
	{
		topo.closed = true;
	}

}

#include <iostream>
#include "EulerInitializer.h"
#include "AuxFunc.h"


template<int d>
void EulerInitializer<d>::Set_Domain(real _length, int _scale, const VectorDi& _prop)
{
	length = _length;
	scale = _scale;
	prop = _prop;
	domain_set = true;
}

template<int d>
void EulerInitializer<d>::Set_Boundary_Width(int x0, int x1, int y0, int y1, int z0, int z1)
{
	bc_width << x0, x1, y0, y1, z0, z1;
	bw_set = true;
}

template<int d>
void EulerInitializer<d>::Set_Boundary_Value(real vx0, real vx1, real vy0, real vy1, real vz0, real vz1)
{
	bc_val << vx0, vx1, vy0, vy1, vz0, vz1;
	bv_set = true;
}

template<int d>
void EulerInitializer<d>::Generate_Parameters(void)
{
	if (domain_set && bw_set && bv_set) {
		Vector3i add3 = bc_width.rowwise().sum();
		VectorDi add_width = AuxFunc::Vi<d>(add3[0], add3[1], add3[2]);
		cell_counts = prop * scale + add_width;
		dx = length / scale;
		VectorD domain_len = cell_counts.template cast<real>() * dx;
		domain_min = domain_len * (-0.5);
		mac_grid.Initialize(cell_counts, dx, domain_min);
		all_set = true;
	}
	else {
		std::cerr << "Error: [EulerInitializer] Generate_Parameters(): parameters not fully set\n";
	}
}

template<int d>
bool EulerInitializer<d>::In_Thick_Boundary_Grid(const VectorDi& cell, int axis, int side, int width)
{
	if (!all_set) std::cerr << "Error: [EulerInitializer] all_set==false\n";
	if (side == 0) return cell[axis] < width;
	else if (side == 1) return cell[axis] >= cell_counts[axis] - width;
	return false;
}

template<int d>
bool EulerInitializer<d>::In_Thick_Boundary_MacGrid(int axis, const VectorDi& face, int chk_axis, int side, int width)
{
	if (!all_set) std::cerr << "Error: [EulerInitializer] all_set==false\n";
	if (side == 0) {
		if (axis == chk_axis) return face[chk_axis] <= width;
		else return face[chk_axis] < width;
	}
	else if (side == 1) {
		return face[chk_axis] >= cell_counts[chk_axis] - width;
	}
	else return false;
}

template<int d>
void EulerInitializer<d>::Fill_Boundary_Condition(BoundaryConditionMacGrid<d>& bc)
{
	if (!all_set) std::cerr << "Error: [EulerInitializer] all_set==false\n";
	iterate_cell(iter, mac_grid.grid) {
		const VectorDi& cell = iter.Coord();
		int in_cnt = 0;
		for (int axis = 0; axis < d; axis++) {
			for (int side = 0; side < 2; side++) {
				if (In_Thick_Boundary_Grid(cell, axis, side, bc_width(axis, side))) in_cnt++;
			}
		}
		if (in_cnt > 0) bc.Set_Psi_D(cell, (ushort)CellType::Solid);
	}
	iterate_face(axis, iter, mac_grid) {
		const VectorDi& face = iter.Coord();
		int in_cnt = 0, _chk_axis, _side;
		for (int chk_axis = 0; chk_axis < d; chk_axis++) {
			for (int side = 0; side < 2; side++) {
				if (In_Thick_Boundary_MacGrid(axis, face, chk_axis, side, bc_width(chk_axis, side))) {
					in_cnt++;
					_chk_axis = axis;
					_side = side;
				}
			}
		}
		if (in_cnt > 0) {
			if (in_cnt == 1 && axis == _chk_axis) {
				bc.Set_Psi_N(axis, face, bc_val(_chk_axis, _side));
			}
			else bc.Set_Psi_N(axis, face, 0);
		}
	}
}



template class EulerInitializer<2>;
template class EulerInitializer<3>;
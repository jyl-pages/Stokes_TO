#include "ProjectionMeta.h"

template<int d>
ProjectionMeta<d>::ProjectionMeta(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, FaceField<real, d>* _alpha, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc, const SolverType& _mode)
	:Base(_mac_grid,_velocity,_type,_bc,_mode)
{ Initialize(_alpha); }

template<int d>
ProjectionMeta<d>::ProjectionMeta(MacGrid<d>& _mac_grid, FaceField<real, d>& _velocity, const SolverType& _mode)
	:Base(&_mac_grid,&_velocity,nullptr,nullptr,_mode)
{ Initialize(nullptr); }

template<int d>
ProjectionMeta<d>::ProjectionMeta(MacGrid<d>& _mac_grid, FaceField<real, d>& _velocity, FaceField<real, d>& _alpha, Field<ushort, d>& _type, BoundaryConditionMacGrid<d>& _bc, const SolverType& _mode)
	:Base(&_mac_grid,&_velocity,&_type,&_bc,_mode)
{ Initialize(&_alpha); }

template<int d>
ProjectionMeta<d>::~ProjectionMeta()
{
	if (alpha != nullptr && own_alpha) delete alpha;
}

template<int d>
inline void ProjectionMeta<d>::Initialize(FaceField<real, d>* _alpha)
{
	const VectorDi& counts = mac_grid->grid.cell_counts;
	if (alpha != nullptr && own_alpha)delete alpha;
	if (_alpha == nullptr) { alpha = new FaceField<real, d>(counts, (real)1); own_alpha = true; }
	else { alpha = _alpha; own_alpha = false; }
}

template<int d>
void ProjectionMeta<d>::Update_b()
{
	int b_size = (int)matrix_to_grid.size();
#pragma omp parallel for
	for (auto r = 0; r < b_size; r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		real div = (real)0;
		for (int axis = 0; axis < d; axis++) {
			div += use_alpha_for_update_b ?
				((*alpha)(axis, cell + VectorDi::Unit(axis)) * (*velocity)(axis, cell + VectorDi::Unit(axis)) - (*alpha)(axis, cell) * (*velocity)(axis, cell))
				: ((*velocity)(axis, cell + VectorDi::Unit(axis)) - (*velocity)(axis, cell));
		}
		div_u[r] = -div;
	}	////Attention: use negative div here. We are solving -lap p=-div u
}

template<int d>
real ProjectionMeta<d>::Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx) const
{
	VectorDi nb_cell = Grid<d>::Nb_C(fluid_cell, nbidx);
	int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(fluid_cell, nbidx, axis, face);
	if (bc->Is_Psi_N(axis, face))return 0;
	if (Is_Fluid_Cell(nb_cell)) return -(*alpha)(axis, face);
	return 0;
}

template<int d>
real ProjectionMeta<d>::Diag_Face_Term(const int& axis, const VectorDi& face) const
{
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Fluid_Cell(cell[0])) std::swap(cell[0], cell[1]);
	if (!Is_Fluid_Cell(cell[0])) return (real)0;
	if (bc->Is_Psi_N(axis, face)) return (real)0;
	return (*alpha)(axis, face);
}

template<int d>
real ProjectionMeta<d>::Velocity_Offset(const int& axis, const VectorDi& face) const
{
	if (bc->Is_Psi_N(axis, face)) return 0;
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	real cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = (Is_Valid_Cell(cell[i]) && Is_Fluid_Cell(cell[i])) ? p[grid_to_matrix(cell[i])] : (real)0;
	return -(cell_p[1] - cell_p[0]) * (use_alpha_for_correction ? (*alpha)(axis, face) : (real)1);
}

//////////////////////////////////////////////////////////////////////////
//Deprecated, planning to remove in next versions
//Only ibm_rigid_coupling use them (Mengdi)
////TODO: merge these customized div/cell/vel input by using the referenced variables

template<int d> void ProjectionMeta<d>::Allocate_System(std::function<bool(const int)>& fluid_cell)
{
	if (is_A_initialized && !update_A)return;
	Build_Grid_Cell_Matrix_Bijective_Mapping(mac_grid->grid, fluid_cell, grid_to_matrix, matrix_to_grid);

	////Setup A and div_u
	int n = (int)matrix_to_grid.size();
	A.resize(n, n); p.resize(n); p.fill((real)0); div_u.resize(n); div_u.fill((real)0);
}


template<int d> void ProjectionMeta<d>::Build(std::function<bool(const int)>& fluid_cell)
{
	Timer<real> timer; timer.Reset();
	Allocate_System(fluid_cell);
	if (verbose)timer.Elapse_And_Output_And_Reset("Allocate A");
	Update_A(fluid_cell);
	if (verbose)timer.Elapse_And_Output_And_Reset("Assemble A");
	Update_b();
	if (verbose)timer.Elapse_And_Output_And_Reset("Assemble b");
}

template<int d> void ProjectionMeta<d>::Update_A(std::function<bool(const int)>& fluid_cell)
{
	if (is_A_initialized && !update_A)return;
	if (solver_mode == SolverType::CPX_GPU) {
		Prepare_CPX_System();
	}
	else {
		Timer<real> timer; timer.Reset();
		Array<TripletT> elements;
		for (auto r = 0; r < matrix_to_grid.size(); r++) {
			const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
			////off-diagonal elements
			for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
				VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
				int nb_cell_idx = mac_grid->grid.Cell_Index(nb_cell);
				if (Is_Valid_Cell(nb_cell) && fluid_cell(nb_cell_idx)) {
					int c = grid_to_matrix(nb_cell);
					int axis = 0; int side = 0; mac_grid->grid.Nb_C_Axis_And_Side(i, axis, side);
					VectorDi face = cell + VectorDi::Unit(axis) * side;
					real a = (*alpha)(axis, face);
					elements.push_back(TripletT((int)r, (int)c, (real)-a));
				}
			}
			////diagonal elements
			real dia_coef = (real)0;
			for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
				real coef = Dia_Coef(cell, i); dia_coef += coef;
			}
			elements.push_back(TripletT((int)r, (int)r, dia_coef));
		}
		if (verbose)timer.Elapse_And_Output_And_Reset("Update A elements");
		A.setFromTriplets(elements.begin(), elements.end());
		A.makeCompressed();
		if (verbose)timer.Elapse_And_Output_And_Reset("Assemble A to sp_mtx");
	}
	is_A_initialized = true;
}

template<int d> void ProjectionMeta<d>::Correction(std::function<bool(const int)>& fluid_cell, FaceField<real, d>* correct_vel /*= nullptr*/)
{
	if (correct_vel == nullptr) correct_vel = velocity;
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			if ((*bc).Is_Psi_N(axis, face))continue;
			VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
			real cell_p[2];
			for (int i = 0; i < 2; i++) {
				int cell_idx = mac_grid->grid.Cell_Index(cell[i]);
				if (Is_Valid_Cell(cell[i]) && fluid_cell(cell_idx)) cell_p[i] = p[grid_to_matrix(cell[i])];
				else cell_p[i] = 0;
			}
			(*correct_vel)(axis, face) -= (use_alpha_for_correction ? (*alpha)(axis, face) : (real)1) * (cell_p[1] - cell_p[0]);
		}
	}
}


template class ProjectionMeta<2>;
template class ProjectionMeta<3>;
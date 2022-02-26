//////////////////////////////////////////////////////////////////////////
// Project a vector field with free boundary to divergence free on a MAC grid
// Copyright (c) (2018-), Fan Feng, Shuqi Yang, Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ProjectionIrregular.h"

//////////////////////////////////////////////////////////////////////////
////constructors
template<int d> ProjectionIrregular<d>::ProjectionIrregular(MacGrid<d>& _mac_grid, FaceField<real, d>& _velocity, LevelSet<d>& _levelset, BoundaryConditionMacGrid<d>& _bc, const SolverType& _mode)
	:Base(&_mac_grid, &_velocity, nullptr, &_bc, _mode)
{
	verbose = false;
	Initialize(&_levelset);
}

template<int d> void ProjectionIrregular<d>::Initialize(LevelSet<d>* _levelset)
{
	if(_levelset!=nullptr&&own_levelset)delete levelset;
	if(_levelset==nullptr){levelset=new LevelSet<d>(mac_grid->grid);own_levelset=true;}
	else{levelset=_levelset;own_levelset=false;}
}

//////////////////////////////////////////////////////////////////////////
////projection functions
template<int d> void ProjectionIrregular<d>::Update_A()
{
	if (use_update_A_in_parallel)Update_A_In_Parallel();
	else Base::Update_A();
}

////parallelized assembling A on CPU
template<int d> void ProjectionIrregular<d>::Update_A_In_Parallel()
{
	if(is_A_initialized&&!update_A)return;

	if(Is_Fluid_Interior_Cell_Index==nullptr){std::cerr<<"[Error] ProjectionIrregular: Is_Fluid_Interior_Cell_Index is null"<<std::endl;return;}

	if (solver_mode == SolverType::CPX_GPU) {
		Prepare_CPX_System();
	}
	else {

		constexpr int int_nb_n = Grid<d>::Number_Of_Nb_C();
		int mtx_size = (int)matrix_to_grid.size();

		int* ptr = A.outerIndexPtr(); ptr[0] = 0;
#pragma omp parallel for
		for (auto r = 0; r < mtx_size; r++) {
			const int idx = matrix_to_grid[r];
			int nb_n = Fluid_Neighbor_Number_Index(idx);
			ptr[r + 1] = nb_n + 1;
		}

		int nnz = 0; for (int i = 0; i <= A.rows(); i++) { nnz += ptr[i]; ptr[i] = nnz; }

		A.resizeNonZeros(nnz);
		int* col = A.innerIndexPtr();
		real* val = A.valuePtr();
#pragma omp parallel for
		for (auto r = 0; r < mtx_size; r++) {
			const int idx = matrix_to_grid[r];
			if (Is_Fluid_Interior_Cell_Index(idx)) {
				////diagonal
				col[ptr[r]] = r; val[ptr[r]] = (real)int_nb_n;
				////off-diagonal
				for (int i = 0; i < int_nb_n; i++) {
					int nb_idx = mac_grid->grid.Cell_Nb_C_Index(idx, i);
					int c = grid_to_matrix.array[nb_idx];
					col[ptr[r] + i + 1] = c; val[ptr[r] + i + 1] = (real)-1;
				}
			}
			else {
				const VectorDi& cell = mac_grid->grid.Cell_Coord(idx);
				//real fluid_phi = levelset->phi(cell);
				////diagonal
				real dia_coef = (real)0;
				for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
					int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
					dia_coef += Diag_Face_Term(axis, face);
					//dia_coef += Diag_Face_Term(cell, i, fluid_phi);
				}
				col[ptr[r]] = r; val[ptr[r]] = dia_coef;
				////off-diagonal
				int nb_idx = 0;
				for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
					real term = Off_Diag_Term(cell, i);
					VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
					if (Is_Fluid_Cell(nb_cell)) {
						int c = grid_to_matrix(nb_cell);
						col[ptr[r] + nb_idx + 1] = c; val[ptr[r] + nb_idx + 1] = term;
						nb_idx++;
					}
				}
			}
		}
	}
	is_A_initialized=true;
}

template<int d> void ProjectionIrregular<d>::Apply_Jump_Condition_To_b()
{
	if(!use_levelset_interface)return;
	real one_over_dx=(real)1/mac_grid->grid.dx;
	int b_size=(int)matrix_to_grid.size();
	#pragma omp parallel for
	for (int r = 0; r < b_size; r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		if (Is_Fluid_Cell(cell)) {
			for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
				real coef = -Off_Diag_Term(cell, i);
				real theta; bool is_intf; real intf_coef = Intf_Coef(cell, i, theta, is_intf);
				if (is_intf) {
					VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
					VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell) + theta * mac_grid->grid.Center(nb_cell);
					real p_jump = Jump_Condition(intf_pos);
					div_u[r] += coef * intf_coef * p_jump * one_over_dx;
				}
			}
		}
	}
}

////need to specify target_vol and current_vol 
template<int d> void ProjectionIrregular<d>::Apply_Vol_Control_To_b()
{
	if(!use_vol_control)return;
	if(target_vol==(real)-1){std::cerr<<"[Error] ProjectionIrregular: target_vol not set"<<std::endl;return;}

	if(calc_current_vol)current_vol=levelset->Total_Volume();

	real vol_correction=(target_vol-current_vol)/current_vol;
	real cell_vol=pow(mac_grid->grid.dx,d);

	if(verbose)std::cout<<"vol correction: "<<vol_correction<<std::endl;

	int b_size=(int)matrix_to_grid.size();
	#pragma omp parallel for
	for (int r = 0; r < b_size; r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		if (Is_Fluid_Cell(cell)) {
			real vol = levelset->Cell_Fraction(cell) * cell_vol;
			real cell_div = vol_correction;
			div_u[r] += vol_control_ks * mac_grid->grid.dx * cell_div;
		}
	}
}

template<int d> void ProjectionIrregular<d>::Apply_Implicit_Surface_tension(const real dt)
{
	if (Is_Interface_Face_Index == nullptr || Dirac == nullptr) { std::cerr << "[Error] ProjectionIrregular: Is_Interface_Face_Index or Dirac is null" << std::endl; return; }
	macgrid_to_matrix.Resize(mac_grid->grid.cell_counts,-1);
	matrix_to_macgrid.clear();
	Build_MacGrid_Face_Matrix_Bijective_Mapping(*mac_grid, Is_Interface_Face_Index, macgrid_to_matrix, matrix_to_macgrid);
	int n = matrix_to_macgrid.size();

	SparseMatrixT B;
	VectorX u_new;
	VectorX u_old;

	// setup A, x, and b
	B.resize(n, n);
	u_new.resize(n); u_new.fill((real)0);
	u_old.resize(n); u_old.fill((real)0);
	Array<TripletT> elements;

	for (int r = 0; r < n; r++) {
		const int axis = matrix_to_macgrid[r].first;
		const VectorDi& face = mac_grid->face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
		VectorD pos = mac_grid->Face_Center(axis, face);
		real kappa = (*levelset).Curvature(pos);
		VectorD normal = (*levelset).Normal(pos);
		real dia_coef = (real)1;
		u_old[r] = (*velocity)(axis, face) - sigma * Dirac((*levelset).Phi(pos)) * dt * kappa * normal[axis];	//neglect Jacobi and Hessian
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			VectorDi nb_face = Grid<d>::Nb_C(face, i);
			VectorD nb_pos = mac_grid->Face_Center(axis, nb_face);
			if (!mac_grid->Valid_Face(axis, nb_face) || (*bc).Is_Psi_N(axis, nb_face)) continue;
			real a = sigma * Dirac((*levelset).Phi((pos + nb_pos) * (real).5)) * dt * dt / (mac_grid->grid.dx * mac_grid->grid.dx);
			dia_coef += a;
			int c = macgrid_to_matrix(axis, nb_face);
			if (Is_Interface_Face_Index(std::make_pair(axis, mac_grid->face_grids[axis].Node_Index(nb_face)))) {
				elements.push_back(TripletT(r, c, -a));
			}
			else u_old[r] += a * (*velocity)(axis, nb_face);
		}
		elements.push_back(TripletT(r, r, dia_coef));
	}
	B.setFromTriplets(elements.begin(), elements.end()); B.makeCompressed();

#ifdef USE_CUDA
	MultiGridCuda::Preconditioned_Conjugate_Gradient<SparseMatrix<real>, real>(&B, &u_old[0], &u_new[0], 3000, (real)1e-5,/*diagonal precond*/true);	////GPU D-PCG
#else
	KrylovSolver::Params params; KrylovSolver::ICPCG(B, u_new, u_old, params);	////CPU IC-PCG
#endif

#pragma omp parallel for
	for (int r = 0; r < n; r++) {
		const int axis = matrix_to_macgrid[r].first;
		const VectorDi& face = mac_grid->face_grids[axis].Node_Coord(matrix_to_macgrid[r].second);
		(*velocity)(axis, face) = u_new[r];
	}
}

template<int d> void ProjectionIrregular<d>::Build()
{
	Timer<real> timer;timer.Reset();
	Allocate_System();
	if(verbose)timer.Elapse_And_Output_And_Reset("Allocate A");
	Update_A();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble A");
	Update_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble b");
	Apply_Jump_Condition_To_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Apply jump condition to b");
	Apply_Vol_Control_To_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Apply volumn control to b");
}

template<int d> void ProjectionIrregular<d>::Project()
{
	Timer<real> timer;		if(verbose)timer.Reset();
	if(use_surface_tension&&!use_jump_condition){
		Apply_Implicit_Surface_tension(current_dt);
							if(verbose)timer.Elapse_And_Output_And_Reset("implicit surface tension");}
	Build();				if(verbose)timer.Elapse_And_Output_And_Reset("build");
	Solve();				if(verbose)timer.Elapse_And_Output_And_Reset("solve");
	Correction();			if(verbose)timer.Elapse_And_Output_And_Reset("correction");
}




template<int d> real ProjectionIrregular<d>::Intf_Coef(const VectorDi& cell,const int i,real& theta,bool& is_intf) const
{
	VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
	if (mac_grid->grid.Valid_Cell(nb_cell)) {
		real phi0 = (*levelset).phi(cell); real phi1 = (*levelset).phi(nb_cell);
		return Intf_Coef(phi0, phi1, theta, is_intf);
	}
	return (real)1;
}

template<int d> real ProjectionIrregular<d>::Intf_Coef(const real& phi0,const real& phi1,real& theta,bool& is_intf) const
{
	real max_scale=(real)1e3;is_intf=false;
	if (LevelSet<d>::Interface(phi0, phi1)) {
		is_intf = true;
		theta = LevelSet<d>::Theta(phi0, phi1);
		real scale = std::min((real)1 / (real)theta, max_scale); return scale;
	}
	return (real)1;
}

template<int d>
int ProjectionIrregular<d>::Fluid_Neighbor_Number_Index(const int& fluid_cell_idx) const
{
	int nb_n = 0;
	if (Is_Fluid_Interior_Cell_Index(fluid_cell_idx))nb_n = Grid<d>::Number_Of_Nb_C();
	else {	////boundary
		const VectorDi& cell = mac_grid->grid.Cell_Coord(fluid_cell_idx);
		for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
			VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
			if (Is_Fluid_Cell(nb_cell)) { nb_n++; }
		}
	}
	return nb_n;
}

//////////////////////////////////////////////////////////////////////////
////Physical interface functions
//see: https://wmdcstdio.com/2021/07/11/projection-matrix-terms/
template<int d>
real ProjectionIrregular<d>::Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx) const
{
	VectorDi nb_cell = Grid<d>::Nb_C(fluid_cell, nbidx);
	int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(fluid_cell, nbidx, axis, face);
	if (bc->Is_Psi_N(axis, face))return 0;
	if (Is_Fluid_Cell(nb_cell)) return (real)-1;
	return 0;
}

template<int d>
real ProjectionIrregular<d>::Diag_Face_Term(const int& axis, const VectorDi& face) const
{
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Fluid_Cell(cell[0])) std::swap(cell[0], cell[1]);
	if (!Is_Fluid_Cell(cell[0])) return (real)0;
	real fluid_phi = levelset->phi(cell[0]);
	if (bc->Is_Psi_N(axis, face)) return (real)0;
	if (use_levelset_interface && mac_grid->grid.Valid_Cell(cell[1])) {
		real nb_phi = levelset->phi(cell[1]);
		real theta; bool is_intf; real intf_coef = Intf_Coef(fluid_phi, nb_phi, theta, is_intf);
		if (is_intf) return intf_coef;
	}
	return (real)1;
}

template<int d>
real ProjectionIrregular<d>::Velocity_Offset(const int& axis, const VectorDi& face) const
{
	if ((*bc).Is_Psi_N(axis, face))return 0;
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	real cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = Is_Fluid_Cell(cell[i]) ? p[grid_to_matrix(cell[i])] : (real)0;	////air cell has zero pressure

	if (use_levelset_interface && mac_grid->grid.Valid_Cell(cell[0]) && mac_grid->grid.Valid_Cell(cell[1])) {
		real one_over_dx = (real)1 / mac_grid->grid.dx;
		real phi[2]; for (int i = 0; i < 2; i++)phi[i] = (*levelset).phi(cell[i]);
		int f_idx = (Is_Fluid_Cell(cell[0]) ? 0 : 1); int a_idx = (f_idx == 0 ? 1 : 0);
		real theta; bool is_intf; real intf_coef = Intf_Coef(phi[f_idx], phi[a_idx], theta, is_intf);
		if (is_intf) {		////set the boundary cell pressure to be the pressure jump value multiplying intf coef
			cell_p[f_idx] *= intf_coef;
			if (use_surface_tension && use_jump_condition) {
				VectorD intf_pos = ((real)1 - theta) * mac_grid->grid.Center(cell[f_idx]) + theta * mac_grid->grid.Center(cell[a_idx]);
				real p_jump = Jump_Condition(intf_pos);
				cell_p[a_idx] = p_jump * one_over_dx * intf_coef;
			}
		}
	}

	return -(cell_p[1] - cell_p[0]);
}


template class ProjectionIrregular<2>;
template class ProjectionIrregular<3>;

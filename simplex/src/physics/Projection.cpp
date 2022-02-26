//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free
// Copyright (c) (2018-),Bo Zhu,Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Projection.h"

//////////////////////////////////////////////////////////////////////////
////Constructor
template<int d>
Projection<d>::Projection(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc, const SolverType& _mode)
{
	Initialize(_mac_grid, _velocity, _type, _bc, _mode);
}


template<int d> Projection<d>::~Projection()
{
	if (mac_grid != nullptr && own_grid) delete mac_grid;
	if (velocity != nullptr && own_velocity) delete velocity;
	if (type != nullptr && own_type) delete type;
	if (bc != nullptr && own_bc) delete bc;
}

template<int d> void Projection<d>::Initialize(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc, const SolverType& _mode)
{
	solver_mode = _mode;
	if (solver_mode == SolverType::AUTO) {
		Auto_Select_Mode();
	}

	if (mac_grid != nullptr && own_grid) delete mac_grid;
	if (_mac_grid == nullptr) { mac_grid = new MacGrid<d>(); own_grid = true; }
	else { mac_grid = _mac_grid; own_grid = false; }

	const VectorDi& counts = mac_grid->grid.cell_counts;

	if (velocity != nullptr && own_velocity)delete velocity;
	if (_velocity == nullptr) { velocity = new FaceField<real, d>(counts, (real)0); own_velocity = true; }
	else { velocity = _velocity; own_velocity = false; }

	if (type != nullptr && own_type)delete type;
	if (_type == nullptr) { type = new Field<ushort, d>(counts, (ushort)CellType::Fluid); own_type = true; }
	else { type = _type; own_type = false; }

	if (bc != nullptr && own_bc)delete bc;
	if (_bc == nullptr) { bc = new BoundaryConditionMacGrid<d>(*_mac_grid); own_bc = true; }
	else { bc = _bc; own_bc = false; }
}

template<int d>
void Projection<d>::Auto_Select_Mode(void)
{
	//cpx(cpu)->multigrid(gpu)->krylov(cpu)
	solver_mode = SolverType::KRYLOV_CPU;
#ifdef USE_CUDA
	solver_mode = SolverType::MULTIGRID_AUTO;
#endif
//#ifdef USE_CPX
//	solver_mode = SolverType::CPX_GPU;
//#endif
}


//////////////////////////////////////////////////////////////////////////
////Build linear system
template<int d> void Projection<d>::Allocate_System()
{
	if(is_A_initialized&&!update_A)return;
	std::function<bool(const int)> fluid_cell=[=](const int idx)->bool{return this->Is_Fluid_Cell(mac_grid->grid.Cell_Coord(idx));};
	Build_Grid_Cell_Matrix_Bijective_Mapping(mac_grid->grid,fluid_cell,grid_to_matrix,matrix_to_grid);

	////Setup A and div_u
	int n=(int)matrix_to_grid.size();
	A.resize(n,n);p.resize(n);p.fill((real)0);div_u.resize(n);div_u.fill((real)0);
}

template<int d>
void Projection<d>::Prepare_CPX_System(void)
{
#ifdef USE_CPX
	int n = mac_grid->grid.cell_counts[0];
	//for (int axis = 0; axis < d; axis++) { if (mac_grid->grid.cell_counts[axis] != n) { AuxFunc::Crash_With_Info("Projection::Prepare_CPX_System: cell size not equal"); } }
	if (n % 8) AuxFunc::Crash_With_Info("Projection::Prepare_CPX_System: cell size must can be divided by 8");

	if (!cpx_inited) {
		cpx_poisson.init(mac_grid->grid.cell_counts, mac_grid->grid.dx);
		cpx_poisson.cg.verbose = false;
		cpx_inited = true;
	}

	//using CPX_Descr = typename If<d == 2, PoissonDescriptor, PoissonDescriptor3D>::Type;

	PoissonDescriptor<d>& descr = cpx_poisson.mg_descr[cpx_poisson.l - 1];
	bool* fixed = new bool[descr.size];
	Scalar* vol = new Scalar[descr.fsize];
	memset(fixed, 0, sizeof(bool) * descr.size);
	memset(vol, 0, sizeof(Scalar) * descr.fsize);

	int offset = 0;
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->Number_Of_Faces(axis);
#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->Face_Coord(axis, i);
			int face_ind = -1;
			if constexpr (d == 2) { face_ind = offset + descr.grid.face_ind(face[0], face[1], axis); }
			else if constexpr (d == 3) { face_ind = offset + descr.grid.face_ind(face[0], face[1], face[2], axis); }
			else AuxFunc::Crash_With_Info("Projection<d>::Prepare_CPX_System crashed");
			vol[face_ind] = Diag_Face_Term(axis, face);
		}
		offset += descr.grid.face_size(axis);
	}

	int cell_num = mac_grid->grid.cell_counts.prod();
#pragma omp parallel for
	for (int i = 0; i < cell_num; i++) {
		VectorDi cell = mac_grid->grid.Cell_Coord(i);
		int cell_ind = -1;
		if constexpr (d == 2) cell_ind = descr.grid.cell_ind(cell[0], cell[1]);
		else if constexpr (d == 3) cell_ind = descr.grid.cell_ind(cell[0], cell[1], cell[2]);
		else AuxFunc::Crash_With_Info("Projection<d>::Prepare_CPX_System crashed");
		fixed[cell_ind] = !Is_Fluid_Cell(cell);
	}

	descr.setFixed(fixed);
	descr.setVol(vol);
	descr.toDevice();
	descr.finish();

	cpx_poisson.closed = false;//NOTE: set to not closed here
	delete[] fixed;
	delete[] vol;
#else
	AuxFunc::Crash_With_Info("Projection<d>::Prepare_CPX_System must compile CPX");
#endif USE_CPX
}

template<int d> void Projection<d>::Update_A()
{
	if(is_A_initialized&&!update_A)return;
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
				real term = Off_Diag_Term(cell, i);
				VectorDi nb_cell = Grid<d>::Nb_C(cell, i);
				if (Is_Fluid_Cell(nb_cell)) {
					int c = grid_to_matrix(nb_cell);
					elements.push_back(TripletT((int)r, (int)c, term));
				}
			}
			////diagonal elements
			real dia_coef = (real)0;
			for (int i = 0; i < Grid<d>::Number_Of_Nb_C(); i++) {
				int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
				dia_coef += Diag_Face_Term(axis, face);
			}
			elements.push_back(TripletT((int)r, (int)r, dia_coef));
		}

		if (verbose)timer.Elapse_And_Output_And_Reset("Update A elements");
		A.setFromTriplets(elements.begin(), elements.end());
		A.makeCompressed();
		if (verbose)timer.Elapse_And_Output_And_Reset("Assemble A to sp_mtx");
	}
	is_A_initialized=true;
}

template<int d> void Projection<d>::Update_b()
{
	int b_size=(int)matrix_to_grid.size();
	#pragma omp parallel for
	for (auto r = 0; r < b_size; r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		real div = (real)0;
		for (int axis = 0; axis < d; axis++) {
			div += (*velocity)(axis, cell + VectorDi::Unit(axis)) - (*velocity)(axis, cell);
		}
		////Attention: use negative div here. We are solving -lap p=-div u
		div_u[r] = -div;
	}
}

template<int d> void Projection<d>::Correction()
{
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			(*velocity)(axis, face) += Velocity_Offset(axis, face);
		}
	}
}

template<int d> void Projection<d>::Update_Mat_Id()
{
	if (solver_mode != SolverType::MULTIGRID_AUTO)return;

	multigrid_params.use_irregular_domain=grid_to_matrix.Has(-1);
	if(!multigrid_params.use_irregular_domain)return;

	is_irregular_domain = true;

	mat_id.Resize(mac_grid->grid.cell_counts,0);
#pragma omp parallel for
	for (int i = 0; i < (int)mat_id.array.size(); i++) {
		if (!Is_Fluid_Cell(mac_grid->grid.Cell_Coord(i))) { mat_id.array[i] = -1; }
	}
}

template<int d> void Projection<d>::Build()
{
	Timer<real> timer;timer.Reset();
	Allocate_System();
	if(verbose)timer.Elapse_And_Output_And_Reset("Allocate A");
	Update_A();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble A");
	Update_b();
	if(verbose)timer.Elapse_And_Output_And_Reset("Assemble b");
}

template<int d>
void Projection<d>::Solve_CPX(void)
{
#ifdef USE_CPX
	if (!cpx_inited) AuxFunc::Crash_With_Info("Projection<d>::Solve_CPX not inited");

	int cell_num = mac_grid->grid.cell_counts.prod();
	cell_b.Resize(mac_grid->grid.cell_counts);
#pragma omp parallel for
	for (int idx = 0; idx < cell_num; idx++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(idx);
		int r = grid_to_matrix(cell);
		if (r == -1)cell_b(cell) = 0;
		else cell_b(cell) = div_u[r];
	}
	cpx_poisson.update_b(cell_b);

	cpx_poisson.solve();

#pragma omp parallel for
	for (auto r = 0; r < matrix_to_grid.size(); r++) {
		const VectorDi& cell = mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		p[r] = cpx_poisson.x(cell);
	}
#else
	AuxFunc::Crash_With_Info("Please compile cpx module with solver mode CPX_GPU");
#endif
}

template<int d> void Projection<d>::Solve()
{
	if (solver_mode == SolverType::AUTO) { Auto_Select_Mode(); }
	if (solver_mode == SolverType::KRYLOV_CPU) {
		KrylovSolver::Params params;
		params.verbose = verbose;
		KrylovSolver::ICPCG(A, p, div_u, params);	////CPU IC-PCG
	}
	else if (solver_mode == SolverType::MULTIGRID_AUTO) {
		MultiGrid::Params multigrid_params;
		multigrid_params.use_auto_calculated_levels = true;
		multigrid_params.dof_on_cell = true;
		multigrid_params.block_size = 1;
		multigrid_params.use_color = multigrid_params.use_gpu;
		multigrid_params.use_irregular_domain = true;
		multigrid_params.use_gpu = true;
		Update_Mat_Id();
#ifdef USE_CUDA
		if (multigrid_params.use_gpu) {
			GMGPCG_GPU<d>(A, p, div_u, mac_grid->grid.cell_counts, mat_id, multigrid_params,/*verbose*/false);
		}
		else { GMGPCG_CPU<d>(A, p, div_u, mac_grid->grid.cell_counts, mat_id, multigrid_params); }
#else
		GMGPCG_CPU<d>(A, p, div_u, mac_grid->grid.cell_counts, mat_id, multigrid_params);
#endif
	}
	else if (solver_mode == SolverType::CPX_GPU) { Solve_CPX(); }
//	if (solver_mode == SolverType::AUTO) {
//		Auto_Select_Mode();
//	}
//
//	if (solver_mode == SolverType::KRYLOV_CPU) {
//		KrylovSolver::Params params;
//		params.verbose = verbose;
//		KrylovSolver::ICPCG(A, p, div_u, params);
//	}
//	else if (solver_mode == SolverType::MULTIGRID_AUTO) {
//		multigrid_params.use_auto_calculated_levels = true;
//		multigrid_params.dof_on_cell = true;
//		multigrid_params.block_size = 1;
//		multigrid_params.use_gpu = true;
//		multigrid_params.init_hier_on_gpu = false;	////calculate hier on CPU to avoid GPU memory crash
//		multigrid_params.use_color = true;
//		Update_Mat_Id();
//
//#ifdef USE_CUDA
//		if (multigrid_params.use_gpu) {
//			gmg_solver_gpu.update_A_levels = update_A;
//			gmg_solver_gpu.Initialize(A, mac_grid->grid.cell_counts, multigrid_params, is_irregular_domain ? &mat_id : nullptr);
//			gmg_solver_gpu.Solve(p, div_u);
//		}
//		else {
//			gmg_solver_cpu.update_A_levels = update_A;
//			gmg_solver_cpu.Initialize(A, mac_grid->grid.cell_counts, multigrid_params, is_irregular_domain ? &mat_id : nullptr);
//			gmg_solver_cpu.Solve(p, div_u);
//		}
//#else
//		gmg_solver_cpu.Initialize(A, mac_grid->grid.cell_counts, multigrid_params, &mat_id);
//		gmg_solver_cpu.Solve(p, div_u);
//#endif
//	}
//	else if (solver_mode == SolverType::CPX_GPU) {
//		Solve_CPX();
//	}
}

template<int d> void Projection<d>::Project()
{
	Build();
	Solve();
	Correction();
}

template<int d> void Projection<d>::Clear()
{
	p.fill((real)0);
}

//////////////////////////////////////////////////////////////////////////
////Check functions
template<int d> void Projection<d>::Pressure(Field<real,d>& pressure) const
{
	pressure.Resize(mac_grid->grid.cell_counts,(real)0);
	iterate_cell(iter,mac_grid->grid){const VectorDi& cell=iter.Coord();
		int idx=grid_to_matrix(cell);if(idx==-1)continue;
		pressure(cell)=p[idx];}
}

template<int d> void Projection<d>::Pressure_Gradient(FaceField<real,d>& grad_p) const
{
	grad_p.Fill((real)0);
	for (int axis = 0; axis < d; axis++) {
		int face_num = mac_grid->face_grids[axis].node_counts.prod();
#pragma omp parallel for
		for (int i = 0; i < face_num; i++) {
			VectorDi face = mac_grid->face_grids[axis].Node_Coord(i);
			grad_p(axis, face) = -Velocity_Offset(axis, face);
		}
	}
}

template<int d> void Projection<d>::Divergence(Field<real,d>& div) const
{
	div.Resize(mac_grid->grid.cell_counts,(real)0);
	int b_size=(int)matrix_to_grid.size();
	#pragma omp parallel for
	for(auto r=0;r<b_size;r++){
		const VectorDi& cell=mac_grid->grid.Cell_Coord(matrix_to_grid[r]);
		real divg=(real)0;
		for(int axis=0;axis<d;axis++){divg+=((*velocity)(axis,cell+VectorDi::Unit(axis))-(*velocity)(axis,cell));}
		div(cell)=divg;}
}

//////////////////////////////////////////////////////////////////////////
////Physical interface functions that defines the problem
//see: https://wmdcstdio.com/2021/07/11/projection-matrix-terms/
template<int d>
real Projection<d>::Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx) const
{
	VectorDi nb_cell = Grid<d>::Nb_C(fluid_cell, nbidx);
	int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(fluid_cell, nbidx, axis, face);
	if (bc->Is_Psi_N(axis, face))return 0;
	if (Is_Fluid_Cell(nb_cell)) return -1;
	return 0;
}

template<int d>
real Projection<d>::Diag_Face_Term(const int& axis, const VectorDi& face) const
{
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	if (!Is_Fluid_Cell(cell[0])) std::swap(cell[0], cell[1]);
	if (!Is_Fluid_Cell(cell[0])) return (real)0;
	if (bc->Is_Psi_N(axis, face)) return (real)0;
	return 1;
}

template<int d>
real Projection<d>::Velocity_Offset(const int& axis, const VectorDi& face) const
{
	if (bc->Is_Psi_N(axis, face)) return 0;
	VectorDi cell[2]; for (int i = 0; i < 2; i++)cell[i] = MacGrid<d>::Face_Incident_Cell(axis, face, i);
	real cell_p[2]; for (int i = 0; i < 2; i++)cell_p[i] = (Is_Valid_Cell(cell[i]) && Is_Fluid_Cell(cell[i])) ? p[grid_to_matrix(cell[i])] : (real)0;
	return -(cell_p[1] - cell_p[0]);
}


template class Projection<2>;
template class Projection<3>;
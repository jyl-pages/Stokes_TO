//////////////////////////////////////////////////////////////////////////
// Poisson solver on a irregular domain
// Copyright (c) (2018-), Jinyuan Liu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __PoissonIrregularOld_h__
#define __PoissonIrregularOld_h__
#include "Common.h"
#include "Field.h"
#include "MacGrid.h"
#include "FaceField.h"
#include "Interpolation.h"
#include "KrylovSolver.h"
#include "Hashtable.h"
#include "LevelSet.h"
#include "Timer.h"

#include "GmgPcgSolverCPU.h"
#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif


template<int d,class TC=real> class PoissonIrregularOld
{Typedef_VectorDii(d);
public:
	MacGrid<d> mac_grid;

	////data
	FaceField<TC,d>* alpha=nullptr;							////face alpha
	Field<ushort,d>* type=nullptr;							////cell type
	BoundaryConditionMacGrid<d>* bc=nullptr;				////boundary condition, this is for regular bcs
	LevelSet<d>* levelset=nullptr;							////free boundary levelset
	Field<TC,d> rhs;										////rhs
	Field<TC,d> p;											////solution
	bool own_data=false;

	////options
	bool verbose=false;
	bool initialized=false;

	////irregular boundary conditions
	std::function<real(const VectorD&)> psi_D_irregular=nullptr;
	real eps;

	////matrix implementation
	SparseMatrix<TC> A;
	VectorN<TC> x;			////solution stored in VectorN
	VectorN<TC> b;
	int nc=0;
	Field<int,d> grid_to_matrix;
	Array<int> matrix_to_grid;
	
	////multigrid solver and GPU parallelization
	bool use_multigrid_solver=true;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif

	PoissonIrregularOld(){}
	~PoissonIrregularOld(){if(own_data){if(alpha)delete alpha;if(type)delete type;if(bc)delete bc;if(levelset)delete levelset;}}

	void Initialize(const MacGrid<d>& _mac_grid)
	{
		mac_grid=_mac_grid;
		eps=pow(mac_grid.grid.dx,2);
		alpha=new FaceField<real,d>(mac_grid.grid.cell_counts,(real)1);
		type=new Field<ushort,d>(mac_grid.grid.cell_counts,(ushort)CellType::Fluid);
		bc=new BoundaryConditionMacGrid<d>(mac_grid);
		levelset=new LevelSet<d>();
		levelset->Initialize(mac_grid.grid);
		own_data=true;
		Allocate_Data();
		initialized=true;
	}

	void Initialize(MacGrid<d>& _mac_grid,FaceField<TC,d>* _alpha,Field<ushort,d>* _type,BoundaryConditionMacGrid<d>* _bc,LevelSet<d>* _levelset)
	{
		mac_grid=_mac_grid;
		eps=pow(mac_grid.grid.dx,2);
		alpha=_alpha;type=_type;bc=_bc;levelset=_levelset;own_data=false;
		Allocate_Data();
		initialized=true;
	}

	void Allocate_Data()
	{
		rhs.Resize(mac_grid.grid.cell_counts);rhs.Fill((TC)0);
	}

	//////////////////////////////////////////////////////////////////////////
	////matrix implementation

	void Update_A()
	{
		Timer<real> timer;timer.Reset();

		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			if(Is_Inside_Levelset_eps(cell)){(*type)(cell)=(ushort)CellType::Fluid;} 
			else{(*type)(cell)=(ushort)CellType::Air;}}	////TOIMPL: solid

		////not parallizable
		std::function<bool(const int)> fluid_cell = [=](const int idx)->bool{return this->Is_Fluid_Cell(idx);}; 
		Build_Grid_Cell_Matrix_Bijective_Mapping(mac_grid.grid,fluid_cell,grid_to_matrix,matrix_to_grid);

		nc=(int)matrix_to_grid.size();A.resize(nc,nc);

		Array<TripletT> elements;
		for(auto r=0;r<matrix_to_grid.size();r++){
			const VectorDi& cell=mac_grid.grid.Cell_Coord(matrix_to_grid[r]);
			////off-diagonal elements
			if(!bc->Is_Psi_D(cell))for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
				VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(Is_Valid_Cell(nb_cell)&&Is_Fluid_Cell(nb_cell)&&!bc->Is_Psi_D(nb_cell)){
					int c=grid_to_matrix(nb_cell);
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(i,axis,side);
					VectorDi face=cell+VectorDi::Unit(axis)*side;
					real a=(*alpha)(axis,face);
					elements.push_back(TripletT((int)r,(int)c,(real)-a));}}
			////diagonal elements
			real dia_coef=(real)0;
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
				VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(!Is_Valid_Cell(nb_cell))continue;
				
				real coef=Dia_Coef(cell,i);
				real phi_cell = levelset->phi(cell); 
				real phi_nb = levelset->phi(nb_cell);
				if(psi_D_irregular != nullptr && levelset->Interface(phi_cell, phi_nb)) {
					real theta = levelset->Theta(phi_cell, phi_nb);
					coef /= theta;}
				dia_coef+=coef;} 
			elements.push_back(TripletT((int)r,(int)r,dia_coef));}
		
		if(verbose)timer.Elapse_And_Output_And_Reset("Update A elements");
		A.setFromTriplets(elements.begin(),elements.end()); 
		A.makeCompressed(); 
		if(verbose)timer.Elapse_And_Output_And_Reset("Assemble A to sp_mtx");
	}

	void Update_Rhs_With_BC()
	{
		x.resize(nc);x.fill((TC)0);b.resize(nc);b.fill((TC)0);
		real dx=mac_grid.grid.dx;
		real dx2=(TC)pow(mac_grid.grid.dx,2);

		////Set b
		for(auto r=0;r<matrix_to_grid.size();r++){	
			const VectorDi& cell=mac_grid.grid.Cell_Coord(matrix_to_grid[r]);
			b[r]=-dx2*rhs(cell);

			////Set irregular boundary
			if (psi_D_irregular == nullptr) {continue;}
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
				VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(!Is_Valid_Cell(nb_cell))continue;

				real phi_cell = levelset->phi(cell); 
				real phi_nb = levelset->phi(nb_cell);
				if(levelset->Interface(phi_cell, phi_nb)) {
					real theta = levelset->Theta(phi_cell, phi_nb);
					real boundary_val = psi_D_irregular(mac_grid.grid.Center(nb_cell));
					b[r] += boundary_val / theta;}}}

		////Set Dirichlet boundary rhs
		for(const auto& pD:(*bc).psi_D_values){const int cell_index=pD.first;const real val=pD.second.second;
			////psi_D's rhs
			const VectorDi cell=mac_grid.grid.Cell_Coord(cell_index);

			int r=grid_to_matrix(cell);b[r]=A.coeff(r,r)*val;
			////psi_D nbs' rhs
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(Is_Valid_Cell(nb_cell)&&Is_Fluid_Cell(nb_cell)&&!(*bc).Is_Psi_D(nb_cell)){
					int r=grid_to_matrix(nb_cell);
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(i,axis,side);
					VectorDi face=cell+VectorDi::Unit(axis)*side;real a=(*alpha)(axis,face);b[r]+=a*val;}}}

		////Set Neumann boundary rhs
		for(const auto& pN:(*bc).psi_N_values){
			int axis=pN.first[0];VectorDi face=mac_grid.face_grids[axis].Node_Coord(pN.first[1]);real value=pN.second;
			const VectorDi cell_0=MacGrid<d>::Face_Incident_Cell(axis,face,0);
			const VectorDi cell_1=MacGrid<d>::Face_Incident_Cell(axis,face,1);
			if(Is_Valid_Cell(cell_0)&&Is_Fluid_Cell(cell_0)){
				int r0=grid_to_matrix(cell_0);
				if(r0!=-1)b[r0]+=(*alpha)(axis,face)*value*dx;}
			if(Is_Valid_Cell(cell_1)&&Is_Fluid_Cell(cell_1)){
				int r1=grid_to_matrix(cell_1);
				if(r1!=-1)b[r1]-=(*alpha)(axis,face)*value*dx;}}
	}

	//////////////////////////////////////////////////////////////////////////
	////solve

	void Build()
	{
		Update_A();
		Update_Rhs_With_BC();
	}

	void Solve()
	{
		use_multigrid_solver=false;
		if(use_multigrid_solver){
			GeometricMultiGrid::Params multigrid_params;
			multigrid_params.use_auto_calculated_levels=true;
			multigrid_params.dof_on_cell=true;
			multigrid_params.block_size=1;
			multigrid_params.use_gpu=true;

			#ifdef USE_CUDA
			if(multigrid_params.use_gpu){
				gmg_solver_gpu.update_A_levels=false;
				gmg_solver_gpu.Initialize(A,mac_grid.grid.cell_counts,multigrid_params);
				gmg_solver_gpu.Solve(x,b);
			}
			else{ 
				gmg_solver_cpu.update_A_levels=true;
				gmg_solver_cpu.Initialize(A,mac_grid.grid.cell_counts,multigrid_params);
				gmg_solver_cpu.Solve(x,b);}
			#else
				gmg_solver_cpu.Initialize(A,mac_grid.grid.cell_counts,multigrid_params);
				gmg_solver_cpu.Solve(x,b);
			#endif
		}
		else{KrylovSolver::ICPCG(A,x,b);}

		p.Resize(mac_grid.grid.cell_counts,(TC)0);
		iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(!Is_Fluid_Cell(cell))continue;
			if(bc->Is_Psi_D(cell))p(cell)=bc->Psi_D_Value(cell);
			else p(cell)=x[grid_to_matrix(cell)];}
	}

	void Build_And_Solve()
	{
		Build();
		Solve();
	}

	////Helper functions
	real Dia_Coef(const VectorDi& cell,const int i)
	{int axis;VectorDi face;MacGrid<d>::Cell_Incident_Face(cell,i,axis,face);
	if((*bc).Is_Psi_N(axis,face))return (real)0;else return (*alpha)(axis,face);}

	inline bool Is_Valid_Cell(const VectorDi& cell) const {return mac_grid.grid.Valid_Cell(cell);}
	inline bool Is_Valid_Cell(const int cell_idx) const {return mac_grid.grid.Valid_Cell(mac_grid.grid.Cell_Coord(cell_idx));}
	inline bool Is_Fluid_Cell(const VectorDi& cell) const {return (*type)(cell)==(ushort)CellType::Fluid;}
	inline bool Is_Fluid_Cell(const int cell_idx) const {return (*type).array[cell_idx]==(ushort)CellType::Fluid;}

	inline bool Is_Inside_Levelset_eps(const VectorD& pos) const 
	{real phi=(*levelset).Phi(pos);return phi<=-(real)eps*mac_grid.grid.dx;}
	inline bool Is_Inside_Levelset_eps(const VectorDi& cell) const 
	{real phi=(*levelset).phi(cell);return phi<=-(real)eps*mac_grid.grid.dx;}
	inline bool Is_Inside_Levelset(const VectorD& pos) const 
	{real phi=(*levelset).Phi(pos);return phi<=(real)0;}
	inline bool Is_Inside_Levelset(const VectorDi& cell) const 
	{real phi=(*levelset).phi(cell);return phi<=(real)0;}
};
#endif

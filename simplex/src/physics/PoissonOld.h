//////////////////////////////////////////////////////////////////////////
// Poisson solver on a MAC grid regular domain
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __PoissonOld_h__
#define __PoissonOld_h__
#include "Common.h"
#include "Field.h"
#include "MacGrid.h"
#include "FaceField.h"
#include "Interpolation.h"
#include "KrylovSolver.h"
#include "Timer.h"
#include "Hashtable.h"
#include "BoundaryCondition.h"
#include "GmgPcgSolverCPU.h"
#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif


template<int d,class TC=real> class PoissonOld
{Typedef_VectorDii(d);
public:
	MacGrid<d>& mac_grid;

	//////////////////////////////////////////////////////////////////////////
	////data
	FaceField<TC,d>& alpha;											////alpha
	Field<TC,d>& rhs;												////rhs
	Field<TC,d>& p;													////solution
	BoundaryConditionMacGrid<d>& bc;								////boundary condition
	bool own_data=false;

	//////////////////////////////////////////////////////////////////////////
	////matrix implementation
	SparseMatrix<TC> A;
	VectorN<TC> x;			////solution stored in VectorN
	VectorN<TC> b;
	int nc=0;
	Field<int,d> grid_to_matrix;
	Array<int> matrix_to_grid;

	//////////////////////////////////////////////////////////////////////////
	////multigrid solver and GPU parallelization
	bool use_multigrid_solver=true;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif

	//////////////////////////////////////////////////////////////////////////
	////constructors
	PoissonOld():mac_grid(*(new MacGrid<d>())),alpha(*(new FaceField<TC,d>())),rhs(*(new Field<TC,d>())),p(*(new Field<TC,d>())),bc(*(new BoundaryConditionMacGrid<d>(mac_grid))){own_data=true;}
	PoissonOld(MacGrid<d>& _mac_grid,FaceField<TC,d>& _alpha,Field<TC,d>& _rhs,Field<TC,d>& _p,BoundaryConditionMacGrid<d>& _bc)
		:mac_grid(_mac_grid),alpha(_alpha),rhs(_rhs),p(_p),bc(_bc){own_data=false;}
	~PoissonOld()
	{
		if(own_data){
			auto* mac_grid_ptr=&mac_grid;delete mac_grid_ptr;
			auto* alpha_ptr=&alpha;delete alpha_ptr;
			auto* rhs_ptr=&rhs;delete rhs_ptr;
			auto* p_ptr=&p;delete p_ptr;
			auto* bc_ptr=&bc;delete bc_ptr;}
	}

	//////////////////////////////////////////////////////////////////////////
	////initialization data
	void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts,dx,domain_min);
		alpha.Resize(mac_grid.grid.cell_counts);alpha.Fill((TC)1);
		rhs.Resize(mac_grid.grid.cell_counts);rhs.Fill((TC)0);
	}
	void Initialize(const Grid<d>& _grid){Initialize(_grid.cell_counts,_grid.dx,_grid.domain_min);}

	//////////////////////////////////////////////////////////////////////////
	////assembling matrix 
	void Update_A()
	{
		grid_to_matrix.Resize(mac_grid.grid.cell_counts,-1);matrix_to_grid.clear();
		nc=0;iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(mac_grid.grid.Valid_Cell(cell)){grid_to_matrix(cell)=nc++;
				matrix_to_grid.push_back(mac_grid.grid.Cell_Index(cell));}}	

		A.resize(nc,nc);
		TC dx=mac_grid.grid.dx;real dx2=(TC)pow(mac_grid.grid.dx,2);
		Array<Triplet<TC> > elements;
		for(auto r=0;r<matrix_to_grid.size();r++){	
			const VectorDi& cell=mac_grid.grid.Cell_Coord(matrix_to_grid[r]);
			////nb elements in the row
			if(!bc.Is_Psi_D(cell))for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(mac_grid.grid.Valid_Cell(nb_cell)&&!bc.Is_Psi_D(nb_cell)){int c=grid_to_matrix(nb_cell);
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(i,axis,side);
					VectorDi face=cell+VectorDi::Unit(axis)*side;
					TC a=alpha(axis,face);elements.push_back(Triplet<TC>((int)r,(int)c,-a));}}
			////diagonal element
			TC dia_coef=(TC)0;
			for(int axis=0;axis<d;axis++)for(int i=0;i<2;i++){////for Neumann boundary
				VectorDi face=mac_grid.Cell_Incident_Face(axis,cell,i);
				if(!bc.Is_Psi_N(axis,face)){TC a=alpha(axis,face);dia_coef+=a;}}
			elements.push_back(Triplet<TC>((int)r,(int)r,dia_coef));}
		A.setFromTriplets(elements.begin(),elements.end());A.makeCompressed();			
	}

	void Update_Rhs_With_BC()
	{
		x.resize(nc);x.fill((TC)0);b.resize(nc);b.fill((TC)0);
		real dx=mac_grid.grid.dx;
		real dx2=(TC)pow(mac_grid.grid.dx,2);

		////set b
		for(auto r=0;r<matrix_to_grid.size();r++){	
			const VectorDi& cell=mac_grid.grid.Cell_Coord(matrix_to_grid[r]);
			b[r]=-dx2*rhs(cell);}

		////set Dirichlet boundary rhs
		for(const auto& pD:bc.psi_D_values){const int cell_index=pD.first;const real val=pD.second.second;
			////psi_D's rhs
			const VectorDi cell=mac_grid.grid.Cell_Coord(cell_index);
			int r=grid_to_matrix(cell);b[r]=A.coeff(r,r)*val;
			////psi_D nbs' rhs
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(mac_grid.grid.Valid_Cell(nb_cell)&&!bc.Is_Psi_D(nb_cell)){int r=grid_to_matrix(nb_cell);
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(i,axis,side);
					VectorDi face=cell+VectorDi::Unit(axis)*side;real a=alpha(axis,face);b[r]+=a*val;}}}

		////set Neumann boundary rhs
		for(const auto& pN:bc.psi_N_values){
			int axis=pN.first[0];VectorDi face=mac_grid.face_grids[axis].Node_Coord(pN.first[1]);real value=pN.second;
			const VectorDi cell_0=MacGrid<d>::Face_Incident_Cell(axis,face,0);
			const VectorDi cell_1=MacGrid<d>::Face_Incident_Cell(axis,face,1);
			if(mac_grid.grid.Valid_Cell(cell_0)){
				int r0=grid_to_matrix(cell_0);
				if(r0!=-1)b[r0]+=alpha(axis,face)*value*dx;}
			if(mac_grid.grid.Valid_Cell(cell_1)){
				int r1=grid_to_matrix(cell_1);
				if(r1!=-1)b[r1]-=alpha(axis,face)*value*dx;}}
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
			if(bc.Is_Psi_D(cell))p(cell)=bc.Psi_D_Value(cell);
			else p(cell)=x[grid_to_matrix(cell)];}
	}

	void Build_And_Solve()
	{
		Build();
		Solve();
	}
};
#endif
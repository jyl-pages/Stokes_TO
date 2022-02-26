//////////////////////////////////////////////////////////////////////////
// Linear Grid FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyNonlinearFemGrid_h__
#define __SoftBodyNonlinearFemGrid_h__
#include "Hashtable.h"
#include "LinearFemFunc.h"
#include "BoundaryCondition.h"
#include "Field.h"
#include "GmgPcgSolverCPU.h"
#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif

template<int d> class SoftBodyNonlinearFemGrid
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	Grid<d> grid;
	SparseMatrixT K;
	VectorX u;
	VectorX f;

	Array<ElasticParam> materials;
	Field<short,d> material_id;
	Field<real,d>* variable_coef=nullptr;
	BoundaryConditionGrid<d> bc;
	
	////Newton solver
	VectorX f_elas;
	VectorX e_elas;
	int max_iter=50;
	real rel_tol=(real)1e-5;
	real abs_tol=(real)1e-10;

	////multigrid and parallelization
	bool use_multigrid_solver=true;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif
	Array<int> colored_cell_ptr;
	Array<int> colored_cell_indices;
	
	void Initialize(const Grid<d>& _grid);
	void Allocate_K();
	void Update_K();
	void Set_RHS_Force(VectorX& _f);
	void Set_K_And_RHS_Psi_D(VectorX& _f);
	void Solve_Nonlinear();
	void Solve_Linear(VectorX& _u,const VectorX& _f);
	
	void Add_Material(real youngs,real poisson);
	void Clear_Materials(){materials.clear();}

	void Set_Fixed(const VectorDi& node);
	void Set_Displacement(const VectorDi& node,const VectorD& dis);
	void Set_Force(const VectorDi& node,const VectorD& force);

	void Compute_Cell_Displacement(const VectorX& u,const VectorDi& cell,VectorX& cell_u) const;
	void Compute_Cell_Displacement(const VectorX& u,const Array<int>& cell_node_matrix_indices,VectorX& cell_u) const;

protected:
	int Node_Index_In_K(const VectorDi& node) const;
};

#endif

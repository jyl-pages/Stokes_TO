//////////////////////////////////////////////////////////////////////////
// Linear Grid FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyLinearFemGrid_h__
#define __SoftBodyLinearFemGrid_h__
#include "Hashtable.h"
#include "LinearFemFunc.h"
#include "BoundaryCondition.h"
#include "Field.h"
#include "GmgPcgSolverCPU.h"
#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif

template<int d> class SoftBodyLinearFemGrid
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	Grid<d> grid;
	SparseMatrixT K;
	VectorX u;
	VectorX f;

	Array<MatrixX> E;
	Array<MatrixX> Ke;
	Field<short,d> material_id;
	Field<real,d>* variable_coef=nullptr;
	BoundaryConditionGrid<d> bc;
	
	////multigrid and parallelization
	bool use_multigrid_solver=true;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif

	Array<int> colored_cell_ptr;
	Array<int> colored_cell_indices;
	
	//SoftBodyLinearFemGrid();

	virtual void Initialize(const Grid<d> _grid);
	void Allocate_K();
	virtual void Update_K_And_f();
	virtual void Solve();

	void Add_Material(real youngs,real poisson);
	//void Add_Material(MatrixX E);
	void Clear_Materials(){Ke.clear();}

	virtual void Set_Fixed(const VectorDi& node);
	void Set_Displacement(const VectorDi& node,const VectorD& dis);
	virtual void Set_Force(const VectorDi& node,const VectorD& force);

	void Compute_Cell_Displacement(const VectorX& u,const VectorDi& cell,VectorX& cell_u) const;
	void Compute_Cell_Displacement(const VectorX& u,const Array<int>& cell_node_matrix_indices,VectorX& cell_u) const;
	virtual void Compute_Elastic_Compliance(Field<real, d>& energy);
	void Compute_Elastic_Energy(Field<real, d>& energy) const;
    void Compute_Strain(Field<MatrixD,d>& strains) const;
	void Compute_Stress(Field<MatrixD,d>& stresses) const;
	void Compute_Von_Mises_Stress(Field<real,d>& von_mises,const Field<MatrixD,d>* stresses=nullptr) const;
	void Compute_Strain(VectorX& strain,const VectorDi& cell) const;
	void Compute_Strain(MatrixD& strain,const VectorDi& cell) const;
	void Compute_Stress(VectorX& stress,const VectorDi& cell,const MatrixX& E) const;
	void Compute_Stress(MatrixD& stress,const VectorDi& cell,const MatrixX& E) const;
protected:
	int Node_Index_In_K(const VectorDi& node) const;
	void Compute_Cell_Tensor_Helper(const Array<MatrixX>& B,const VectorX& u,Field<MatrixD,d>& tensors) const;
};



#endif

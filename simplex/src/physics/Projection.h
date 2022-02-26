//////////////////////////////////////////////////////////////////////////
// Project a vector field to divergence free on a MAC grid
// Copyright (c) (2018-), Bo Zhu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

//////////////////////////////////////////////////////////////////////////
//// This projection solver currently only takes zero Neumann bc. The velocities on the Neumann boundary need to be set to the correct values before the projection.
//// The calculated pressure value is the real pressure value divided by delta_x.
//////////////////////////////////////////////////////////////////////////

#ifndef __Projection_h__
#define __Projection_h__
#include <functional>
#include "MacGrid.h"
#include "FaceField.h"
#include "Field.h"
#include "BoundaryCondition.h"
#include "GmgPcgSolverCPU.h"
#include "TypeFunc.h"

#ifdef USE_CUDA
#include "GmgPcgSolverGPU.h"
#endif

#ifdef USE_CPX
#include "Poisson.h"
#endif

template<int d> class Projection
{Typedef_VectorDii(d);
public:
	////data
	MacGrid<d>* mac_grid;
	FaceField<real,d>* velocity=nullptr;
	Field<ushort,d>* type=nullptr;
	BoundaryConditionMacGrid<d>* bc=nullptr;
	bool own_grid=false;
	bool own_velocity=false;
	bool own_type=false;
	bool own_bc=false;

	////flags
	bool verbose=true;

	////linear solver
	Field<int,d> grid_to_matrix;
	Array<int> matrix_to_grid;
	SparseMatrixT A;
	VectorX p;				////unknown
	VectorX div_u;			////rhs: the solver solves Ap=div_u
	bool is_A_initialized=false;
	bool update_A=true;

	////multigrid solver
	Field<short, d> mat_id;
	bool is_irregular_domain = false;
	MultiGrid::Params multigrid_params;
	GMGPCG_Solver_CPU<d> gmg_solver_cpu;
	#ifdef USE_CUDA
	GMGPCG_Solver_GPU<real,d> gmg_solver_gpu;
	#endif

	////multigrid solver
#ifdef USE_CPX
	//using CPX_Solver = typename If<d == 2, Poisson, Poisson3D >::Type;
	//CPX_Solver cpx_poisson;
	Poisson<d> cpx_poisson;
	bool cpx_inited = false;
	FaceField<Scalar, d> face_vol;
	Field<int, d> cell_fixed;
	Field<Scalar, d> cell_b;
#endif

protected:
	SolverType solver_mode = SolverType::AUTO;

public:

	////constructors
	Projection(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, Field<ushort, d>* _type = nullptr, BoundaryConditionMacGrid<d>* _bc = nullptr, const SolverType& _mode = SolverType::AUTO);
	~Projection();

	virtual void Initialize(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, Field<ushort, d>* _type, BoundaryConditionMacGrid<d>* _bc, const SolverType& _mode);

	////set attributes
	void Auto_Select_Mode(void);
	void Set_Velocity(FaceField<real,d>& _velocity){if(velocity!=nullptr&&own_velocity)delete velocity;velocity=&_velocity;own_velocity=false;}
	void Set_Type(Field<ushort,d>& _type){if(type!=nullptr&&own_type)delete type;type=&_type;own_type=false;}
	void Set_BC(BoundaryConditionMacGrid<d>& _bc){if(bc!=nullptr&&own_bc)delete bc;bc=&_bc;own_bc=false;}

	////projection functions
	virtual void Allocate_System();
	virtual void Prepare_CPX_System(void);
	virtual void Update_A();
	virtual void Update_b();				////calculate b as div velocity
	virtual void Correction();
	void Update_Mat_Id();			////for irregular domain
	virtual void Build();					////call allocate, update_A, and update_b
	void Solve_CPX(void);
	virtual void Solve();
	virtual void Project();					////call both build, solve, and correction
	void Clear();

	////read data
	void Pressure(Field<real,d>& pressure) const;				////write values of p into pressure
	void Pressure_Gradient(FaceField<real,d>& grad_p) const;	////write values of p into pressure
	void Divergence(Field<real,d>& div) const;					////write values of velocity into div

	////Physical interface functions that defines the problem
	//NOTE: if you want to design a derived class of this, theoretically you only need to implement these 5 functions
	virtual real Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx)const;
	virtual real Diag_Face_Term(const int& axis, const VectorDi& face)const;
	virtual real Velocity_Offset(const int& axis, const VectorDi& face)const;
	virtual bool Is_Valid_Cell(const VectorDi& cell) const {return mac_grid->grid.Valid_Cell(cell);}
	virtual bool Is_Fluid_Cell(const VectorDi& cell) const { return Is_Valid_Cell(cell) && (*type)(cell) == (ushort)CellType::Fluid; }
};

#endif

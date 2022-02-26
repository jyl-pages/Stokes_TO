//////////////////////////////////////////////////////////////////////////
// Project a vector field with free boundary to divergence free on a MAC grid
// Copyright (c) (2018-), Fan Feng, Shuqi Yang, Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __ProjectionIrregular_h__
#define __ProjectionIrregular_h__

#include "Projection.h"
#include "LevelSet.h"

template<int d> class ProjectionIrregular : public Projection<d>
{
	Typedef_VectorDii(d);
	using Base = Projection<d>;
public:
	////data
	LevelSet<d>* levelset = nullptr;
	bool own_levelset = false;

	////flags
	bool use_levelset_interface = true;		////do not access levelset for free boundary if false 
	bool use_surface_tension = false;

	////jump condition (surface tension)
	bool use_jump_condition = false;
	real sigma = (real)1e-5;					////surface tension coefficient for the default pressure jump
	real current_dt = (real)1;				////need to set dt when using the jump condition because dt is absorbed in p when calculating the projection

	////callback functions
	std::function<bool(const int)> Is_Fluid_Interior_Cell_Index = nullptr;
	std::function<real(const VectorD&)> Jump_Condition = std::bind(&ProjectionIrregular<d>::Pressure_Jump, this, std::placeholders::_1);
	std::function<bool(const std::pair<int, int>)> Is_Interface_Face_Index = std::bind(&ProjectionIrregular<d>::Is_Levelset_Interface_Face_Index, this, std::placeholders::_1);
	std::function<real(real)> Dirac = std::bind(&ProjectionIrregular<d>::Levelset_Dirac, this, std::placeholders::_1);;

	////linear solver
	FaceField<int, d> macgrid_to_matrix;
	Array<std::pair<int, int>> matrix_to_macgrid;
	bool use_update_A_in_parallel = false;	////needs to specify the Is_Fluid_Interior_Cell_Index callback before setting this flag

	////divergence control
	bool use_vol_control = false;
	bool calc_current_vol = true;				////calculate the current vol within Apply_Vol_Control_To_b or not
	real target_vol = (real)-1;				////ALWAYS NEEDS TO BE SET EXTERNALLY!
	real current_vol = (real)-1;				////needs to be set externally or calculated within Apply_Vol_Control_To_b when calc_current_vol=true			
	real vol_control_ks = (real)1e2;


	// default narrow band width
	int narrow_band_cell_num = 5;
	int dirac_band_cell_num = 3;


public:
	////constructors
	ProjectionIrregular(MacGrid<d>& _mac_grid, FaceField<real, d>& _velocity, LevelSet<d>& _levelset, BoundaryConditionMacGrid<d>& _bc, const SolverType& _mode = SolverType::AUTO);
	virtual void Initialize(LevelSet<d>* _levelset);

	////set attributes
	void Set_Level_Set(LevelSet<d>& _levelset) { if (levelset != nullptr && own_levelset)delete levelset; levelset = &_levelset; own_levelset = false; }

	////projection functions
	virtual void Update_A();
	void Update_A_In_Parallel();
	void Apply_Jump_Condition_To_b();
	void Apply_Vol_Control_To_b();
	void Apply_Implicit_Surface_tension(const real dt);
	virtual void Build();
	virtual void Project();


	////interface coef returns clamped 1/theta
	real Intf_Coef(const real& phi0, const real& phi1, real& theta, bool& is_intf) const;				////levelset-related
	//////////////////////////////////////////////////////////////////////////
	//[WARNING] Deprecated functions, planning to remove in next versions
	real Intf_Coef(const VectorDi& cell,/*nb idx*/const int i, real& theta, bool& is_intf) const;		////levelset-related
	//////////////////////////////////////////////////////////////////////////

	int Fluid_Neighbor_Number_Index(const int& fluid_cell_idx)const;

	//////////////////////////////////////////////////////////////////////////
	////default callback functions: all levelset-related

	inline real Pressure_Jump(const VectorD& pos) const { real curvature = (*levelset).Curvature(pos); return current_dt * sigma * curvature; }

	inline bool Is_Levelset_Interface_Face(const int axis, const VectorDi& face) const
	{
		VectorD pos = mac_grid->Face_Center(axis, face);
		if ((*bc).Is_Psi_N(axis, face)) return false;
		real phi = (*levelset).Phi(pos);
		return (phi > -(real)narrow_band_cell_num * mac_grid->grid.dx && phi < (real)narrow_band_cell_num* mac_grid->grid.dx);
	}

	inline bool Is_Levelset_Interface_Face_Index(const std::pair<int, int> face_idx) const
	{
		int axis = face_idx.first;
		VectorDi face = mac_grid->face_grids[axis].Node_Coord(face_idx.second);
		VectorD pos = mac_grid->Face_Center(axis, face);
		if ((*bc).Is_Psi_N(axis, face)) return false;
		real phi = (*levelset).Phi(pos);
		return (phi > -(real)narrow_band_cell_num * mac_grid->grid.dx && phi < (real)narrow_band_cell_num* mac_grid->grid.dx);
	}

	inline real Levelset_Dirac(const real phi) const
	{
		if (phi < -(real)dirac_band_cell_num * mac_grid->grid.dx) return 0;
		else if (phi > (real)dirac_band_cell_num * mac_grid->grid.dx) return 0;
		else return 0.5 * (1.0 + cos(pi * phi / (real)dirac_band_cell_num / mac_grid->grid.dx)) / ((real)dirac_band_cell_num * mac_grid->grid.dx);
	}

	////Physical interface functions that defines the problem
	virtual real Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx)const;
	virtual real Diag_Face_Term(const int& axis, const VectorDi& face)const;
	virtual real Velocity_Offset(const int& axis, const VectorDi& face)const;
	//Is_Valid_Cell: same as Projection<d>
	virtual bool Is_Fluid_Cell(const VectorDi& cell) const { return mac_grid->grid.Valid_Cell(cell) && ((*levelset).phi(cell) < (real)0); }
};
#endif
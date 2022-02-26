/////////////////////////////////////////////////////////////////////////
// Projection with coefficient face field alpha
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#pragma once

#include "Projection.h"

template<int d> class ProjectionMeta : public Projection<d>
{
	Typedef_VectorDii(d);
	using Base = Projection<d>;
public:
	////data
	FaceField<real, d>* alpha = nullptr;
	bool own_alpha = false;

	////flags
	bool use_alpha_for_correction = false;
	bool use_alpha_for_update_b = true;

	//Initialization
	ProjectionMeta(MacGrid<d>* _mac_grid, FaceField<real, d>* _velocity, FaceField<real, d>* _alpha, Field<ushort, d>* _type = nullptr, BoundaryConditionMacGrid<d>* _bc = nullptr, const SolverType& _mode = SolverType::AUTO);
	ProjectionMeta(MacGrid<d>& _mac_grid, FaceField<real, d>& _velocity, const SolverType& _mode = SolverType::AUTO);
	ProjectionMeta(MacGrid<d>& _mac_grid, FaceField<real, d>& _velocity, FaceField<real, d>& _alpha, Field<ushort, d>& _type, BoundaryConditionMacGrid<d>& _bc, const SolverType& _mode = SolverType::AUTO);
	~ProjectionMeta();
	virtual void Initialize(FaceField<real, d>* _alpha);

	//Set attributes
	void Set_Alpha(FaceField<real, d>& _alpha) { if (alpha != nullptr && own_alpha)delete alpha; alpha = &_alpha; own_alpha = false; }

	////projection functions
	virtual void Update_b();				////calculate b as div velocity

	////Physical interface functions that defines the problem
	virtual real Off_Diag_Term(const VectorDi& fluid_cell, const int& nbidx)const;
	virtual real Diag_Face_Term(const int& axis, const VectorDi& face)const;
	virtual real Velocity_Offset(const int& axis, const VectorDi& face)const;
	//Is_Valid_Cell: same as Projection
	//Is_Fluid_Cell: same as Projection


	//////////////////////////////////////////////////////////////////////////
	//[WARNING] Deprecated functions, planning to remove in next versions	
	////TODO: These function will be deprecated after using function object for fluid_cell
	////customized fluid cell access
	real Dia_Coef(const VectorDi& cell, const int i)			////return the alpha value on the incident face
	{
		int axis; VectorDi face; MacGrid<d>::Cell_Incident_Face(cell, i, axis, face); if ((*bc).Is_Psi_N(axis, face))return (real)0; else return (*alpha)(axis, face);
	}
	void Allocate_System(std::function<bool(const int)>& fluid_cell);
	void Build(std::function<bool(const int)>& fluid_cell);
	void Update_A(std::function<bool(const int)>& fluid_cell);
	void Correction(std::function<bool(const int)>& fluid_cell, FaceField<real, d>* correct_vel = nullptr);
};


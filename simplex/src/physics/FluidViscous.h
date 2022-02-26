//////////////////////////////////////////////////////////////////////////
// Viscous fluid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidViscous_h__
#define __FluidViscous_h__
#include <set>
#include "FluidEuler.h"
#include "Viscosity.h"

template<int d> class FluidViscous : public FluidEuler<d>
{Typedef_VectorDii(d);using Base=FluidEuler<d>;
public:
	using Base::mac_grid;using Base::velocity;using Base::g;using Base::type;using Base::alpha;
	using Base::Initialize_Fields;using Base::Enforce_Incompressibility;using Base::Enforce_Boundary_Conditions;
	using Base::Apply_Body_Forces;using Base::Update_Cell_Types;

	real nu=(real)1;
	BoundaryConditionMacGridViscosity<d> bc_vis=BoundaryConditionMacGridViscosity<d>(mac_grid);

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts,dx,domain_min);
		Initialize_Fields();
	}

	virtual void Advance(const real dt)
	{
		Apply_Viscosity_And_Advection(dt);
		Apply_Body_Forces(dt);
		Enforce_Incompressibility();
	}

	void Apply_Viscosity_And_Advection(const real dt)
	{
		FaceField<real,d> vel_copy=velocity;
		Viscosity<d> viscosity(mac_grid,velocity,alpha,type,bc_vis);
		viscosity.Solve(dt,nu);
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		Advection::Semi_Lagrangian(dt,mac_grid,velocity,mac_grid,vel_copy);
		velocity=vel_copy;
		Update_Cell_Types();
	}
};
#endif
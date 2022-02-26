//////////////////////////////////////////////////////////////////////////
// Incompressible fluid solver on an Eulerian grid, supporting free interface and surface tension
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidEulerImmersedBoundary_h__
#define __FluidEulerImmersedBoundary_h__
#include "FluidEuler.h"
#include "LevelSet.h"

template<int d> class FluidEulerImmersedBoundary : public FluidEuler<d>
{Typedef_VectorDii(d);using Base=FluidEuler<d>;
public:
	using Base::mac_grid;using Base::velocity;using Base::alpha;using Base::type;using Base::bc;using Base::g;using Base::use_body_force;
	FaceField<real,d> velocity_flow;
	FaceField<real,d> velocity_body;
	LevelSet<d> levelset;
	real epsilon;
	real one_over_epsilon;

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		Base::Initialize(cell_counts,dx,domain_min);
		Initialize_IB();
	}

	virtual void Initialize_IB()
	{
		velocity_flow.Resize(mac_grid.grid.cell_counts,(real)0);
		velocity_body.Resize(mac_grid.grid.cell_counts,(real)0);
		levelset.Initialize(mac_grid.grid);
		epsilon=mac_grid.grid.dx*(real)4;
		one_over_epsilon=(real)1/epsilon;
	}

	virtual void Advance(const real dt)
	{
		Advection(dt);
		Apply_Body_Forces(dt);
		Update_Alpha();
		Update_Meta_Velocity();
		Enforce_Incompressibility();
	}

	virtual void Advection(const real dt)
	{
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity_flow;
		Advection::Semi_Lagrangian(dt,ghost_grid,ghost_velocity,mac_grid,velocity_flow);
		Update_Cell_Types();
	}

	virtual void Apply_Body_Forces(const real dt)
	{
		if(!use_body_force)return;
		iterate_face(axis,iter,mac_grid){const VectorDi& face=iter.Coord();
			velocity_flow.face_fields[axis](face)+=g[axis]*dt;}
	}

	virtual void Update_Cell_Types()
	{
		iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			real phi=levelset.phi(cell);
			if(phi<-epsilon)type(cell)=(ushort)CellType::Solid;
			else if(phi>epsilon)type(cell)=(ushort)CellType::Air;
			else type(cell)=(ushort)CellType::IB;}	
	}

	void Update_Alpha()
	{
		iterate_face(axis,iter,mac_grid){const VectorDi& face=iter.Coord();
			const VectorD& pos=mac_grid.Face_Center(axis,face);
			real phi=levelset.Phi(pos);
			alpha(axis,face)=Kernel_Fluid(phi);}
	}

	virtual void Update_Meta_Velocity()
	{
		iterate_face(axis,iter,mac_grid){const VectorDi& face=iter.Coord();
			VectorD vf=velocity_flow(axis,face);VectorD vb=velocity_body(axis,face);
			velocity(axis,face)=alpha(axis,face)*vf+((real)1-alpha(axis,face))*vb;}
	}

	virtual void Enforce_Incompressibility()
	{
		Enforce_Boundary_Conditions();
		Projection<d> projection(mac_grid,velocity,alpha,type,bc);
		projection.use_alpha_for_correction=true;
		projection.Project();
	}

protected:
	real Kernel_Fluid(real phi) const
	{
		if(phi<-epsilon)return 0;
		else if(phi>epsilon)return (real)1;
		else return (real).5*((real)1+phi*one_over_epsilon+one_over_pi*sin(phi*one_over_epsilon*pi));
	}

	bool Is_Valid_Cell(const VectorDi& cell) const {return mac_grid.grid.Valid_Cell(cell);}
	bool Is_Fluid_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Fluid;}
	bool Is_IB_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::IB;}
};
#endif
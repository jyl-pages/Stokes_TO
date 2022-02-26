//////////////////////////////////////////////////////////////////////////
// advection-reflection fluid
// Copyright (c) (2021-), Fan Feng
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidAdvectionReflection_h__
#define __FluidAdvectionReflection_h__
template<int d> class FluidAdvectionReflection : public FluidEuler<d>
{
	Typedef_VectorDii(d); using Base = FluidEuler<d>;
public:
	using Base::velocity;
	using Base::mac_grid;
	using Base::projection;
	using Base::use_maccormack;
	using Base::FluidEuler;

	FaceField<real, d> velocity_temp;

	virtual void Advance(const real dt)
	{
		Advection(dt/(real)2);
		Apply_Body_Forces(dt / (real)2);
		velocity_temp = velocity; //u_tilde
		Enforce_Incompressibility();

		//u_hat^(1/2)=2u^(1/2)-u_tilde(1/2)
		velocity *= 2; 
		velocity -= velocity_temp;

		Advection(dt / (real)2);
		Enforce_Incompressibility();
	}
};
#endif
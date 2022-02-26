//////////////////////////////////////////////////////////////////////////
// Rigid body
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Rigid_h__
#define __Rigid_h__
#include "Common.h"
#include "Particles.h"
#include "AuxFunc.h"

template<int d> class Rigid{};

template<> class Rigid<2>
{static const int d=2;Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	real mass=(real)1;
	real one_over_mass=(real)1;
	real inertia_body=(real)1;
	real one_over_inertia_body=(real)1;
	real one_over_inertia=(real)1;

	VectorD position=VectorD::Zero();
	real orientation=(real)0;
	VectorD velocity=VectorD::Zero();
	real angular_momentum=(real)0;	////not used in 2D
	real angular_velocity=(real)0;
	VectorD force=VectorD::Zero();
	real torque=(real)0;

	void Advance(const real dt)
	{
		position+=dt*velocity;
		orientation+=dt*angular_velocity;
		velocity+=dt*force*one_over_mass;
		angular_velocity+=dt*torque*one_over_inertia;
	}

	void Clear_Force_And_Torque(){force=VectorD::Zero();torque=(real)0;}

	void Add_Force_In_World_Space(const VectorD& f,const VectorD& f_pos){VectorD f_r=f_pos-position;force+=f;torque+=AuxFunc::Cross(f_r,f)[0];}
	void Add_Force_In_Local_Space(const VectorD& f,const VectorD& f_pos){force+=Rotation2(orientation)*f;torque+=AuxFunc::Cross(f_pos,f)[0];}
	void Add_Body_Force_In_World_Space(const VectorD& f){force+=f;}
	void Add_Body_Force_In_Local_Space(const VectorD& f){force+=Rotation2(orientation)*f;}
	void Add_Torque_In_World_Space(const real tq){torque+=tq;}
	void Add_Torque_In_Local_Space(const real tq){torque+=tq;}
	void Add_Torque_In_World_Space(const VectorD& tq){}	////TEMP: for template compilation only, should never be called
	void Add_Torque_In_Local_Space(const VectorD& tq){}	////TEMP: for template compilation only, should never be called

	VectorD World_Space_Vector(const VectorD& object_space_vector) const {return Rotation2(orientation)*object_space_vector;}
	VectorD World_Space_Point(const VectorD& object_space_point) const {return position+Rotation2(orientation)*object_space_point;}
	VectorD Object_Space_Vector(const VectorD& world_space_vector) const {return Rotation2(-orientation)*world_space_vector;}
	VectorD Object_Space_Point(const VectorD& world_space_point) const {return Rotation2(-orientation)*(world_space_point-position);}

	MatrixD Inertia_Body() const {return MatrixD::Identity()*inertia_body;}	////TEMP: 2d inertia is a scalar, this is for template compilation only

	Transform2 Transform() const
	{Transform2 transform=Transform2::Identity();transform.pretranslate(position).rotate(orientation);return transform;}
	Transform2f Transformf() const
	{Transform2f transform=Transform2f::Identity();transform.pretranslate(position.cast<float>()).rotate((float)orientation);return transform;}

	//////////////////////////////////////////////////////////////////////////
	////particle initialization
	Array<VectorD> r0;

	void Initialize_From_Particles(const Particles<d>& particles,const int obj_idx_start,const int obj_idx_end)
	{
		mass=(real)0;
		position=VectorD::Zero();
		inertia_body=(real)0;
		int n=obj_idx_end-obj_idx_start;
		r0.resize(n);

		for(int i=obj_idx_start;i<obj_idx_end;i++){
			mass+=particles.M(i);
			position+=particles.M(i)*particles.X(i);}
		one_over_mass=(real)1/mass;
		position/=mass;

		for(int i=obj_idx_start,j=0;i<obj_idx_end;i++,j++){
			r0[j]=particles.X(i)-position;
			inertia_body+=particles.M(i)*pow(r0[j].norm(),2);}
		one_over_inertia_body=(real)1/inertia_body;
		one_over_inertia=one_over_inertia_body;
	}
	void Sync_Particles(Particles<d>& particles,const int obj_idx_start,const int obj_idx_end)
	{
		for(int i=obj_idx_start,j=0;i<obj_idx_end;i++,j++){
			VectorD r=Rotation2(orientation)*r0[j];
			particles.X(i)=position+r;
			particles.V(i)=velocity+AuxFunc::Cross(angular_velocity,r);}	
	}
};

template<> class Rigid<3>
{static const int d=3;Typedef_VectorDii(d);Typedef_MatrixD(d);Typedef_TransformD(d);
public:
	real mass=(real)1;
	real one_over_mass=(real)1;
	MatrixD inertia_body=MatrixD::Identity();
	MatrixD one_over_inertia_body=MatrixD::Identity();
	MatrixD one_over_inertia=MatrixD::Identity();

	VectorD position=VectorD::Zero();
	MatrixD orientation=MatrixD::Identity();
	
	VectorD velocity=VectorD::Zero();
	VectorD angular_momentum=VectorD::Zero();
	VectorD angular_velocity=VectorD::Zero();
	VectorD force=VectorD::Zero();
	VectorD torque=VectorD::Zero();

	void Advance(const real dt)
	{
		position+=dt*velocity;
		velocity+=dt*force*one_over_mass;

		orientation+=AuxFunc::Skew(angular_velocity)*dt*orientation;
		AngleAxis rot(orientation);orientation=rot.toRotationMatrix();	////normalizing the rotation matrix
		angular_momentum+=dt*torque;
		one_over_inertia=orientation*one_over_inertia_body*orientation.transpose();
		angular_velocity=one_over_inertia*angular_momentum;
	}

	void Clear_Force_And_Torque(){force=VectorD::Zero();torque=VectorD::Zero();}

	void Add_Force_In_World_Space(const VectorD& f,const VectorD& f_pos){VectorD f_r=f_pos-position;force+=f;torque+=f_r.cross(f);}
	void Add_Force_In_Local_Space(const VectorD& f,const VectorD& f_pos){force+=orientation*f;torque+=orientation*(f_pos.cross(f));}
	void Add_Body_Force_In_World_Space(const VectorD& f){force+=f;}
	void Add_Body_Force_In_Local_Space(const VectorD& f){force+=orientation*f;}
	void Add_Torque_In_World_Space(const VectorD& tq){torque+=tq;}
	void Add_Torque_In_Local_Space(const VectorD& tq){torque+=orientation*tq;}

	VectorD World_Space_Vector(const VectorD& object_space_vector) const {return orientation*object_space_vector;}
	VectorD World_Space_Point(const VectorD& object_space_point) const {return position+orientation*object_space_point;}

	VectorD Object_Space_Vector(const VectorD& world_space_vector) const {return orientation.transpose()*world_space_vector;}
	VectorD Object_Space_Point(const VectorD& world_space_point) const {return orientation.transpose()*(world_space_point-position);}

	void Set_Inertia(const MatrixD& _inertia_body)
	{inertia_body=_inertia_body;one_over_inertia_body=_inertia_body.inverse();one_over_inertia=orientation*one_over_inertia_body*orientation.transpose();}
	void Set_Inertia(const VectorD& _inertia_body)
	{MatrixD inertia_body_mat=MatrixD::Zero();for(int i=0;i<d;i++)inertia_body_mat.coeffRef(i,i)=_inertia_body[i];Set_Inertia(inertia_body_mat);}
	MatrixD Inertia_Body() const {return inertia_body;}

	void Set_Mass(const real _mass){mass=_mass;one_over_mass=(real)1./mass;}

	Transform3 Transform() const
	{Transform3 transform=Transform3::Identity();transform.pretranslate(position).rotate(orientation);return transform;}
	Transform3f Transformf() const 
	{Transform3f transform=Transform3f::Identity();transform.pretranslate(position.cast<float>()).rotate(orientation.cast<float>());return transform;}

	//////////////////////////////////////////////////////////////////////////
	////particle initialization
	Array<VectorD> r0;

	void Initialize_From_Particles(const Particles<d>& particles,const int obj_idx_start,const int obj_idx_end)
	{
		mass=(real)0;
		position=VectorD::Zero();
		MatrixD inertia_body=MatrixD::Zero();
		int n=obj_idx_end-obj_idx_start;
		r0.resize(n);

		for(int i=obj_idx_start;i<obj_idx_end;i++){
			mass+=particles.M(i);
			position+=particles.M(i)*particles.X(i);}
		one_over_mass=(real)1/mass;
		position/=(real)n;

		for(int i=obj_idx_start,j=0;i<obj_idx_end;i++,j++){
			r0[j]=particles.X(i)-position;
			inertia_body+=-particles.M(i)*AuxFunc::Skew(r0[j])*AuxFunc::Skew(r0[j]);}

		Set_Inertia(inertia_body);
		angular_velocity=one_over_inertia*angular_momentum;
	}
	
	void Sync_Particles(Particles<d>& particles,const int obj_idx_start,const int obj_idx_end)
	{
		for(int i=obj_idx_start,j=0;i<obj_idx_end;i++,j++){
			VectorD r=orientation*r0[j];
			particles.X(i)=position+r;
			particles.V(i)=velocity+angular_velocity.cross(r);}	
	}
};

#endif
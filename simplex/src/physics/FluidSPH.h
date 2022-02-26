//////////////////////////////////////////////////////////////////////////
// SPH Fluid
// Copyright (c) (2018-), Xiangxin Kong, Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidSPH_h__
#define __FluidSPH_h__
#include "Field.h"
#include "Particles.h"
#include "SparseFunc.h"
#include "SpatialHashing.h"
#include "MacGrid.h"
#include "BoundaryCondition.h"
#include "GeometryPrimitives.h"
#include "Kernels.h"

template<int d> class SPHParticles : public Particles<d>
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Particles<d>;
public:
    using Base::Size;using Base::Resize;
	using Base::X;using Base::V;using Base::F;using Base::M;
	using Base::C;	///curvature

	SPHParticles(){New_Attributes();}

	Declare_Attribute(real,ND,nrho,Rebind_ND);						////number density
	Declare_Attribute(real,Vol,vol,Rebind_Vol);						////volume
	Declare_Attribute(real,P,p,Rebind_P);							////pressure
	Declare_Attribute(VectorD,SF,sf,Rebind_SF);						////surface force
	Declare_Attribute(VectorD,SN,sn,Rebind_SN);						////surface normal

	Declare_Attribute_Inherent_Func(nrho,vol,p,sf,sn);	////does not need the attributes from the base class
};

template<int d,int bs=128> class FluidSPH
{Typedef_VectorDii(d);
public:
	SPHParticles<d> particles;
	Array<ArrayF<int,bs*4> > nbs;
	SpatialHashing<d,bs> spatial_hashing;
	std::shared_ptr<KernelSPH<d> > kernel;
	std::function<void(const real)> Collision;

	real length_scale;	////typical distance between two neighboring particles
	real den_0;			////rest density	
	real nden_0;		////rest number density
	real V_0;			////rest particle volume
	real mass_0;		////rest particle mass
	int avg_nb_num;		////averaged #particles
	real h;				////supporting radius			
	real kp;			////density-pressure coefficient
	real vis;			////viscosity
	real st;			////surface tension

	VectorD g=VectorD::Unit(1)*(real)-1.;
	bool use_body_force=false;
	bool use_fixed_dt=false;
	bool use_central_gravity=false;
	bool verbose=false;

public:
	void Initialize(const real _length_scale,const real _rho0=(real)1000.)
	{
		length_scale=_length_scale;
		den_0=_rho0;
		avg_nb_num=16*(int)pow(2,d-1);	////1d - 16; 2d - 32; 3d - 64

		if constexpr (d==1){
			V_0=length_scale;
			h=avg_nb_num*V_0;}
		else if constexpr (d==2){
			V_0=pi*pow((real).5*length_scale,2);
			h=sqrt((real)avg_nb_num*V_0/pi);}
		else if constexpr (d==3){
			V_0=(real)4/(real)3*pi*pow((real).5*length_scale,3);
			h=pow((real)3*(real)avg_nb_num*V_0/((real)4*pi),one_third);}

		mass_0=den_0*V_0;
		kp=(real)1e3;
		vis=den_0*(real)1e-2;
		st=(real)1;

		kernel=std::make_shared<KernelSPH<d> >(h);
		spatial_hashing.Initialize(h);

		////update nden_0
		nbs.resize(particles.Size());
		Update_Neighbors();
		int pn=particles.Size();
		nden_0=(real)0;
		for(int i=0;i<pn;i++){
			real nden=(real)0;
			int nb_n=nbs[i][0];
			for(int k=1;k<=nb_n;k++){
				int j=nbs[i][k];
				real dis_ij=(particles.X(i)-particles.X(j)).norm();
				nden+=kernel->W_Poly6(dis_ij);}
			particles.M(i)=den_0/nden;
			nden_0+=nden;}
		nden_0/=(real)pn;

		if(verbose)Print_Params();
	}

	void Update_Neighbors()
	{
		spatial_hashing.Update_Voxels(particles.XRef());
		int pn=particles.Size();
		for(int i=0;i<pn;i++){
			spatial_hashing.Find_Nbs(particles.X(i),particles.XRef(),h,nbs[i]);}	
	}

	virtual void Advance(const real dt,const real time=0)
	{	
		////update nbs
		Update_Neighbors();
		const int pn=particles.Size();

		////update number density and pressure
		for(int i=0;i<pn;i++){
			real nden=(real)0;
			int nb_n=nbs[i][0];
			for(int k=1;k<=nb_n;k++){
				int j=nbs[i][k];
				real dis_ij=(particles.X(i)-particles.X(j)).norm();
				nden+=kernel->W_Quintic(dis_ij);}
			particles.ND(i)=nden;
			particles.Vol(i)=(real)1/nden;
			particles.P(i)=kp*(nden/nden_0-(real)1.);}

		////update forces
		for(int i=0;i<pn;i++){
			real one_over_m=(real)1/particles.M(i);
			VectorD f=VectorD::Zero();
			int nb_n=nbs[i][0];
			for(int k=1;k<=nb_n;k++){
				int j=nbs[i][k];
				VectorD r_ij=particles.X(i)-particles.X(j);
				real r2=r_ij.squaredNorm();
				real one_over_r2=(r2==(real)0?(real)0:(real)1/r2);
				VectorD v_ij=particles.V(i)-particles.V(j);
				VectorD f_p=-(particles.P(i)*pow(particles.Vol(i),2)+particles.P(j)*pow(particles.Vol(j),2))*kernel->Grad_Spiky(r_ij);
				VectorD f_v=vis*particles.Vol(i)*particles.Vol(j)*v_ij*one_over_r2*r_ij.dot(kernel->Grad_Spiky(r_ij));
				f+=(f_p+f_v);}
			VectorD a_g=use_central_gravity?-particles.X(i).normalized():g;
			particles.F(i)=one_over_m*f+a_g;}

		////time integration
		for(int i=0;i<pn;i++){
			particles.V(i)+=particles.F(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}

		////collision
		if(Collision!=nullptr)Collision(dt);

		if(verbose)Print_Statistics();
	}

	real CFL() const
	{
		if(use_fixed_dt) return (real)1;
		else{
			real max_abs=(real)1e-5;
			for(int i=0;i<particles.Size();i++){
				real v_mag=particles.V(i).norm();
				if(v_mag>max_abs)max_abs=v_mag;}
			return length_scale/max_abs;}
	}

	void Print_Params()
	{
		std::cout<<"length_scale: "<<length_scale<<std::endl;
		std::cout<<"den_0: "<<den_0<<std::endl;
		std::cout<<"nden_0: "<<nden_0<<std::endl;
		std::cout<<"V_0: "<<V_0<<std::endl;
		std::cout<<"mass_0: "<<mass_0<<std::endl;
		std::cout<<"avg_nb_num: "<<avg_nb_num<<std::endl;
		std::cout<<"h: "<<h<<std::endl;
		std::cout<<"kp: "<<kp<<std::endl;
		std::cout<<"vis: "<<vis<<std::endl;
		std::cout<<"st: "<<st<<std::endl;
	}

	void Print_Statistics()
	{
		int pn=particles.Size();

		int avg_nb_num=0;
		int max_nb_num=-1;
		int min_nb_num=std::numeric_limits<int>::max();
		real avg_nden=(real)0;
		real max_nden=(real)-1;
		real min_nden=std::numeric_limits<real>::max();
		for(int i=0;i<pn;i++){
			avg_nb_num+=nbs[i][0];
			if(nbs[i][0]>max_nb_num)max_nb_num=nbs[i][0];
			if(nbs[i][0]<min_nb_num)min_nb_num=nbs[i][0];
			avg_nden+=particles.ND(i);
			if(particles.ND(i)>max_nden)max_nden=particles.ND(i);
			if(particles.ND(i)<min_nden)min_nden=particles.ND(i);}

		std::cout<<"nbs_num avg: "<<avg_nb_num/pn<<", max: "<<max_nb_num<<", min: "<<min_nb_num<<std::endl;
		std::cout<<"nden init: "<<nden_0<<", avg: "<<avg_nden/(real)pn<<", max: "<<max_nden<<", min: "<<min_nden<<std::endl;
	}

	//virtual void Update_Color_Field_And_Surface_Normal()
	//{
	//	for (int i=0; i<particles.Size(); i++)
	//	{
	//		real i_color=(real)0;
	//		VectorD i_sn=VectorD::Zero();
	//		for (int k=1; k<=nbs[i][0]; k++)
	//		{
	//			int j=nbs[i][k];
	//			VectorD r_ij=particles.X(i)-particles.X(j);
	//			real dist_ij=r_ij.norm();
	//			i_sn += particles.M(j)/particles.D(j)*kernel->Gradient_poly6(r_ij);
	//		}
	//		particles.C(i)=i_color;
	//		particles.SN(i)=i_sn;
	//	}
	//}

	//virtual void Update_Curvature_And_Surface_Force()
	//{
	//	for (int i=0; i<particles.Size(); i++){
	//		real i_curvature=(real)0;
	//		real surface_norm=particles.SN(i).norm();
	//		VectorD i_force=VectorD::Zero();
	//		for (int k=1; k<=nbs[i][0]; k++){
	//			int j=nbs[i][k];
	//			VectorD r_ij=particles.X(i)-particles.X(j);
	//			real dist_ij=r_ij.norm();				
	//			if(surface_norm != (real)0)
	//				i_curvature -= particles.M(j)/particles.D(j)*kernel->Wpoly6(dist_ij)/surface_norm;
	//		}
	//		particles.Curvature(i)=i_curvature;

	//		
	//		if (surface_norm > surface_threshold){
	//			real coef=surface_tension*i_curvature;
	//			i_force=coef *particles.SN(i);
	//			particles.F(i) -= i_force;
	//		}				
	//	}
	//}


	//void Enforce_Boundary_Conditions()
	//{
	//	for (int i=0; i<particles.Size(); i++)
	//	{
	//		for (int j=0; j<env_objects.size(); j++)
	//		{
	//			real phi=env_objects[j]->Phi(particles.X(i));
	//			if (phi<particles.R(i))
	//			{
	//				VectorD normal=env_objects[j]->Normal(particles.X(i));
	//				particles.F(i) += normal*kd*(particles.R(i)-phi)*particles.D(i);
	//			}
	//		}
	//	}
	//}
};
#endif
//////////////////////////////////////////////////////////////////////////
// Soft body mass spring
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyMassSpring_h__
#define __SoftBodyMassSpring_h__
#include "Hashtable.h"
#include "Particles.h"
#include "SparseFunc.h"

template<int d> class SoftBodyMassSpring
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	real ks_0=(real)1e5;
	real kd_0=(real)1e3;
protected:
	Particles<d>* particles_ptr=nullptr;
public:	
	Particles<d>& particles;
	bool own_data=false;

	Array<Vector2i> springs;
	Array<real> rest_length;

	bool use_varying_k=false;
	Array<real> ks;
	Array<real> kd;
	Array<VectorD> ext_force;

	Hashtable<int,VectorD> psi_D_values;	////displacements
	Hashtable<int,VectorD> psi_N_values;	////forces
	VectorD g=VectorD::Unit(1)*(real)-1.;
	bool use_body_force=true;
	
	bool use_implicit=true;
	SparseMatrixT K;
	VectorX u,b;

	SoftBodyMassSpring(Particles<d>* _p=nullptr)
		:particles_ptr(_p==nullptr?new Particles<d>():_p),particles(*particles_ptr){if(_p==nullptr)own_data=true;}
	~SoftBodyMassSpring(){if(own_data)delete particles_ptr;}

	virtual void Initialize()	////assuming springs and particles are initialized
	{
		int spring_size=(int)springs.size();
		rest_length.resize(spring_size);
		
		#pragma omp parallel for
		for(int i=0;i<spring_size;i++){const Vector2i& s=springs[i];
			rest_length[i]=(particles.X(s[0])-particles.X(s[1])).norm();}

		if(!use_varying_k){
			ks.resize(spring_size,ks_0);
			kd.resize(spring_size,kd_0);}
		
		ext_force.resize(particles.Size(),VectorD::Zero());

		if(use_implicit)Initialize_Implicit_K();
	}

	virtual void Advance(const real dt)
	{
		if(!use_implicit){
			real cfl=(real).005;int n=(int)(1./cfl);	////TOFIX: fix cfl
			for(int i=0;i<n;i++)Advance_Midpoint(dt*cfl);}
		else Advance_Implicit(dt);
		//Advance_Explicit_Euler(dt);
		//Advance_Midpoint(dt);
	}

	void Update_Forces(const Array<VectorD>& X,const Array<VectorD>& V,Array<VectorD>& F)
	{
		for(int i=0;i<particles.Size();i++){F[i]=ext_force[i];}
		for(int i=0;i<(int)springs.size();i++){
			VectorD f=Spring_Force(X,V,i);
			F[springs[i][0]]+=f;F[springs[i][1]]-=f;}
		if(use_body_force)for(int i=0;i<particles.Size();i++){
			F[i]+=particles.M(i)*g;}
		for(auto p:psi_N_values){int idx=p.first;const VectorD& f=p.second;F[idx]+=f;}
	}

	VectorD Spring_Force(const Array<VectorD>& X,const Array<VectorD>& V,const int s_idx)	////return f_ij=f_s+f+d
	{
		int i=springs[s_idx][0];int j=springs[s_idx][1];
		VectorD dir=X[j]-X[i];
		real length=dir.norm();
		dir/=length;
		VectorD rel_dir=V[j]-V[i];
		VectorD f_s=ks[s_idx]*(length-rest_length[s_idx])*dir;
		VectorD f_d=kd[s_idx]*rel_dir.dot(dir)*dir;
		return f_s+f_d;
	}

	static VectorD Spring_Force(const VectorD& Xi,const VectorD& Xj,const VectorD& Vi,const VectorD& Vj,const real rest_length,const real ks,const real kd)	////return f_ij=f_s+f+d
	{
		VectorD dir=Xj-Xi;
		real length=dir.norm();
		dir/=length;
		VectorD rel_dir=Vj-Vi;
		VectorD f_s=ks*(length-rest_length)*dir;
		VectorD f_d=kd*rel_dir.dot(dir)*dir;
		return f_s+f_d;
	}

	void Set_Psi_D(const int p,const VectorD v=VectorD::Zero()){psi_D_values[p]=v;}
	bool Is_Psi_D(const int p){return psi_D_values.find(p)!=psi_D_values.end();}
	void Set_Psi_N(const int p,const VectorD f){psi_N_values[p]=f;}
	bool Is_Psi_N(const int p){return psi_N_values.find(p)!=psi_N_values.end();}

	void Enforce_Boundary_Conditions(Array<VectorD>& V,Array<VectorD>& F)
	{
		for(auto p:psi_D_values){int idx=p.first;const VectorD& v=p.second;V[idx]=v;F[idx]=VectorD::Zero();}
	}

	virtual real CFL(){return (real).002;}

	void Advance_Explicit_Euler(const real dt)
	{
		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());
		for(int i=0;i<particles.Size();i++){
			particles.V(i)+=particles.F(i)/particles.M(i)*dt;
			particles.X(i)+=particles.V(i)*dt;}	
	}

	void Advance_Midpoint(const real dt)
	{

		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());	////update f_n and store it in particles.F
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());
		Array<VectorD> v_half(particles.Size(),VectorD::Zero());
		Array<VectorD> x_half(particles.Size(),VectorD::Zero());
		for(int i=0;i<particles.Size();i++){
			v_half[i]=particles.V(i)+dt*(real).5*particles.F(i);	////x_half=x_n+half_dt*v_n
			x_half[i]=particles.X(i)+dt*(real).5*particles.V(i);}	////v_half=v_n+half_dt*f_n

		Update_Forces(x_half,v_half,particles.FRef());			////update f_half and store it in particles.F

		Enforce_Boundary_Conditions(v_half,particles.FRef());
		for(int i=0;i<particles.Size();i++){
			particles.X(i)+=dt*v_half[i];						////x_{n+1}=x_n+dt*v_half
			particles.V(i)+=dt*particles.F(i)/particles.M(i);}	////v_{n+1}=v_n+dt*f_half
	}

	////Implicit time integration
	void Advance_Implicit(const real dt)
	{
		Update_Forces(particles.XRef(),particles.VRef(),particles.FRef());
		Enforce_Boundary_Conditions(particles.VRef(),particles.FRef());

		Update_Implicit_K(particles.XRef(),particles.VRef(),particles.FRef(),dt);
		for(int i=0;i<particles.Size();i++){
			VectorD v=particles.V(i);
			for(int j=0;j<d;j++)u[i*d+j]=v[j];}	////set initial guess to be the velocity from the last time step
		SparseSolver::Conjugate_Gradient(K,u,b);

		for(int i=0;i<particles.Size();i++){
			VectorD v;for(int j=0;j<d;j++)v[j]=u[i*d+j];
			particles.V(i)=v;
			particles.X(i)+=particles.V(i)*dt;}
	}

	void Initialize_Implicit_K()
	{
		int n=d*particles.Size();
		K.resize(n,n);u.resize(n);u.fill((real)0);b.resize(n);b.fill((real)0);
		Array<TripletT> elements;
		for(int s=0;s<(int)springs.size();s++){int i=springs[s][0];int j=springs[s][1];
			Add_Triplet_Helper(i,i,elements);
			Add_Triplet_Helper(i,j,elements);
			Add_Triplet_Helper(j,i,elements);
			Add_Triplet_Helper(j,j,elements);}
		K.setFromTriplets(elements.begin(),elements.end());
		K.makeCompressed();	
	}

	void Update_Implicit_K(const Array<VectorD>& X,const Array<VectorD>& V,const Array<VectorD>& F,const real dt)
	{
		////Clear K
		K.setZero();
		
		////Set K diagonal blocks and rhs
		for(int i=0;i<particles.Size();i++){
			for(int ii=0;ii<d;ii++){K.coeffRef(i*d+ii,i*d+ii)+=particles.M(i);}
			VectorD rhs=particles.M(i)*V[i]+dt*F[i];Set_Block(b,i,rhs);}

		////Set K blocks with spring and damping Jacobians
		MatrixD Ks,Kd;real dt_sq=pow(dt,2);
		for(int s=0;s<(int)springs.size();s++){int i=springs[s][0];int j=springs[s][1];
			Compute_Ks_Block(X,s,Ks);
			Ks*=-dt_sq;
			Add_Block_Helper(K,i,j,Ks);
			Compute_Kd_Block(X,s,Kd);
			Kd*=-dt;
			Add_Block_Helper(K,i,j,Kd);

			VectorD rhs_d_i=Kd*(V[i]-V[j]);
			if(!Is_Psi_D(i))Add_Block(b,i,rhs_d_i);
			if(!Is_Psi_D(j))Add_Block(b,j,-rhs_d_i);}
	}

	////Spring force derivative
	void Compute_Ks_Block(const Array<VectorD>& X,const int s,MatrixD& Ks)
	{
		int i=springs[s][0];int j=springs[s][1];
		VectorD x_ij=X[j]-X[i];
		real length=x_ij.norm();
		real length_0=rest_length[s];
		Ks=ks[s]*((length_0/length-(real)1)*MatrixD::Identity()-(length_0/pow(length,3))*x_ij*x_ij.transpose());
	}

	////Damping force derivative
	void Compute_Kd_Block(const Array<VectorD>& X,const int s,MatrixD& Kd)
	{
		int i=springs[s][0];int j=springs[s][1];
		VectorD n_ij=(X[j]-X[i]).normalized();
		Kd=-kd[s]*n_ij*n_ij.transpose();
	}

protected:
	void Add_Triplet_Helper(const int i,const int j,Array<TripletT>& elements)
	{
		for(int ii=0;ii<d;ii++)for(int jj=0;jj<d;jj++)elements.push_back(TripletT(i*d+ii,j*d+jj,(real)0));
	}

	void Add_Block_Helper(SparseMatrixT& K,const int i,const int j,const MatrixD& Ks)
	{
		SparseFunc::Add_Block<d,MatrixD>(K,i,i,Ks);
		SparseFunc::Add_Block<d,MatrixD>(K,j,j,Ks);
		if(!Is_Psi_D(i)&&!Is_Psi_D(j)){
			SparseFunc::Add_Block<d,MatrixD>(K,i,j,-Ks);
			SparseFunc::Add_Block<d,MatrixD>(K,j,i,-Ks);}
	}

	void Set_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]=bi[ii];}

	void Add_Block(VectorX& b,const int i,const VectorD& bi)
	{for(int ii=0;ii<d;ii++)b[i*d+ii]+=bi[ii];}
};

#endif
//////////////////////////////////////////////////////////////////////////
// Nonlinear Tetrahedral FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SoftBodyNonlinearFemTet.h"
#include "KrylovSolver.h"
#include "MeshFunc.h"
#include "AuxFunc.h"
#include "NonlinearFemFunc.h"
#include "Timer.h"

template<int d> void SoftBodyNonlinearFemTet<d>::Initialize(VolumetricMesh<d>& _mesh)
{
	mesh=std::make_shared<VolumetricMesh<d> >(particles.XPtr());
	particles.Resize((int)_mesh.Vertices().size());
	*mesh=_mesh;
	Add_Material((real)1,(real).45);

	int vtx_n=Vtx_Num();int dof_n=vtx_n*d;int ele_n=Ele_Num();
	K.resize(dof_n,dof_n);u.resize(dof_n);u.fill((real)0);f.resize(dof_n);f.fill((real)0);material_id.resize(ele_n,0);

	////initialize X0, Dm_inv, and M
	Set_Rest_Shape(particles.XRef());
	M.resize(dof_n);for(int i=0;i<M.diagonal().size();i++)M.diagonal()[i]=(real)0;
	Dm_inv.resize(ele_n,MatrixD::Zero());vol.resize(ele_n);
	for(int i=0;i<ele_n;i++){
		ArrayF<VectorD,d+1> vtx;for(int j=0;j<d+1;j++)vtx[j]=X0[E()[i][j]];
		real volume=(real)0;NonlinearFemFunc<d>::D_Inv_And_Vol(vtx,Dm_inv[i],volume);vol[i]=volume;
		for(int j=0;j<d+1;j++){int v=E()[i][j];
			for(int k=0;k<d;k++){M.diagonal()[v*d+k]+=density*volume/(real)(d+1);}}}

	AuxFunc::Fill(F(),VectorD::Zero());
	AuxFunc::Fill(V(),VectorD::Zero());

	////initialize implicit variables
	v.resize(dof_n);v.fill((real)0);
	a.resize(dof_n);a.fill((real)0);
	bf.resize(dof_n);bf.fill((real)0);
}

template<int d> void SoftBodyNonlinearFemTet<d>::Add_Material(real youngs,real poisson)
{materials.push_back(ElasticParam(youngs,poisson));}

template<int d> void SoftBodyNonlinearFemTet<d>::Set_Fixed(const int node)
{bc.psi_D_values[node]=VectorD::Zero();}

template<int d> void SoftBodyNonlinearFemTet<d>::Set_Displacement(const int node,const VectorD& dis)
{bc.psi_D_values[node]=dis;}

template<int d> void SoftBodyNonlinearFemTet<d>::Set_Force(const int node,const VectorD& force)
{bc.forces[node]=force;}

template<int d> void SoftBodyNonlinearFemTet<d>::Add_Force(const int node,const VectorD& force)
{bc.forces[node]+=force;}

template<int d> void SoftBodyNonlinearFemTet<d>::Clear_Force()
{bc.forces.clear();}

template<int d> void SoftBodyNonlinearFemTet<d>::Set_Rest_Shape(const Array<VectorD>& _X0)
{X0=_X0;}

template<int d> void SoftBodyNonlinearFemTet<d>::Advance(const real dt,const real time)
{
	if(use_explicit){
		if(!use_omp)Advance_Explicit(dt,time);
		else Advance_Explicit_Omp(dt,time);}
}

template<int d> void SoftBodyNonlinearFemTet<d>::Advance_Explicit(const real dt,const real time)
{
	Timer<real> timer;
	timer.Reset();

	AuxFunc::Fill(F(),VectorD::Zero());
	const int vtx_num=Vtx_Num();
	
	////body force, damping force, boundary
	if(use_body_force){
		for(int i=0;i<vtx_num;i++)F()[i]+=Mass(i)*g;}

	if(use_damping){
		for(int i=0;i<vtx_num;i++)F()[i]-=Mass(i)*damping*V()[i];}

	for(auto& iter:bc.forces){F()[iter.first]+=iter.second;}

	if(Collision!=nullptr)Collision(dt);
	if(Kinematic_Boundary_Condition!=nullptr)Kinematic_Boundary_Condition(dt,time);

	////elastic forces
	const int ele_num=Ele_Num();

	for(int i=0;i<ele_num;i++){
		ArrayF<VectorD,d+1> vtx;
		for(int j=0;j<d+1;j++)vtx[j]=X()[E()[i][j]];
		MatrixD ds;NonlinearFemFunc<d>::D(vtx,ds);
		MatrixD defo_grad=ds*Dm_inv[i];
		real mu=materials[material_id[i]].Mu<d>();
		real lambda=materials[material_id[i]].Lambda<d>();
		
		////material models, with 1st PK

		//////linear model
		//MatrixD grad_u=defo_grad-MatrixD::Identity();
		//MatrixD pk_stress=mu*(grad_u+grad_u.transpose())+lambda*grad_u.trace()*MatrixD::Identity();		

		//////co-rotated model
		MatrixD R,S;AuxFunc::Polar_Decomposition<d>(defo_grad,R,S);
		MatrixD pk_stress=(real)2*mu*(defo_grad-R)+lambda*(R.transpose()*defo_grad-MatrixD::Identity()).trace()*R;
		
		//////Neohookean model, has some problem
		//real J=defo_grad.determinant();MatrixD FinvT=defo_grad.inverse().transpose();
		//MatrixD pk_stress=mu*(defo_grad-mu*FinvT)+lambda*log(J)*FinvT;

		////////noninv co-rotated
		//real epsilon=(real)1e-3;
		//VectorD D;MatrixD U,V;AuxFunc::Svd_Rot_UV<d>(defo_grad,D,U,V);
		//for(int j=0;j<d;j++)if(D[j]<epsilon)D[j]=epsilon;
		//MatrixD defo_grad_noninv=MatrixD::Zero();
		//for(int j=0;j<d;j++)defo_grad_noninv.coeffRef(j,j)=D[j];
		//MatrixD pk_stress=(real)2*mu*(defo_grad_noninv-MatrixD::Identity())+lambda*(defo_grad_noninv-MatrixD::Identity()).trace()*MatrixD::Identity();
		//pk_stress=U*pk_stress*V.transpose();

		MatrixD ff=vol[i]*pk_stress*Dm_inv[i].transpose();

		//MatrixD cauchy_stress=(real)1/abs(defo_grad.determinant())*pk_stress*defo_grad.transpose();		////1st-PK -> Cauchy
		//MatrixD awn;Area_Weighted_Normals(X(),i,awn);
		//MatrixD ff=cauchy_stress*awn;

		{for(int j=0;j<d;j++){
			F()[E()[i][j]]-=ff.col(j);
			F()[E()[i][d]]+=ff.col(j);}}
	}

	////time integration
	for(int i=0;i<vtx_num;i++){
		V()[i]+=F()[i]/Mass(i)*dt;
		if(bc.psi_D_values.find(i)!=bc.psi_D_values.end())V()[i]=VectorD::Zero();
		X()[i]+=V()[i]*dt;}

	//timer.Elapse_And_Output("Sim");
}

//template<int d> void SoftBodyNonlinearFemTet<d>::Advance_Implicit(const real dt,const real time)
//{
//	VectorX old_u=u;VectorX old_v=v;
//	//if(use_corotation)Update_Matrices_Stiffness();
//	//if(use_heat){Advance_Heat(dt);Update_Matrices_Heat();}
//	//Update_RHS_With_Force_And_Heat(f);
//	SparseMatrixT A=K;
//	T coef=(T)4/(dt*dt)+(T)2*damping/dt;
//	for(int i=0;i<A.rows();i++)A.coeffRef(i,i)+=coef*M.diagonal()[i];
//	VectorX b=f-bf;
//	T coef2=(T)4/dt+damping;
//	for(int i=0;i<b.size();i++)b(i)+=M.diagonal()[i]*(coef*u[i]+coef2*v[i]+a[i]);
//	Modify_Matrix_And_RHS_With_Boundary_Conditions(A,b);
//	Solve_Linear_System(A,u,b);
//	for(int i=0;i<u.size();i++){v[i]=(u[i]-old_u[i])/dt;a[i]=(v[i]-old_v[i])/dt;}
//
//	////time integration
//	for(int i=0;i<Vtx_Num();i++){
//		V()[i]+=F()[i]/Mass(i)*dt;
//		if(bc.psi_D_values.find(i)!=bc.psi_D_values.end())V()[i]=VectorD::Zero();
//		X()[i]+=V()[i]*dt;}
//}

template<> void SoftBodyNonlinearFemTet<2>::Strain_To_Stress(const Matrix2& strain,const MatrixX& E,Matrix2& stress)
{
	Vector3 strain_vec={strain(0,0),strain(1,1),strain(0,1)};
	Vector3 stress_vec=E*strain_vec;
	stress<<stress_vec[0],stress_vec[2], stress_vec[2],stress_vec[1];
}
template<> void SoftBodyNonlinearFemTet<3>::Strain_To_Stress(const Matrix3& strain,const MatrixX& E,Matrix3& stress)
{/*TOIMPL*/}

template<> void SoftBodyNonlinearFemTet<2>::Area_Weighted_Normals(const Array<VectorD>& X,const int e,Matrix2& normals)
{
	const Vector3i& ele=E()[e];
	Vector2 v1=X[ele[2]]-X[ele[1]];
	Vector2 v2=X[ele[0]]-X[ele[2]];
	normals<<-v1[1],-v2[1], v1[0],v2[0];
}
template<> void SoftBodyNonlinearFemTet<3>::Area_Weighted_Normals(const Array<VectorD>& X,const int e,Matrix3& normals)
{/*TOIMPL*/}

//////////////////////////////////////////////////////////////////////////
////omp accelerated functions

inline void Add_Element_Force_To_Vertices(Array<Vector2>& F,const Vector3i& e,const Matrix2& ff)
{
	Vector2& f0=F[e[0]];
	Vector2& f1=F[e[1]];
	Vector2& f2=F[e[2]];
	Vector2 c0=ff.col(0);
	Vector2 c1=ff.col(1);
	Vector2 c2=c0+c1;
	#pragma omp atomic
	f0[0]-=c0[0];
	#pragma omp atomic
	f0[1]-=c0[1];
	#pragma omp atomic
	f1[0]-=c1[0];
	#pragma omp atomic
	f1[1]-=c1[1];
	#pragma omp atomic
	f2[0]+=c2[0];
	#pragma omp atomic
	f2[1]+=c2[1];
}

inline void Add_Element_Force_To_Vertices(Array<Vector3>& F,const Vector4i& e,const Matrix3& ff)
{ 
	Vector3& f0=F[e[0]];
	Vector3& f1=F[e[1]];
	Vector3& f2=F[e[2]];
	Vector3& f3=F[e[3]];
	Vector3 c0=ff.col(0);
	Vector3 c1=ff.col(1);
	Vector3 c2=ff.col(2);
	Vector3 c3=c0+c1+c2;
	#pragma omp atomic
	f0[0]-=c0[0];
	#pragma omp atomic
	f0[1]-=c0[1];
	#pragma omp atomic
	f0[2]-=c0[2];
	#pragma omp atomic
	f1[0]-=c1[0];
	#pragma omp atomic
	f1[1]-=c1[1];
	#pragma omp atomic
	f1[2]-=c1[2];
	#pragma omp atomic
	f2[0]-=c2[0];
	#pragma omp atomic
	f2[1]-=c2[1];
	#pragma omp atomic
	f2[2]-=c2[2];
	#pragma omp atomic
	f3[0]+=c3[0];
	#pragma omp atomic
	f3[1]+=c3[1];
	#pragma omp atomic
	f3[2]+=c3[2];
}

template<int d> void SoftBodyNonlinearFemTet<d>::Advance_Explicit_Omp(const real dt,const real time)
{
	Timer<real> timer;
	timer.Reset();

	const int vtx_num=Vtx_Num();
	const int ele_num=Ele_Num();

	#pragma omp parallel for
	for(int i=0;i<vtx_num;i++){
		auto& f=F()[i];
		f=VectorD::Zero();
		if(use_body_force)f+=Mass(i)*g;
		if(use_damping)f-=Mass(i)*damping*V()[i];}

	for(auto& iter:bc.forces){F()[iter.first]+=iter.second;}
	if(Collision!=nullptr)Collision(dt);
	if(Kinematic_Boundary_Condition!=nullptr)Kinematic_Boundary_Condition(dt,time);

	////elastic forces
	#pragma omp parallel for
	for(int i=0;i<ele_num;i++){
		ArrayF<VectorD,d+1> vtx;
		for(int j=0;j<d+1;j++)vtx[j]=X()[E()[i][j]];
		MatrixD ds;NonlinearFemFunc<d>::D(vtx,ds);
		MatrixD defo_grad=ds*Dm_inv[i];
		real mu=materials[material_id[i]].Mu<d>();
		real lambda=materials[material_id[i]].Lambda<d>();
		MatrixD R,S;AuxFunc::Polar_Decomposition<d>(defo_grad,R,S);
		MatrixD pk_stress=(real)2*mu*(defo_grad-R)+lambda*(R.transpose()*defo_grad-MatrixD::Identity()).trace()*R;
		MatrixD ff=vol[i]*pk_stress*Dm_inv[i].transpose();
		Add_Element_Force_To_Vertices(F(),E()[i],ff);}

	////time integration
	#pragma omp parallel for
	for(int i=0;i<vtx_num;i++){
		V()[i]+=F()[i]/Mass(i)*dt;
		if(bc.psi_D_values.find(i)!=bc.psi_D_values.end())V()[i]=VectorD::Zero();
		X()[i]+=V()[i]*dt;}

	//timer.Elapse_And_Output("Sim");
}

template class SoftBodyNonlinearFemTet<2>;
template class SoftBodyNonlinearFemTet<3>;

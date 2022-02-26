//////////////////////////////////////////////////////////////////////////
// Linear Tetrahedral FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SoftBodyLinearFemTet.h"
#include "KrylovSolver.h"
#include "MeshFunc.h"

template<int d> void SoftBodyLinearFemTet<d>::Initialize(VolumetricMesh<d>& _mesh)
{
	mesh=&_mesh;
	Add_Material((real)1,(real).45);

	int vtx_n=(int)mesh->Vertices().size();int dof_n=vtx_n*d;
	K.resize(dof_n,dof_n);u.resize(dof_n);u.fill((real)0);f.resize(dof_n);f.fill((real)0);
	material_id.resize(mesh->Elements().size(),0);
}

template<int d> void SoftBodyLinearFemTet<d>::Add_Material(real youngs,real poisson)
{materials.push_back(ElasticParam(youngs,poisson));}

template<int d> void SoftBodyLinearFemTet<d>::Allocate_K()
{
	Array<Vector2i> edges;MeshFunc::Get_Edges(*mesh,edges);
	Array<int> vertices;MeshFunc::Get_Vertices(*mesh,vertices);
	Array<TripletT> elements;

	for(int i=0;i<(int)vertices.size();i++){
		int r=vertices[i];int c=vertices[i];
		for(int rr=r*d;rr<(r+1)*d;rr++)for(int cc=c*d;cc<(c+1)*d;cc++){
			elements.push_back(TripletT(rr,cc,(real)0));}}
	for(int i=0;i<(int)edges.size();i++){
		const Vector2i& e=edges[i];int r=e[0];int c=e[1];
		for(int rr=r*d;rr<(r+1)*d;rr++)for(int cc=c*d;cc<(c+1)*d;cc++){
			elements.push_back(TripletT(rr,cc,(real)0));}
		r=e[1];c=e[0];
		for(int rr=r*d;rr<(r+1)*d;rr++)for(int cc=c*d;cc<(c+1)*d;cc++){
			elements.push_back(TripletT(rr,cc,(real)0));}}

	K.setFromTriplets(elements.begin(),elements.end());
	K.makeCompressed();	
}

template<int d> void SoftBodyLinearFemTet<d>::Update_K_And_f()
{
	//////Update K
	for(int i=0;i<(int)mesh->elements.size();i++){
		const VectorEi& e=mesh->elements[i];int mat=material_id[i];MatrixX Ke;
		ArrayF<VectorD,d+1> tet;for(int j=0;j<d+1;j++)tet[j]=mesh->Vertices()[e[j]];
		LinearFemFunc<d>::Tet_Stiffness_Matrix(materials[mat].youngs_modulus,materials[mat].poisson_ratio,tet,Ke);
		LinearFemFunc<d>::Add_Tet_Stiffness_Matrix(K,Ke,e);}

	////Update rhs
	f.fill((real)0);
	for(auto& b:bc.forces){int node=b.first;VectorD force=b.second;
		for(int axis=0;axis<d;axis++){int idx=node*d+axis;f[idx]+=force[axis];}}

	////Update bc
	for(auto& b:bc.psi_D_values){int node=b.first;VectorD dis=b.second;
		for(int axis=0;axis<d;axis++){int idx=node*d+axis;
			LinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(K,f,idx,dis[axis]);}}
}

template<int d> void SoftBodyLinearFemTet<d>::Solve()
{
	u.fill((real)0);
	KrylovSolver::ICPCG(K,u,f);
}

//template<int d> void SoftBodyLinearFemTet<d>::Advance(const real dt,const real time)
//{
//	damping=(T)2e-3;
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
//}

template<int d> void SoftBodyLinearFemTet<d>::Set_Fixed(const int node)
{
	bc.psi_D_values[node]=VectorD::Zero();
}

template<int d> void SoftBodyLinearFemTet<d>::Set_Displacement(const int node,const VectorD& dis)
{
	bc.psi_D_values[node]=dis;
}

template<int d> void SoftBodyLinearFemTet<d>::Set_Force(const int node,const VectorD& force)
{
	bc.forces[node]=force;
}

template<int d> void SoftBodyLinearFemTet<d>::Add_Force(const int node,const VectorD& force)
{
	bc.forces[node]+=force;
}

template<int d> void SoftBodyLinearFemTet<d>::Compute_Strain(Array<MatrixD>& strains) const
{
	Compute_Tet_Tensor_Helper(u,strains);
}

template<int d> void SoftBodyLinearFemTet<d>::Compute_Stress(Array<MatrixD>& stresses) const
{
	Array<MatrixX> E;E.resize(materials.size());
	for(int i=0;i<(int)materials.size();i++){
		LinearFemFunc<d>::Strain_Stress_Matrix_Linear(materials[i].youngs_modulus,materials[i].poisson_ratio,E[i]);}
    Compute_Tet_Tensor_Helper(u,stresses,&E);
}

template<int d> void SoftBodyLinearFemTet<d>::Compute_Von_Mises_Stress(Array<real>& von_mises,const Array<MatrixD>* stresses) const
{
	const Array<MatrixD>* str_ptr=nullptr;Array<MatrixD> str;
	if(stresses==nullptr){Compute_Stress(str);str_ptr=&str;}else{str_ptr=stresses;}
	von_mises.resize(str_ptr->size());
	for(int i=0;i<(int)str_ptr->size();i++){von_mises[i]=LinearFemFunc<d>::Von_Mises_Stress((*str_ptr)[i]);}
}

template<int d> void SoftBodyLinearFemTet<d>::Compute_Tet_Tensor_Helper(const VectorX& u,Array<MatrixD>& tensors,Array<MatrixX>* E) const
{
	tensors.resize((int)mesh->elements.size());
	for(int i=0;i<(int)mesh->elements.size();i++){
		const VectorEi& e=mesh->elements[i];int mat=material_id[i];MatrixX B;
		ArrayF<VectorD,d+1> tet;for(int j=0;j<d+1;j++)tet[j]=(mesh->Vertices()[e[j]]);
		LinearFemFunc<d>::Tet_Strain_Displacement_Matrix(tet,B);
		VectorX tet_u((d+1)*d);for(int ii=0;ii<d+1;ii++)for(int jj=0;jj<d;jj++)tet_u[ii*d+jj]=u[e[ii]*d+jj];
		VectorX tet_tensor;if(E){tet_tensor=(*E)[mat]*B*tet_u;}else tet_tensor=B*tet_u;
		LinearFemFunc<d>::Symmetric_Tensor_Vector_To_Matrix(tet_tensor,tensors[i]);}
}

template class SoftBodyLinearFemTet<2>;
template class SoftBodyLinearFemTet<3>;

//////////////////////////////////////////////////////////////////////////
// Linear Tetrahedral FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SoftBodyLinearFemHex.h"
#include "KrylovSolver.h"
#include "MeshFunc.h"

template<int d> void SoftBodyLinearFemHex<d>::Initialize(HexMesh<d>& _mesh)
{
	mesh=&_mesh;
	Add_Material((real)1,(real).45);

	int vtx_n=(int)mesh->Vertices().size();
	material_id.resize(vtx_n,0);
	////initialize dof
	int n=vtx_n*d;
	K.resize(n,n);u.resize(n);u.fill((real)0);f.resize(n);f.fill((real)0);
}

template<int d> void SoftBodyLinearFemHex<d>::Add_Material(real youngs,real poisson)
{materials.push_back(ElasticParam(youngs,poisson));}

template<int d> void SoftBodyLinearFemHex<d>::Allocate_K()
{
	Hashset<Vector2i> element_hashset;
	for(int i=0;i<(int)mesh->elements.size();i++){
		const VectorEi& e=mesh->elements[i];
		for(int j=0;j<ei;j++){int r=e[j];
			for(int k=0;k<ei;k++){int c=e[k];
				for(int rr=r*d;rr<(r+1)*d;rr++)for(int cc=c*d;cc<(c+1)*d;cc++){
					element_hashset.insert(Vector2i(rr,cc));}}}}
	Array<TripletT> elements;
	for(const auto& v:element_hashset){
		elements.push_back(TripletT(v[0],v[1],(real)0));}

	K.setFromTriplets(elements.begin(),elements.end());
	K.makeCompressed();	
}

template<int d> void SoftBodyLinearFemHex<d>::Update_K_And_f()
{
	//////Update K
	for(int i=0;i<(int)mesh->elements.size();i++){
		const VectorEi& e=mesh->elements[i];int mat=material_id[i];MatrixX Ke;
		ArrayF<VectorD,ei> hex;for(int j=0;j<ei;j++)hex[j]=mesh->Vertices()[e[j]];
		Array<int> idx(ei,0);for(int j=0;j<ei;j++)idx[j]=e[j];
		LinearFemFunc<d>::Hex_Stiffness_Matrix(materials[mat].youngs_modulus,materials[mat].poisson_ratio,hex,Ke);
		LinearFemFunc<d>::Add_Cell_Stiffness_Matrix(K,Ke,idx);}

	////Update rhs
	f.fill((real)0);
	for(auto& b:bc.forces){int node=b.first;VectorD force=b.second;
		for(int axis=0;axis<d;axis++){int idx=node*d+axis;f[idx]+=force[axis];}}

	////Update bc
	for(auto& b:bc.psi_D_values){int node=b.first;VectorD dis=b.second;
		for(int axis=0;axis<d;axis++){int idx=node*d+axis;
			LinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(K,f,idx,dis[axis]);}}
}

template<int d> void SoftBodyLinearFemHex<d>::Solve()
{
	u.fill((real)0);
	KrylovSolver::ICPCG(K,u,f);
}

template<int d> void SoftBodyLinearFemHex<d>::Set_Fixed(const int node)
{
	bc.psi_D_values[node]=VectorD::Zero();
}

template<int d> void SoftBodyLinearFemHex<d>::Set_Displacement(const int node,const VectorD& dis)
{
	bc.psi_D_values[node]=dis;
}

template<int d> void SoftBodyLinearFemHex<d>::Set_Force(const int node,const VectorD& force)
{
	bc.forces[node]=force;
}

template<int d> void SoftBodyLinearFemHex<d>::Add_Force(const int node,const VectorD& force)
{
	bc.forces[node]+=force;
}

template class SoftBodyLinearFemHex<2>;
template class SoftBodyLinearFemHex<3>;

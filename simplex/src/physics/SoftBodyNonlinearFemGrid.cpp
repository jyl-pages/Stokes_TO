//////////////////////////////////////////////////////////////////////////
// Non-Linear Grid FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SoftBodyNonlinearFemGrid.h"
#include "KrylovSolver.h"
#include "ColorGrid.h"
#include "AuxFunc.h"
#include "LinearFemFunc.h"
#include "NonlinearFemFunc.h"

using namespace AuxFunc;

template<int d> void SoftBodyNonlinearFemGrid<d>::Initialize(const Grid<d>& _grid)
{
	grid=_grid;
	materials.clear();
	Add_Material((real)1,(real).45);
	material_id.Resize(grid.cell_counts,0);

	int n=grid.node_counts.prod()*d;
	K.resize(n,n);u.resize(n);u.fill((real)0);f.resize(n);f.fill((real)0);
	f_elas.resize(n);f_elas.fill((real)0);e_elas.resize(n);
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Add_Material(real youngs,real poisson)
{
	ElasticParam mat(youngs,poisson);
	materials.push_back(mat);
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Solve_Nonlinear()
{
	VectorX f_ext(f.size());
	Set_RHS_Force(f_ext);
	real r0=AuxFunc::Linf(f_ext);

	Allocate_K();
	int i=0;

	do{
		i++;
		Update_K();
		f=f_ext-f_elas;
		Set_K_And_RHS_Psi_D(f);

		real r=AuxFunc::Linf(f); 
		std::cout<<"r= "<<r<<std::endl;
		if(r<abs_tol||r<rel_tol*r0){
			std::cout<<"Newton converges in "<<i<<" iterations"<<std::endl;break;}
		
		VectorX delta_u(u.size());delta_u.fill((real)0);
		Solve_Linear(delta_u,f);
		real alpha=(real)1.;
		u+=alpha*delta_u;
	}while(i<max_iter);
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Allocate_K()
{
    Array<TripletT> elements;
    iterate_node(iter,grid){const VectorDi& node=iter.Coord();int r=grid.Node_Index(node);
        for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++){VectorDi nb=Grid<d>::Nb_R(node,i);
			if(grid.Valid_Node(nb)){int c=Node_Index_In_K(nb);
				for(int rr=r*d;rr<(r+1)*d;rr++)for(int cc=c*d;cc<(c+1)*d;cc++){elements.push_back(TripletT(rr,cc,(real)0));}}}}
	K.setFromTriplets(elements.begin(),elements.end());
	K.makeCompressed();	

	ColorGrid<d>::Color(grid.cell_counts,colored_cell_ptr,colored_cell_indices);
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Update_K()
{
	f_elas.fill((real)0);

	////Update K, f_elas, and e_elas, using the current u
	int color_n=ColorGrid<d>::Number_Of_Colors();
	for(int c=0;c<color_n;c++){
		#pragma omp parallel for
		for(int i=colored_cell_ptr[c];i<colored_cell_ptr[c+1];i++){
			VectorDi cell=grid.Cell_Coord(colored_cell_indices[i]);int mat_id=material_id(cell);
			ArrayF2P<VectorDi,d> cell_nodes;grid.Cell_Incident_Nodes(cell,cell_nodes);
			Array<int> cell_node_indices(Grid<d>::Number_Of_Cell_Incident_Nodes(),0);
			for(int j=0;j<cell_node_indices.size();j++)cell_node_indices[j]=Node_Index_In_K(cell_nodes[j]);

			{real mu=materials[mat_id].Mu<d>();real lambda=materials[mat_id].Lambda<d>();MatrixX Ke_nonlinear;
			VectorX cell_f_nonlinear((int)pow(2,d)*d);cell_f_nonlinear.fill((real)0);real cell_e_nonlinear=(real)0;
			VectorX cell_u((int)pow(2,d)*d);Compute_Cell_Displacement(u,cell,cell_u);
			NonlinearFemFunc<d>::Cell_Stiffness_Matrix_And_f_Nonlinear(mu,lambda,grid.dx,cell_u,Ke_nonlinear,&cell_f_nonlinear,&cell_e_nonlinear);
			LinearFemFunc<d>::Add_Cell_Stiffness_Matrix(K,variable_coef?Ke_nonlinear*(*variable_coef)(cell):Ke_nonlinear,cell_node_indices);
			LinearFemFunc<d>::Add_Cell_Vector(f_elas,cell_f_nonlinear,cell_node_indices);
			e_elas[colored_cell_indices[i]]=cell_e_nonlinear;}}}
}

////Update rhs with force bc
template<int d> void SoftBodyNonlinearFemGrid<d>::Set_RHS_Force(VectorX& _f)
{
	_f.fill((real)0);
	for(auto& b:bc.forces){VectorDi node=b.first;VectorD force=b.second;
		for(int axis=0;axis<d;axis++){int idx=Node_Index_In_K(node)*d+axis;_f[idx]+=force[axis];}}
}

////Update K and rhs with displacement bc
template<int d> void SoftBodyNonlinearFemGrid<d>::Set_K_And_RHS_Psi_D(VectorX& _f)
{
	for(auto& b:bc.psi_D_values){VectorDi node=b.first;VectorD dis=b.second;
		for(int axis=0;axis<d;axis++){int idx=Node_Index_In_K(node)*d+axis;
			LinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(K,_f,idx,dis[axis]);}}
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Solve_Linear(VectorX& _u,const VectorX& _f)
{
	_u.fill((real)0);

	multigrid_params.use_auto_calculated_levels=true;
	multigrid_params.dof_on_cell=false;
	multigrid_params.block_size=d;
	multigrid_params.use_gpu=true;
	
	#ifdef USE_CUDA
	if(multigrid_params.use_gpu){
		//GeometricMultiGrid::GMGPCG_GPU<d>(K,_u,_f,grid.node_counts,multigrid_params);
		gmg_solver_gpu.update_A_levels=true;
		gmg_solver_gpu.Initialize(K,grid.node_counts,multigrid_params,&material_id);
		gmg_solver_gpu.Solve(_u,_f);
	}
	else{ 
		gmg_solver_cpu.update_A_levels=true;
		gmg_solver_cpu.Initialize(K,grid.node_counts,multigrid_params,&material_id);
		gmg_solver_cpu.Solve(_u,_f);}
	#else
		gmg_solver_cpu.Initialize(K,grid.node_counts,multigrid_params,&material_id);
		gmg_solver_cpu.Solve(_u,_f);
	#endif
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Set_Fixed(const VectorDi& node)
{
	bc.psi_D_values[node]=VectorD::Zero();
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Set_Displacement(const VectorDi& node,const VectorD& dis)
{
	bc.psi_D_values[node]=dis;
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Set_Force(const VectorDi& node,const VectorD& force)
{
	bc.forces[node]=force;
}

////helper functions
template<int d> int SoftBodyNonlinearFemGrid<d>::Node_Index_In_K(const VectorDi& node) const
{return grid.Node_Index(node);}

template<int d> void SoftBodyNonlinearFemGrid<d>::Compute_Cell_Displacement(const VectorX& u,const VectorDi& cell,VectorX& cell_u) const
{
	int number_of_cell_nodes=Grid<d>::Number_Of_Cell_Incident_Nodes();
	cell_u.resize(number_of_cell_nodes*d);
	for(int i=0;i<number_of_cell_nodes;i++){
		VectorDi nb_node=Grid<d>::Cell_Incident_Node(cell,i);int nb_node_mtx_idx=Node_Index_In_K(nb_node);
		for(int j=0;j<d;j++)cell_u(i*d+j)=u(nb_node_mtx_idx*d+j);}
}

template<int d> void SoftBodyNonlinearFemGrid<d>::Compute_Cell_Displacement(const VectorX& u,const Array<int>& cell_node_matrix_indices,VectorX& cell_u) const
{
	int number_of_cell_nodes=Grid<d>::Number_Of_Cell_Incident_Nodes();
	cell_u.resize(number_of_cell_nodes*d);
	for(int i=0;i<number_of_cell_nodes;i++){
		int nb_node_mtx_idx=cell_node_matrix_indices[i];
		for(int j=0;j<d;j++)cell_u(i*d+j)=u(nb_node_mtx_idx*d+j);}
}

template class SoftBodyNonlinearFemGrid<2>;
template class SoftBodyNonlinearFemGrid<3>;

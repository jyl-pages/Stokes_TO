//////////////////////////////////////////////////////////////////////////
// Linear Grid FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "SoftBodyLinearFemGrid.h"
#include "ColorGrid.h"

//template<int d> SoftBodyLinearFemGrid<d>::SoftBodyLinearFemGrid():colored_cell_ptr(Pow(2,d)+1)
//{for(int i=0;i<colored_cell_ptr.size();i++)colored_cell_ptr[i]=i;}

template<int d> void SoftBodyLinearFemGrid<d>::Initialize(const Grid<d> _grid)
{
	colored_cell_ptr.resize(Pow(2, d) + 1);
	for (int i = 0; i < colored_cell_ptr.size(); i++)colored_cell_ptr[i] = i;

	grid=_grid;
	Add_Material((real)1,(real).45);
	material_id.Resize(grid.cell_counts,0);

	int n=grid.node_counts.prod()*d;
	K.resize(n,n);u.resize(n);u.fill((real)0);f.resize(n);f.fill((real)0);
}

template<int d> void SoftBodyLinearFemGrid<d>::Add_Material(real youngs,real poisson)
{
	MatrixX E0;LinearFemFunc<d>::Strain_Stress_Matrix_Linear(youngs,poisson,E0);E.push_back(E0);
	MatrixX Ke0;LinearFemFunc<d>::Cell_Stiffness_Matrix(youngs,poisson,grid.dx,Ke0);Ke.push_back(Ke0);
}

template<int d> void SoftBodyLinearFemGrid<d>::Allocate_K()
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

template<int d> void SoftBodyLinearFemGrid<d>::Update_K_And_f()
{
	////Update K
	int color_n=ColorGrid<d>::Number_Of_Colors();
	for(int c=0;c<color_n;c++){
		#pragma omp parallel for
		for(int i=colored_cell_ptr[c];i<colored_cell_ptr[c+1];i++){
			VectorDi cell=grid.Cell_Coord(colored_cell_indices[i]);int mat_id=material_id(cell);
			ArrayF2P<VectorDi,d> cell_nodes;grid.Cell_Incident_Nodes(cell,cell_nodes);
			Array<int> cell_node_indices(Grid<d>::Number_Of_Cell_Incident_Nodes(),0);
			for(int j=0;j<cell_node_indices.size();j++)cell_node_indices[j]=Node_Index_In_K(cell_nodes[j]);
			LinearFemFunc<d>::Add_Cell_Stiffness_Matrix(K,variable_coef?Ke[mat_id]*(*variable_coef)(cell):Ke[mat_id],cell_node_indices);}}

	////Update rhs
	f.fill((real)0);
	for(auto& b:bc.forces){VectorDi node=b.first;VectorD force=b.second;
		for(int axis=0;axis<d;axis++){int idx=Node_Index_In_K(node)*d+axis;f[idx]+=force[axis];}}

	////Update bc
	for(auto& b:bc.psi_D_values){VectorDi node=b.first;VectorD dis=b.second;
		for(int axis=0;axis<d;axis++){int idx=Node_Index_In_K(node)*d+axis;
			LinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(K,f,idx,dis[axis]);}}
}

template<int d> void SoftBodyLinearFemGrid<d>::Solve()
{
	u.fill((real)0);

	multigrid_params.use_auto_calculated_levels=true;
	multigrid_params.dof_on_cell=false;
	multigrid_params.block_size=d;
	multigrid_params.use_gpu=true;
	multigrid_params.init_hier_on_gpu=true;	////calculate hier on CPU to avoid GPU memory crash

	#ifdef USE_CUDA
	if(multigrid_params.use_gpu){
		//GeometricMultiGrid::GMGPCG_GPU<d>(K,u,f,grid.node_counts,multigrid_params);
		gmg_solver_gpu.update_A_levels=true;

		gmg_solver_gpu.Initialize(K,grid.node_counts,multigrid_params,&material_id);
		gmg_solver_gpu.Solve(u,f);
	}
	else{ 
		gmg_solver_cpu.update_A_levels=true;
		gmg_solver_cpu.Initialize(K,grid.node_counts,multigrid_params,&material_id);
		gmg_solver_cpu.Solve(u,f);}
	#else
		gmg_solver_cpu.Initialize(K,grid.node_counts,multigrid_params,&material_id);
		gmg_solver_cpu.Solve(u,f);
	#endif
}

template<int d> void SoftBodyLinearFemGrid<d>::Set_Fixed(const VectorDi& node)
{
	bc.psi_D_values[node]=VectorD::Zero();
}

template<int d> void SoftBodyLinearFemGrid<d>::Set_Displacement(const VectorDi& node,const VectorD& dis)
{
	bc.psi_D_values[node]=dis;
}

template<int d> void SoftBodyLinearFemGrid<d>::Set_Force(const VectorDi& node,const VectorD& force)
{
	bc.forces[node]=force;
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Strain(Field<MatrixD,d>& strains) const
{	
    Array<MatrixX> B(1);VectorD point=VectorD::Zero();
    LinearFemFunc<d>::Cell_Strain_Displacement_Matrix(point,grid.dx,B[0]);
	Compute_Cell_Tensor_Helper(B,u,strains);
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Stress(Field<MatrixD,d>& stresses) const
{
    MatrixX B;VectorD point=VectorD::Zero();
    LinearFemFunc<d>::Cell_Strain_Displacement_Matrix(point,grid.dx,B);

	Array<MatrixX> S;S.resize((int)E.size());
	for(int i=0;i<(int)E.size();i++){S[i]=E[i]*B;}
    Compute_Cell_Tensor_Helper(S,u,stresses);
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Von_Mises_Stress(Field<real,d>& von_mises,const Field<MatrixD,d>* stresses) const
{
	von_mises.Resize(grid.cell_counts);
	const Field<MatrixD,d>* str_ptr=nullptr;Field<MatrixD,d> str;
	if(stresses==nullptr){Compute_Stress(str);str_ptr=&str;}else{str_ptr=stresses;}
	for(auto i=0;i<(int)str_ptr->array.size();i++){
		von_mises.array[i]=LinearFemFunc<d>::Von_Mises_Stress(str_ptr->array[i]);}
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Strain(VectorX& strain,const VectorDi& cell) const
{
	VectorX cell_u;Compute_Cell_Displacement(u,cell,cell_u);
	MatrixX B;VectorD point=VectorD::Zero();
	LinearFemFunc<d>::Cell_Strain_Displacement_Matrix(point,grid.dx,B);
	strain=B*cell_u;
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Strain(MatrixD& strain,const VectorDi& cell) const
{
	VectorX strain_vec;Compute_Strain(strain_vec,cell);
	LinearFemFunc<d>::Symmetric_Tensor_Vector_To_Matrix(strain_vec,strain);
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Stress(/*rst*/VectorX& stress,const VectorDi& cell,const MatrixX& E) const
{
	VectorX strain_vec;Compute_Strain(strain_vec,cell);stress=E*strain_vec;
	if(variable_coef)stress*=(*variable_coef)(cell);
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Stress(/*rst*/MatrixD& stress,const VectorDi& cell,const MatrixX& E) const
{
	VectorX stress_vec;Compute_Stress(stress_vec,cell,E);
	LinearFemFunc<d>::Symmetric_Tensor_Vector_To_Matrix(stress_vec,stress);
	if(variable_coef)stress*=(*variable_coef)(cell);
}

////helper functions
template<int d> int SoftBodyLinearFemGrid<d>::Node_Index_In_K(const VectorDi& node) const
{return grid.Node_Index(node);}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Cell_Tensor_Helper(const Array<MatrixX>& B,const VectorX& u,Field<MatrixD,d>& tensors) const
{
	tensors.Resize(grid.cell_counts);
	iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();int mat_id=material_id(cell);
		if(mat_id!=-1){
			VectorX cell_u;Compute_Cell_Displacement(u,cell,cell_u);VectorX cell_tensor=B[mat_id]*cell_u;
			LinearFemFunc<d>::Symmetric_Tensor_Vector_To_Matrix(cell_tensor,tensors(cell));}
		else{tensors(cell)=MatrixD::Zero();}}
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Cell_Displacement(const VectorX& u,const VectorDi& cell,VectorX& cell_u) const
{
	int number_of_cell_nodes=Grid<d>::Number_Of_Cell_Incident_Nodes();
	cell_u.resize(number_of_cell_nodes*d);
	for(int i=0;i<number_of_cell_nodes;i++){
		VectorDi nb_node=Grid<d>::Cell_Incident_Node(cell,i);int nb_node_mtx_idx=Node_Index_In_K(nb_node);
		for(int j=0;j<d;j++)cell_u(i*d+j)=u(nb_node_mtx_idx*d+j);}
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Cell_Displacement(const VectorX& u,const Array<int>& cell_node_matrix_indices,VectorX& cell_u) const
{
	int number_of_cell_nodes=Grid<d>::Number_Of_Cell_Incident_Nodes();
	cell_u.resize(number_of_cell_nodes*d);
	for(int i=0;i<number_of_cell_nodes;i++){
		int nb_node_mtx_idx=cell_node_matrix_indices[i];
		for(int j=0;j<d;j++)cell_u(i*d+j)=u(nb_node_mtx_idx*d+j);}
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Elastic_Compliance(Field<real,d>& energy)
{
	energy.Resize(grid.cell_counts,(real)0);
	//iterate_cell(iter,grid){VectorDi cell=iter.Coord();int mat_id=material_id(cell);if(mat_id==-1)continue;
	int cell_n=grid.Number_Of_Cells();
	#pragma omp parallel for
	for(int i=0;i<cell_n;i++){VectorDi cell=grid.Cell_Coord(i);
		int mat_id=material_id(cell);if(mat_id==-1)continue;
		VectorX cell_u;Compute_Cell_Displacement(u,cell,cell_u);
			const MatrixX& K0=Ke[mat_id];
			energy(cell)=cell_u.dot(K0*cell_u);
			//if(variable_coef){energy(cell)*=(*variable_coef)(cell);}
			}
}

template<int d> void SoftBodyLinearFemGrid<d>::Compute_Elastic_Energy(Field<real, d>& energy) const
{
	energy.Resize(grid.cell_counts, (real)0);
	//iterate_cell(iter,grid){VectorDi cell=iter.Coord();int mat_id=material_id(cell);if(mat_id==-1)continue;
	int cell_n = grid.Number_Of_Cells();
#pragma omp parallel for
	for (int i = 0; i < cell_n; i++) {
		VectorDi cell = grid.Cell_Coord(i);
		int mat_id = material_id(cell); if (mat_id == -1)continue;
		VectorX cell_u; Compute_Cell_Displacement(u, cell, cell_u);
		const MatrixX& K0 = Ke[mat_id];
		energy(cell) = cell_u.dot(K0*cell_u);
		if (variable_coef) { energy(cell) *= (*variable_coef)(cell); }
	}
}

template class SoftBodyLinearFemGrid<2>;
template class SoftBodyLinearFemGrid<3>;

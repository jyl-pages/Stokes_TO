//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid
// Copyright (F) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include <numeric>
#include <functional>
#include "AuxFunc.h"
#include "SparseFunc.h"
#include "Grid.h"
#include "Field.h"
#include "Timer.h"

#include "GeometricMultiGrid.h"
//#ifdef USE_CUDA
//#include "MultiGridCuda.h"
//#endif
#include "ColorGrid.h"

namespace GeometricMultiGrid
{
using namespace MultiGrid;

//////////////////////////////////////////////////////////////////////////
////Build R for regular domains
//////////////////////////////////////////////////////////////////////////

template<int d> Vector<int,d> Coarsened(const Vector<int,d>& counts,const int coarsen_factor)
{Vector<int,d> coarsened;for(int i=0;i<d;i++){coarsened[i]=counts[i]/coarsen_factor;if(coarsened[i]==0)coarsened[i]=1;}return coarsened;}

////Four-point interpolation: For multigrid matrices discretized on grid cells
inline real Interpolation_Coefficient_Four_Points(const Vector1i& offset)
{ArrayF<real,4> coef={(real).125,(real).375,(real).375,(real).125};return coef[offset[0]];}
inline real Interpolation_Coefficient_Four_Points(const Vector2i& offset)
{ArrayF<real,4> coef={(real).125,(real).375,(real).375,(real).125};return coef[offset[0]]*coef[offset[1]];}
inline real Interpolation_Coefficient_Four_Points(const Vector3i& offset)
{ArrayF<real,4> coef={(real).125,(real).375,(real).375,(real).125};return coef[offset[0]]*coef[offset[1]]*coef[offset[2]];}
template<int d,class T_SPARSE> void Build_Restriction_Matrix_Four_Points(Vector<int,d>& cell_counts_fine,const int number_of_components,
	const int coarsen_factor,const bool use_psi_P,/*rst*/T_SPARSE& R);

inline int Number_Of_Nb_2RV_1D(const int coord,const int count){int n=4;if(coord-1<0)n--;if(coord+1>=count)n-=2;else if(coord+2>=count)n--;return n;}
inline int Number_Of_Nb_2RV(const Vector1i& coord,const Vector1i& counts)
{return Number_Of_Nb_2RV_1D(coord[0],counts[0]);}
inline int Number_Of_Nb_2RV(const Vector2i& coord,const Vector2i& counts)
{return Number_Of_Nb_2RV_1D(coord[0],counts[0])*Number_Of_Nb_2RV_1D(coord[1],counts[1]);}
inline int Number_Of_Nb_2RV(const Vector3i& coord,const Vector3i& counts)
{return Number_Of_Nb_2RV_1D(coord[0],counts[0])*Number_Of_Nb_2RV_1D(coord[1],counts[1])*Number_Of_Nb_2RV_1D(coord[2],counts[2]);}

////Three-point interpolation: For multigrid matrices discretized on grid nodes
inline real Interpolation_Coefficient_Three_Points(const Vector1i& offset)
{return (real).5*pow((real).5,abs(offset[0]));}
inline real Interpolation_Coefficient_Three_Points(const Vector2i& offset)
{return (real).25*pow((real).5,abs(offset[0])+abs(offset[1]));}
inline real Interpolation_Coefficient_Three_Points(const Vector3i& offset)
{return (real).125*pow((real).5,abs(offset[0])+abs(offset[1])+abs(offset[2]));}
template<int d,class T_SPARSE> void Build_Restriction_Matrix_Three_Points(Vector<int,d>& cell_counts_fine,const int number_of_components,
	const int coarsen_factor,const bool use_psi_P,/*rst*/T_SPARSE& R);

template<int d,class T_SPARSE>
void Build_Restriction_Matrix_Four_Points(
	Vector<int,d>& cell_counts_fine,const int number_of_components,
	const int coarsen_factor,const bool use_psi_P,/*rst*/T_SPARSE& R)
{	Typedef_VectorDii(d);
	Grid<d> grid_fine(cell_counts_fine,(real)1);
	Grid<d> grid_coarse=grid_fine.Coarsened(coarsen_factor);
	VectorDi coarsen_factor_vec=grid_fine.cell_counts.cwiseQuotient(grid_coarse.cell_counts);
	VectorDi zero=VectorDi::Zero();

	int m=grid_coarse.Number_Of_Cells();
	int n=grid_fine.Number_Of_Cells();
	int rows=(m*number_of_components);
	int cols=(n*number_of_components);
	R.resize(rows,cols);
	int* ptr=R.outerIndexPtr();
	ptr[0]=0;

	////Allocate R
	#pragma omp parallel for
	for(int i=0;i<m;i++){
		VectorDi cell=grid_coarse.Cell_Coord(i);
		VectorDi cell_fine=coarsen_factor_vec.cwiseProduct(cell);
		int n=use_psi_P?Pow(4,d):Number_Of_Nb_2RV(cell_fine,grid_fine.cell_counts);
		for(int j=0;j<number_of_components;j++){int r=i*number_of_components+j;ptr[r+1]=n;}}
	int nnz=0;for(int i=0;i<=R.rows();i++){nnz+=ptr[i];ptr[i]=nnz;}
	R.resizeNonZeros(nnz);
	real* val=R.valuePtr();
	int* col=R.innerIndexPtr();

	////Fill R
	Grid<d> nb_grid(VectorDi::Ones()*4);
	#pragma omp parallel for
	for(int R_i=0;R_i<m;R_i++){
		VectorDi cell=grid_coarse.Cell_Coord(R_i);int c=0;
		int nb_num=Pow(4,d);
		VectorDi cell_fine=coarsen_factor_vec.cwiseProduct(cell);
		iterate_cell(iter,nb_grid){VectorDi offset=iter.Coord();	
			VectorDi nb_fine=(cell_fine+offset-VectorDi::Ones());
			if(use_psi_P){nb_fine=grid_fine.Periodic_Cell(nb_fine);}
			else{if(AuxFunc::Outside<d>(grid_fine.cell_counts,nb_fine))continue;}
				
			real coef=Interpolation_Coefficient_Four_Points(offset);
			int R_j=grid_fine.Cell_Index(nb_fine);
			for(int nc=0;nc<number_of_components;nc++){
				int ii=number_of_components*R_i+nc;
				int jj=number_of_components*R_j+nc;
				int index=ptr[ii]+c;val[index]=coef;col[index]=jj;}c++;}}
}

template<int d,class T_SPARSE>
void Build_Restriction_Matrix_Three_Points(
	Vector<int,d>& cell_counts_fine,const int number_of_components,
	const int coarsen_factor,const bool use_psi_P,/*rst*/T_SPARSE& R)
{	Typedef_VectorDii(d);
	Grid<d> grid_fine(cell_counts_fine,(real)1);
	Grid<d> grid_coarse=grid_fine.Coarsened(coarsen_factor);
	VectorDi zero=VectorDi::Zero();
	VectorDi coarsen_factor_vec=grid_fine.cell_counts.cwiseQuotient(grid_coarse.cell_counts);

	int m=grid_coarse.Number_Of_Nodes();
	int n=grid_fine.Number_Of_Nodes();
	int rows=(m*number_of_components);
	int cols=(n*number_of_components);
	R.resize(rows,cols);
	int* ptr=R.outerIndexPtr();
	ptr[0]=0;

	#pragma omp parallel for
	for(int i=0;i<m;i++){
		VectorDi node=grid_coarse.Node_Coord(i);
		VectorDi node_fine=coarsen_factor_vec.cwiseProduct(node);
		int n=use_psi_P?Grid<d>::Number_Of_Nb_R():Grid<d>::Number_Of_Nb_RV(node_fine,grid_fine.node_counts);
		for(int j=0;j<number_of_components;j++){int r=i*number_of_components+j;ptr[r+1]=n;}}

	int nnz=0;for(int i=0;i<=R.rows();i++){nnz+=ptr[i];ptr[i]=nnz;}
	R.resizeNonZeros(nnz);
	real* val=R.valuePtr();
	int* col=R.innerIndexPtr();
	
	#pragma omp parallel for
	for(int R_i=0;R_i<m;R_i++){
		VectorDi node=grid_coarse.Node_Coord(R_i);int c=0;
		int nb_num=Grid<d>::Number_Of_Nb_R();
		for(int i=0;i<nb_num;i++){
			VectorDi node_fine=coarsen_factor_vec.cwiseProduct(node);
			VectorDi offset=Grid<d>::Nb_R(zero,i);
			VectorDi nb_fine=(node_fine+offset);
			if(use_psi_P){nb_fine=grid_fine.Periodic_Node(nb_fine);}
			else{if(AuxFunc::Outside<d>(grid_fine.node_counts,nb_fine))continue;}

			real coef=Interpolation_Coefficient_Three_Points(offset);
			int R_j=grid_fine.Node_Index(nb_fine);
			for(int nc=0;nc<number_of_components;nc++){
				int ii=number_of_components*R_i+nc;
				int jj=number_of_components*R_j+nc;
				int index=ptr[ii]+c;val[index]=coef;col[index]=jj;}c++;}}
}

#define Inst_Helper(d,T_SPARSE) \
template void Build_Restriction_Matrix_Four_Points<d,T_SPARSE>(Vector<int,d>&,const int,const int,const bool,T_SPARSE&); \
template void Build_Restriction_Matrix_Three_Points<d,T_SPARSE>(Vector<int,d>&,const int,const int,const bool,T_SPARSE&);
Inst_Helper(2,SparseMatrixT);
Inst_Helper(3,SparseMatrixT);
#undef Inst_Helper

//////////////////////////////////////////////////////////////////////////
////Build R for irregular domains
//////////////////////////////////////////////////////////////////////////
template<int d> void Coarsen_Cell_Material(const Grid<d>& grid_fine,const Grid<d>& grid_coarse,
	const int number_of_components,const Field<short,d>& material_id_fine,/*result*/Field<short,d>& material_id_coarse)
{	Typedef_VectorDii(d);using namespace AuxFunc;
	int m=grid_coarse.Number_Of_Cells();
	material_id_coarse.Resize(grid_coarse.cell_counts);
	VectorDi coarsen_factor_vec=Cwise_Dvid(grid_fine.cell_counts,grid_coarse.cell_counts);

	#pragma omp parallel for
	for(int i=0;i<m;i++){
		VectorDi cell=grid_coarse.Cell_Coord(i);int cell_mat_id_coarse=-1;
		VectorDi fine_cell_start=Cwise_Prod(cell,coarsen_factor_vec);
		VectorDi fine_cell_end=fine_cell_start+coarsen_factor_vec-VectorDi::Ones();
		iterate_range(iter,fine_cell_start,fine_cell_end){const VectorDi& fine_cell=iter.Coord();
			if(material_id_fine(fine_cell)!=-1){cell_mat_id_coarse=material_id_fine(fine_cell);break;}}
		material_id_coarse(cell)=cell_mat_id_coarse;}	
}

template<int d> int Number_Of_Valid_Node_Nbs(const Vector<int,d>& coord,const Field<int,d>& grid_node_to_matrix)
{
	int n=0;for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++){
		const Vector<int,d>& nb=Grid<d>::Nb_R(coord,i);
		if(Grid<d>::Valid(nb,grid_node_to_matrix.counts)&&grid_node_to_matrix(nb)!=-1)n++;}return n;
}

template<int d> int Number_Of_Valid_Cell_Nbs(const Vector<int,d>& coord,const Field<int,d>& grid_cell_to_matrix)
{
	Vector<int,d> start=Vector<int,d>::Ones()*(-1);Vector<int,d> end=Vector<int,d>::Ones()*2;
	int n=0;iterate_range(iter,start,end){const Vector<int,d>& nb=iter.Coord()+coord;
		if(Grid<d>::Valid(nb,grid_cell_to_matrix.counts)&&grid_cell_to_matrix(nb)!=-1)n++;}return n;
}

template<int d,class T_SPARSE> 
void Build_Restriction_Matrix_Three_Points_Irregular_Domain(
	Vector<int,d>& cell_counts_fine,const Array<short>& material_id_array_fine,
	const int number_of_components,const int coarsen_factor,
	/*result*/Array<short>& material_id_array_coarse,/*rst*/T_SPARSE& R)
{	Typedef_VectorDii(d);using namespace AuxFunc;
	Grid<d> grid_fine(cell_counts_fine,(real)1);
	Grid<d> grid_coarse=grid_fine.Coarsened(coarsen_factor);
	VectorDi coarsen_factor_vec=grid_fine.cell_counts.cwiseQuotient(grid_coarse.cell_counts);

	Field<short,d> material_id_fine(grid_fine.cell_counts,material_id_array_fine);
	Field<short,d> material_id_coarse(grid_coarse.cell_counts);

	Field<int,d> grid_to_matrix_fine;Array<int> matrix_to_grid_fine;
	std::function<bool(const int)> valid_fine_cell=[&material_id_fine](const int idx)->bool{return material_id_fine.array[idx]!=-1;};
	Build_Grid_Node_With_Incident_Cell_Matrix_Bijective_Mapping(grid_fine,valid_fine_cell,grid_to_matrix_fine,matrix_to_grid_fine);
	Coarsen_Cell_Material(grid_fine,grid_coarse,number_of_components,material_id_fine,material_id_coarse);

	Field<int,d> grid_to_matrix_coarse;Array<int> matrix_to_grid_coarse;
	std::function<bool(const int)> valid_coarse_cell=[&material_id_coarse](const int idx)->bool{return material_id_coarse.array[idx]!=-1;};
	Build_Grid_Node_With_Incident_Cell_Matrix_Bijective_Mapping(grid_coarse,valid_coarse_cell,grid_to_matrix_coarse,matrix_to_grid_coarse);
	material_id_array_coarse=material_id_coarse.array;

	int m=(int)matrix_to_grid_coarse.size();
	int n=(int)matrix_to_grid_fine.size();
	int rows=(m*number_of_components);
	int cols=(n*number_of_components);
	R.resize(rows,cols);
	int* ptr=R.outerIndexPtr();
	ptr[0]=0;
	VectorDi zero=VectorDi::Zero();

	#pragma omp parallel for
	for(int i=0;i<m;i++){int node_coarse_idx=matrix_to_grid_coarse[i];
		VectorDi node_coarse=grid_coarse.Node_Coord(node_coarse_idx);
		VectorDi node_fine=coarsen_factor_vec.cwiseProduct(node_coarse);
		int valid_nb_n=Number_Of_Valid_Node_Nbs<d>(node_fine,grid_to_matrix_fine);
		assert(valid_nb_n>0);
		for(int j=0;j<number_of_components;j++){
			int r=i*number_of_components+j;ptr[r+1]=valid_nb_n;}}

	int nnz=0;for(int i=0;i<=m;i++){nnz+=ptr[i];ptr[i]=nnz;}

	R.resizeNonZeros(nnz);
	real* val=R.valuePtr();
	int* col=R.innerIndexPtr();

	#pragma omp parallel for
	for(int i=0;i<m;i++){int node_coarse_idx=matrix_to_grid_coarse[i];
		VectorDi node_coarse=grid_coarse.Node_Coord(node_coarse_idx);int c=0;
		int nb_num=Grid<d>::Number_Of_Nb_R();
		for(int j=0;j<nb_num;j++){
			VectorDi node_fine=coarsen_factor_vec.cwiseProduct(node_coarse);
			VectorDi offset=Grid<d>::Nb_R(zero,j);
			VectorDi nb_fine=node_fine+offset;
			if(Outside<d>(grid_fine.node_counts,nb_fine)||grid_to_matrix_fine(nb_fine)==-1)continue;
			real coef=Interpolation_Coefficient_Three_Points(offset);
			int node_fine_mtx_idx=grid_to_matrix_fine(nb_fine);
			for(int nc=0;nc<number_of_components;nc++){
				int ii=number_of_components*i+nc;
				int jj=number_of_components*node_fine_mtx_idx+nc;
				int index=ptr[ii]+c;val[index]=coef;col[index]=jj;}
			c++;}}	
}

template<int d,class T_SPARSE> 
void Build_Restriction_Matrix_Four_Points_Irregular_Domain(
	Vector<int,d>& cell_counts_fine,const Array<short>& material_id_array_fine,
	const int number_of_components,const int coarsen_factor,
	/*result*/Array<short>& material_id_array_coarse,/*rst*/T_SPARSE& R)
{	Typedef_VectorDii(d);using namespace AuxFunc;
	Grid<d> grid_fine(cell_counts_fine,(real)1);
	Grid<d> grid_coarse=grid_fine.Coarsened(coarsen_factor);
	VectorDi coarsen_factor_vec=grid_fine.cell_counts.cwiseQuotient(grid_coarse.cell_counts);

	Field<short,d> material_id_fine(grid_fine.cell_counts,material_id_array_fine);
	Field<short,d> material_id_coarse(grid_coarse.cell_counts);

	Field<int,d> grid_to_matrix_fine;Array<int> matrix_to_grid_fine;
	std::function<bool(const int)> valid_fine_cell=[&material_id_fine](const int idx)->bool{return material_id_fine.array[idx]!=-1;};
	Build_Grid_Cell_Matrix_Bijective_Mapping(grid_fine,valid_fine_cell,grid_to_matrix_fine,matrix_to_grid_fine);
	Coarsen_Cell_Material(grid_fine,grid_coarse,number_of_components,material_id_fine,material_id_coarse);

	Field<int,d> grid_to_matrix_coarse;Array<int> matrix_to_grid_coarse;
	std::function<bool(const int)> valid_coarse_cell=[&material_id_coarse](const int idx)->bool{return material_id_coarse.array[idx]!=-1;};
	Build_Grid_Cell_Matrix_Bijective_Mapping(grid_coarse,valid_coarse_cell,grid_to_matrix_coarse,matrix_to_grid_coarse);
	material_id_array_coarse=material_id_coarse.array;

	int m=(int)matrix_to_grid_coarse.size();
	int n=(int)matrix_to_grid_fine.size();
	int rows=(m*number_of_components);
	int cols=(n*number_of_components);
	R.resize(rows,cols);
	int* ptr=R.outerIndexPtr();
	ptr[0]=0;
	VectorDi zero=VectorDi::Zero();

	////Allocate R
	#pragma omp parallel for
	for(int i=0;i<m;i++){int cell_coarse_idx=matrix_to_grid_coarse[i];
		VectorDi cell_coarse=grid_coarse.Cell_Coord(cell_coarse_idx);
		VectorDi cell_fine=coarsen_factor_vec.cwiseProduct(cell_coarse);
		int valid_nb_n=Number_Of_Valid_Cell_Nbs<d>(cell_fine,grid_to_matrix_fine);
		assert(valid_nb_n>0);
		for(int j=0;j<number_of_components;j++){
			int r=i*number_of_components+j;ptr[r+1]=valid_nb_n;}}
	int nnz=0;for(int i=0;i<=m;i++){nnz+=ptr[i];ptr[i]=nnz;}

	R.resizeNonZeros(nnz);
	real* val=R.valuePtr();
	int* col=R.innerIndexPtr();

	////Fill R
	Grid<d> nb_grid(VectorDi::Ones()*4);
	#pragma omp parallel for
	for(int R_i=0;R_i<m;R_i++){int cell_coarse_idx=matrix_to_grid_coarse[R_i];
		VectorDi cell_coarse=grid_coarse.Cell_Coord(cell_coarse_idx);int c=0;
		int nb_num=Pow(4,d);
		VectorDi cell_fine=coarsen_factor_vec.cwiseProduct(cell_coarse);
		iterate_cell(iter,nb_grid){VectorDi offset=iter.Coord();
			VectorDi nb_fine=(cell_fine+offset-VectorDi::Ones());
			if(Outside<d>(grid_fine.cell_counts,nb_fine)||grid_to_matrix_fine(nb_fine)==-1)continue;
			real coef=Interpolation_Coefficient_Four_Points(offset);
			int R_j=grid_to_matrix_fine(nb_fine);
			for(int nc=0;nc<number_of_components;nc++){
				int ii=number_of_components*R_i+nc;
				int jj=number_of_components*R_j+nc;
				int index=ptr[ii]+c;val[index]=coef;col[index]=jj;}c++;}}
}

#define Inst_Helper(d,T_SPARSE) \
template void Coarsen_Cell_Material<d>(const Grid<d>&,const Grid<d>&,const int,const Field<short,d>&,Field<short,d>&); \
template int Number_Of_Valid_Node_Nbs<d>(const Vector<int,d>&,const Field<int,d>&); \
template void Build_Restriction_Matrix_Three_Points_Irregular_Domain<d,T_SPARSE>(Vector<int,d>&,const Array<short>&,const int,const int,Array<short>&,T_SPARSE&); \
template void Build_Restriction_Matrix_Four_Points_Irregular_Domain<d,T_SPARSE>(Vector<int,d>&,const Array<short>&,const int,const int,Array<short>&,T_SPARSE&);
Inst_Helper(2,SparseMatrixT);
Inst_Helper(3,SparseMatrixT);
#undef Inst_Helper

//////////////////////////////////////////////////////////////////////////
////Member functions of GeometricMultiGridSystem
//////////////////////////////////////////////////////////////////////////

template<int d,class T,class T_SPARSE,class T_VECTOR>
void GeometricMultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Initialize(const T_SPARSE& _A,const Vector<int,d>& counts,
	const Params _params/*=Params()*/,const Field<short,d>* material_id/*=nullptr*/)
{
	params=_params;

	////initialize mg_system parameters that need to be recorded from the passed-in params
	one_over_intp=(T)1/params.relax_coef;
	dof_on_cell=params.dof_on_cell;
	block_size=params.block_size;

	////initialize levels
	if(params.use_auto_calculated_levels){
		int lv_coarsest_system=(int)std::ceil(std::log2(_A.nonZeros())-params.coarsest_size_log2);
		int lv_coarsest_grid=(int)std::ceil(std::log2(AuxFunc::Max(counts)))-1;
		levels=std::min(lv_coarsest_system,lv_coarsest_grid);
		levels=std::max(1,levels);
		std::cout<<"Auto mg level with #nonzero: "<<_A.nonZeros()<<", #levels: "<<levels<<" ["<<lv_coarsest_system<<", "<<lv_coarsest_grid<<"]"<<std::endl;}
	else levels=params.levels;
	
	////re-check use_irregular_domain to make sure it is necessary
	use_irregular_domain=(params.use_irregular_domain&&material_id!=nullptr&&material_id->Has(-1));

	A=new T_SPARSE*[levels];
	R=new T_SPARSE*[levels-1];
	P=new T_SPARSE*[levels-1];
	b=new T_VECTOR*[levels];for(int i=0;i<levels;i++)b[i]=nullptr;
	x=new T_VECTOR*[levels];for(int i=0;i<levels;i++)x[i]=nullptr;
	r=new T_VECTOR*[levels];for(int i=0;i<levels;i++)r[i]=nullptr;
	////b, x, r are allocated but not initialized
	own_data_b.resize(levels);AuxFunc::Fill(own_data_b,true);
	own_data_x.resize(levels);AuxFunc::Fill(own_data_x,true);

	own_data_A.resize(levels);
	if(init_A_levels)AuxFunc::Fill(own_data_A,true);
	else AuxFunc::Fill(own_data_A,false);

	////cell_counts for level-0, if dof are on nodes, counts is node_counts
	Vector<int,d> cell_counts_on_level=dof_on_cell?counts:counts-Vector<int,d>::Ones();	
	cell_counts.push_back(cell_counts_on_level);

	////A for level-0, set A[0] as the input pointer (no ownership)
	A[0]=const_cast<T_SPARSE*>(&_A);own_data_A[0]=false;
	int n0=(int)A[0]->rows();
	
	////resize x,b,r for level-0
	Resize(x[0],n0);
	Resize(b[0],n0);
	Resize(r[0],n0);
	
	////initialize irregular domain for levels-0
	if(use_irregular_domain){
		mat_id=new Field<short,d>*[levels];
		own_mat_id.resize(levels);AuxFunc::Fill(own_mat_id,true);
		mat_id[0]=const_cast<Field<short,d>*>(material_id);own_mat_id[0]=false;}

	for(int i=0;i<levels-1;i++){
		R[i]=new T_SPARSE();
		////build R between level i and i+1
		if(!use_irregular_domain){
			if(dof_on_cell) Build_Restriction_Matrix_Four_Points<d,T_SPARSE>(cell_counts_on_level,params.block_size,params.coarsen_factor,params.use_psi_P,*R[i]);
			else Build_Restriction_Matrix_Three_Points<d,T_SPARSE>(cell_counts_on_level,params.block_size,params.coarsen_factor,params.use_psi_P,*R[i]);}
		else{
			mat_id[i+1]=new Field<short,d>();
			mat_id[i+1]->counts=Coarsened<d>(cell_counts_on_level,params.coarsen_factor);
			if(dof_on_cell) Build_Restriction_Matrix_Four_Points_Irregular_Domain<d,T_SPARSE>
				(cell_counts_on_level,mat_id[i]->array,params.block_size,params.coarsen_factor,mat_id[i+1]->array,*R[i]);
			else Build_Restriction_Matrix_Three_Points_Irregular_Domain<d,T_SPARSE>
				(cell_counts_on_level,mat_id[i]->array,params.block_size,params.coarsen_factor,mat_id[i+1]->array,*R[i]);}
		////update P as R transpose
		P[i]=new T_SPARSE();
		*P[i]=R[i]->transpose();
		////update A_i+1=R*A_i*P
		if(init_A_levels){
			A[i+1]=new T_SPARSE();
			*A[i+1]=(*R[i])*(*A[i])*(*P[i])*one_over_intp;}
		////coarsening for next level
		cell_counts_on_level=Coarsened<d>(cell_counts_on_level,params.coarsen_factor);
		cell_counts.push_back(cell_counts_on_level);
		////resize for next level
		int n=A_size(i+1);
		Resize(x[i+1],n);
		Resize(b[i+1],n);
		Resize(r[i+1],n);}

	if(params.use_color){Initialize_Color();}
}

template<int d,class T,class T_SPARSE,class T_VECTOR>
void GeometricMultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Initialize_Color()
{
	color_n.resize(levels);
	color_ptr.resize(levels);
	color.resize(levels);
	if(use_irregular_domain){
		for(int i=0;i<levels;i++){
			Vector<int,d> counts=dof_on_cell?(cell_counts[i]):(cell_counts[i]+Vector<int,d>::Ones());
			color_n[i]=ColorGrid<d>::Number_Of_Colors();
			ColorGrid<d>::Color(counts,mat_id[i]->array,color_ptr[i],color[i]);}}
	else{
		for(int i=0;i<levels;i++){
			Vector<int,d> counts=dof_on_cell?(cell_counts[i]):(cell_counts[i]+Vector<int,d>::Ones());
			color_n[i]=ColorGrid<d>::Number_Of_Colors();
			ColorGrid<d>::Color(counts,color_ptr[i],color[i]);}}
}

template<int d,class T,class T_SPARSE,class T_VECTOR>
void GeometricMultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Update_A_Levels()
{
	////Assuming A[0] is updated
	T one_over_intp=(T)1/params.relax_coef;
	for(int i=0;i<levels-1;i++){
		*A[i+1]=(*R[i])*(*A[i])*(*P[i])*one_over_intp;}
}

template<int d,class T,class T_SPARSE,class T_VECTOR>
void GeometricMultiGrid::GeometricMultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Clear()
{
	MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Clear();
	if(mat_id){for(int i=0;i<levels;i++)if(own_mat_id[i]&&mat_id[i]!=nullptr)delete mat_id[i];delete[] mat_id;}
}

#define Inst_Helper(d,real) \
template class GeometricMultiGridSystem<d,real,SparseMatrix<real>,VectorN<real> >;
Inst_Helper(2,real);
Inst_Helper(3,real);
#undef Inst_Helper
};

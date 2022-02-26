//////////////////////////////////////////////////////////////////////////
// Multigrid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "MultiGrid.h"
#include <iostream>

namespace MultiGrid
{
////Helper functions
template<class T> void Resize(T*& val,int n){if(val!=nullptr)delete val;val=new T[(size_type)n];}
template<class T> void Resize(Array<T>*& val,int n){if(val==nullptr)val=new Array<T>();val->resize((size_type)n);}
template<class T> void Resize(VectorN<T>*& val,int n){if(val==nullptr)val=new VectorN<T>();val->resize((size_type)n);}
template<class T> void Delete(T*& val){delete [] val;}
template<class T> void Delete(Array<T>*& val){delete val;}
template<class T> void Delete(VectorN<T>*& val){delete val;}
template<class T> void Set_Zero(T*& val,int n){for(int i=0;i<n;i++)val[i]=(T)0;}
template<class T> void Set_Zero(Array<T>*& val,int n){for(int i=0;i<n;i++)(*val)[i]=(T)0;}
template<class T> void Set_Zero(VectorN<T>*& val,int n){for(int i=0;i<n;i++)(*val)[i]=(T)0;}

template<class T> void Gauss_Seidel_Smoothing(const SparseMatrix<T>& A,VectorN<T>& x,const VectorN<T>& b,const int iter_num=1,const bool reversed_order=false)
{
	for(int iter=0;iter<iter_num;iter++)for(int ii=0;ii<A.rows();ii++){
		int i=(reversed_order?(int)A.rows()-1-ii:ii);T delta=b(i);
		for(InnerIterator<T> iter(A,i);iter;++iter){
			int j=(int)iter.col();T A_ij=iter.value();delta-=A_ij*x(j);}
		delta/=A.coeff(i,i);x(i)+=delta;}
}

#define Inst_Helper(T) \
template void Resize(T*&,int); \
template void Resize(Array<T>*&,int); \
template void Resize(VectorN<T>*&,int); \
template void Delete(T*&); \
template void Delete(Array<T>*&); \
template void Delete(VectorN<T>*&); \
template void Set_Zero(T*&,int n); \
template void Set_Zero(Array<T>*&,int); \
template void Set_Zero(VectorN<T>*&,int); \
template void Gauss_Seidel_Smoothing<T>(const SparseMatrix<T>&,VectorN<T>&,const VectorN<T>&,const int,const bool);
Inst_Helper(float);
Inst_Helper(double);
#undef Inst_Helper

//////////////////////////////////////////////////////////////////////////
////Multigrid data structure
////initialize with a uniform grid
template<int d,class T,class T_SPARSE,class T_VECTOR>
void MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Set_X_And_B(T_VECTOR& _x,const T_VECTOR& _b)
{
	if(x[0]!=&_x&&own_data_x[0])Delete(x[0]);x[0]=&_x;own_data_x[0]=false;
	if(b[0]!=&_b&&own_data_b[0])Delete(b[0]);b[0]=const_cast<T_VECTOR*>(&_b);own_data_b[0]=false;
}

template<int d,class T,class T_SPARSE,class T_VECTOR>
void MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::Clear()
{
	if(A){for(int i=0;i<levels;i++)if(own_data_A[i]&&A[i]!=nullptr)delete A[i];delete[] A;}
	if(P){for(int i=0;i<levels-1;i++)delete P[i];delete[] P;}
	if(R){for(int i=0;i<levels-1;i++)delete R[i];delete[] R;}
	if(AP){for(int i=0;i<levels-1;i++)delete AP[i];delete[] AP;}

	if(!on_device){
		if(b){for(int i=0;i<levels;i++)if(own_data_b[i]&&b[i]!=nullptr)Delete(b[i]);delete[] b;}
		if(x){for(int i=0;i<levels;i++)if(own_data_x[i]&&x[i]!=nullptr)Delete(x[i]);delete[] x;}
		if(r){for(int i=0;i<levels;i++)if(r[i]!=nullptr)Delete(r[i]);delete[] r;}}
	//else{
		//if(b){for(int i=0;i<levels;i++)Free_On_Device(b[i]);delete[] b;}
		//if(x){for(int i=0;i<levels;i++)Free_On_Device(x[i]);delete[] x;}
		//if(r){for(int i=0;i<levels;i++)Free_On_Device(r[i]);delete[] r;}}	
}

//////////////////////////////////////////////////////////////////////////
////This function is used for the current CPU multigrid solver (GMGPCG_Solver_CPU)
////TODO: use multicolor GS to accelerate
template<int d,class T,class T_SPARSE,class T_VECTOR>
void MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::V_Cycle()
{
	for(int i=0;i<=levels-2;i++){
		////smooth A_i*x_i=b_i
		if(i>0)Set_Zero(x[i],(int)A[i]->rows());
		Gauss_Seidel_Smoothing(*A[i],*x[i],*b[i],params.smooth_num,false);
		////r_i=b_i-A_i*x_i
		(*r[i])=(*b[i])-(*A[i])*(*x[i]);
		////b_{i+1}=R_i*r_i
		(*b[i+1])=(*R[i])*(*r[i]);}

	////solve A_{n-1}x_{n-1}=b_{n-1} with a direct solver
	SparseSolver::LU(*A[levels-1],*x[levels-1],*b[levels-1]);
			
	for(int i=levels-2;i>=0;i--){
		////x_i=x_i+P_i*x_{i+1}
		(*x[i])=(*x[i])+(*P[i])*(*x[i+1]);
		////smooth A_i*x_i=b_i
		Gauss_Seidel_Smoothing(*A[i],*x[i],*b[i],params.smooth_num,true);}
}

#define Inst_Helper(d,T) \
template class MultiGridSystem<d,T,SparseMatrix<T>,VectorN<T> >;
Inst_Helper(2,real);
Inst_Helper(3,real);
#undef Inst_Helper
};
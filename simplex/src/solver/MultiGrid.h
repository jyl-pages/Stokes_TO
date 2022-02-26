//////////////////////////////////////////////////////////////////////////
// Multigrid for sparse matrix
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __MultiGrid_h__
#define __MultiGrid_h__
#include "Common.h"
#include "SparseFunc.h"

namespace MultiGrid{

//////////////////////////////////////////////////////////////////////////
////Multigrid parameters
class Params
{public:
	int levels=3;
	int v_cycle_num=1;
	int des_iter_num=1;
	int asc_iter_num=1;
	int max_cg_iter_num=1000;
	int smooth_num=1; 
	int block_size=1;			////number of DoFs for each element
	int coarsen_factor=2;
	real tolerance=(real)1e-5;
	real jacobian_weight=(real)2/(real)3;
	real relax_coef=(real)1;
	int coarsest_size_log2=12;	////4096 for coarsest A by default, effective only when use_auto_calculated_levels is true

	////customizable parameters
	bool dof_on_cell=true;		////whether DoF is on cell or on node
	bool use_auto_calculated_levels=true;
	bool use_color=false;
	bool use_psi_P=false;
	bool use_irregular_domain=false;
	bool use_gpu=false;
	bool init_hier_on_gpu=true;	////this flag is only meaningful for cpu
};

//////////////////////////////////////////////////////////////////////////
////CPU multigrid data structure base class
template<int d,class T,class T_SPARSE,class T_VECTOR> class MultiGridSystem
{public:
	Params params;
	////matrix hierarchies
	int levels=1;
	T_SPARSE** A=nullptr;
	T_SPARSE** P=nullptr;
	T_SPARSE** R=nullptr;
	T_SPARSE** AP=nullptr;
	T_VECTOR** b=nullptr;
	T_VECTOR** x=nullptr;
	T_VECTOR** r=nullptr;
	bool on_device=false;
	Array<bool> own_data_A;
	Array<bool> own_data_x;
	Array<bool> own_data_b;

	////blocks and colors
	int block_size=1;
	Array<int> color_n;				////number of colors
	Array<Array<int> > color_ptr;	////pointing to the first index of each color, size = number of colors; the outer array is for all levels
	Array<Array<int> > color;		////all indices grouped for each color, starting from color 0; the outer array is for all levels

	MultiGridSystem(){}
	~MultiGridSystem(){Clear();}

	virtual void Set_X_And_B(T_VECTOR& _x,const T_VECTOR& _b);
	virtual void Clear();
	virtual void V_Cycle();
};

//////////////////////////////////////////////////////////////////////////
////Helper functions
template<class T> void Resize(T*& val,int n);
template<class T> void Resize(Array<T>*& val,int n);
template<class T> void Resize(VectorN<T>*& val,int n);
template<class T> void Delete(T*& val);
template<class T> void Delete(Array<T>*& val);
template<class T> void Delete(VectorN<T>*& val);
template<class T> void Set_Zero(T*& val,int n);
template<class T> void Set_Zero(Array<T>*& val,int n);
template<class T> void Set_Zero(VectorN<T>*& val,int n);
};
#endif

//////////////////////////////////////////////////////////////////////////
// Geometric Multigrid on CPU
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __GeometricMultiGrid_h__
#define __GeometricMultiGrid_h__
#include "MultiGrid.h"
#include "KrylovSolver.h"
#include "Field.h"
#include "Timer.h"

namespace GeometricMultiGrid{
using namespace MultiGrid;

//////////////////////////////////////////////////////////////////////////
////CPU geometric multigrid data structure
template<int d,class T,class T_SPARSE,class T_VECTOR> class GeometricMultiGridSystem : public MultiGridSystem<d,T,T_SPARSE,T_VECTOR>
{using Base=MultiGridSystem<d,T,T_SPARSE,T_VECTOR>;
public:
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::params;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::levels;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::A;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::P;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::R;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::AP;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::b;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::x;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::r;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::own_data_A;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::own_data_x;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::own_data_b;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::block_size;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::color_n;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::color_ptr;
	using MultiGridSystem<d,T,T_SPARSE,T_VECTOR>::color;
	
	////parameters
	Array<Vector<int,d> > cell_counts;
	T one_over_intp=(T)1;
	Field<short,d>** mat_id=nullptr;
	Array<bool> own_mat_id;

	////flags
	bool dof_on_cell=false;				////If not, dofs are on grid nodes.
	bool init_A_levels=true;			////This flag specifies whether to calculate the A hierarchy on CPU or not. It calculates A on CPU if true. If not, the A levels will be directly calculated on GPU.
	bool use_irregular_domain=false;

	virtual void Clear();
	virtual void V_Cycle(){Base::V_Cycle();}	////TODO: parallelize G-S using multiple colors

	////Initialize the MG hierarchy on CPU
	void Initialize(const T_SPARSE& _A,const Vector<int,d>& counts,const Params _params=Params(),const Field<short,d>* material_id=nullptr);
	////Initialize color for each level on CPU
	void Initialize_Color();
	////update levels of A with constant R and P on CPU: A_i+1=R_i*A_i*P_i
	void Update_A_Levels();		
	////helper function for system without a host storage for A
	int A_size(int level) const {if(level==0)return (int)A[0]->rows();else return (int)R[level-1]->rows();}
};

//////////////////////////////////////////////////////////////////////////
////Multigrid preconditioned Krylov solver
template<int d,class T,class KA,class KX> class KrylovPreGMG
{public:
	const KA* ka=nullptr;
	Vector<int,d> counts;
	GeometricMultiGridSystem<d,T,KA,KX> multigrid;
	Params params;

	KrylovPreGMG(const KA* _ka,const Vector<int,d>& _counts,const Params& _params,const Field<short,d>* mat_id=nullptr)
		:ka(_ka),counts(_counts),params(_params){multigrid.Initialize(*ka,counts,params,mat_id);}

	void Initialize() 
	{
		multigrid.Initialize(*ka,counts,params);
	}

	template<class KXW> void Precond(const KXW& r,KXW& rst)
	{
		multigrid.Set_X_And_B(*rst.kx,*r.kx);
		multigrid.V_Cycle();
	}
};
};
#endif

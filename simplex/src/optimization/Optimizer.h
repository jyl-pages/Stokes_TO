//#####################################################################
// Optimizer interface
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __Optimizer_h__
#define __Optimizer_h__

#include <float.h>          // FLT_MAX
class Optimizer
{
public:
	int n_var;	////var num
	int n_cons;	////constraint num
	bool use_var_L;
	bool use_var_U;
	bool write_opt_substep;

	////bound values
	real var_lb;		////lower bound for var
	real var_ub;		////upper bound for var
	real cons_lb;		////lower bound for constraints
	real cons_ub;		////upper bound for constraints

	bool write_intmed;
	real intmed_obj;
	Array<real> intmed_var;
	Array<real> intmed_grad;
	Array<real> intmed_cons;

	Optimizer():n_var(0),n_cons(0),use_var_L(false),use_var_U(false),write_opt_substep(true),
		var_lb((real)-FLT_MAX),var_ub((real)FLT_MAX),cons_lb((real)-FLT_MAX),cons_ub((real)FLT_MAX),write_intmed(true),intmed_obj((real)0){}
	
	virtual void Initialize_Optimizer(){}
	virtual void Optimize(){}
	virtual void Write_Substep(const int frame){}

	virtual void Set_Var_Lower_Bound(const real val){var_lb=val;use_var_L=true;}
	virtual void Set_Var_Upper_Bound(const real val){var_ub=val;use_var_U=true;}
	virtual void Set_Cons_Lower_Bound(const real val){cons_lb=val;}
	virtual void Set_Cons_Upper_Bound(const real val){cons_ub=val;}
};
#endif
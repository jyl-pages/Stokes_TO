//#####################################################################
// Optimizer IpOpt interface
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifdef USE_IPOPT
#ifndef __OptimizerIpOpt_h__
#define __OptimizerIpOpt_h__
#include "Common.h"
#include "Optimizer.h"

struct IpoptProblemInfo;

class IntermediateDataIpOpt
{public:
	int alg_mod;int iter_count;real obj_value;real inf_pr;real inf_du;
	real mu;real d_norm;real regularization_size;real alpha_du;real alpha_pr;int ls_trials;
	IntermediateDataIpOpt(){}
	IntermediateDataIpOpt(int _alg_mod,int _iter_count,real _obj_value,real _inf_pr,real _inf_du,real _mu,real _d_norm,real _regularization_size,real _alpha_du,real _alpha_pr,int _ls_trials);
};

class FinalDataIpOpt
{public:
	FinalDataIpOpt();
	int n_var;
	int n_cons;
	real obj;
	Array<real> var;
	Array<real> lmda_cons;
	Array<real> lmda_var_L;
	Array<real> lmda_var_U;

	void Copy_From(const int _n_var,const int _n_cons,const real _obj,const real* _var,const real* _lmda_cons,const real* _lmda_var_L,const real* _lmda_var_U);
	bool Copy_To(const int _n_var,const int _n_cons,real* _var,real* _lmda_cons,real* _lmda_var_L,real* _lmda_var_U);
};

class OptimizerIpOpt : public Optimizer
{typedef Optimizer Base;
public:
	real* var;
	real* var_L;
	real* var_U;
	real* cons_L;
	real* cons_U;
	real* lmda_cons;
	real* lmda_var_L;
	real* lmda_var_U;


	bool use_elementwise_bound = false;
	real* var_LB = nullptr;		////preset lower bound
	real* var_UB = nullptr;		////preset upper bound

	int nnz_cons_jac;
	int n_hess;
	real obj;

	IntermediateDataIpOpt intmed_data;
	FinalDataIpOpt final_data;
	bool use_acceptable_tol;
	bool write_final;

	////ipopt params
	real tol;
	int max_iter_num;
	int acceptable_iter_num;
	real acceptable_tol_change;
	bool accept_every_trial_step;
	enum SolveStatus{Uninited=-1,Succeed=0,Acceptable,MaxIter,Failed} solve_status;
	bool use_warm_start;

	OptimizerIpOpt();
	~OptimizerIpOpt();
	virtual void Optimize();
	virtual real Compute_Objective(const real* var){return (real)0;}
	virtual void Compute_Gradient(const real* var,real* grad){}
	virtual void Compute_Constraint(const real* var,real* constraint){if(n_cons==0)return;}
	virtual void Resize_Constraint_Jacobian(int* row_num,int* col_num){if(n_cons==0)return;}
	virtual void Compute_Constraint_Jacobian(const real* var,real* constraint_jacobian){if(n_cons==0)return;}
	virtual void Allocate_Hessian(int* iRow,int* jCol){}
	virtual void Compute_Hessian(const real* var,real* hess_values){}
protected:
	template<class T_VAL> T_VAL* New(const int n){if(n==0)return nullptr;else return (T_VAL*)malloc(sizeof(T_VAL)*n);}
	template<class T_VAL> void Delete(T_VAL* & ptr){if(ptr!=nullptr){free(ptr);ptr=nullptr;}}
	virtual void Allocate_Data();
	virtual void Set_Boundary();
	virtual void Delete_Data();
	virtual void Initialize_IpOpt_Parameters(IpoptProblemInfo* opt);
	void Optimize_IpOpt();
public:

	////ipopt interfaces
	static int Ipopt_Objective(int n_var,real* var,int new_var,real* obj_value,void* user_data);
	static int Ipopt_Gradient(int n_var,real* var,int new_var,real* grad,void* user_data);
	static int Ipopt_Constraint(int n_var,real* var,int new_var,int n_cons,real* cons,void* user_data);
	static int Ipopt_Constraint_Jacobian(int n_var,real* var,int new_var,int n_cons,int nele_jac,int* iRow,int* jCol,real* values,void* user_data);
	static int Ipopt_Hessian(int n_var,real* var,int new_var,real obj_factor,int n_cons,real* lambda,int new_lambda,int n_hess,
		int* iRow,int* jCol,real* values,void* user_data);
	static int Ipopt_Intermediate(int alg_mod,int iter_count,real obj_value,
		real inf_pr,real inf_du,real mu,real d_norm,real regularization_size,real alpha_du,real alpha_pr,int ls_trials,void* user_data);
};
#endif
#endif
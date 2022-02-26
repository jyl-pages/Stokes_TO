//#####################################################################
// Optimizer NlOpt
// Bo Zhu, MIT, 10/16
//#####################################################################
#ifndef __OptimizerNlOpt_h__
#define __OptimizerNlOpt_h__
#ifdef USE_NLOPT
#include "Common.h"
#include "Optimizer.h"
#include "nlopt.h"

class OptimizerNlOpt : public Optimizer
{typedef Optimizer Base;
public:
	real* var;
	real* var_L;
	real* var_U;
	real* cons_tol;
	real tol;
	real c_tol;
	real obj;
	int iter_count;

	bool write_intmed;
	real intmed_obj;
	Array<real> intmed_var;
	Array<real> intmed_grad;
	Array<real> intmed_cons;

	OptimizerNlOpt():Base(),var(nullptr),var_L(nullptr),var_U(nullptr),cons_tol(nullptr),tol((real)1e-3),c_tol((real)1e-3),
		obj((real)FLT_MAX),iter_count(0),write_intmed(true),intmed_obj((real)0){}
	~OptimizerNlOpt(){Delete_Data();}

	virtual void Optimize(){Optimize_NlOpt();}

	virtual real Compute_Objective(const real* var){return (real)0;}
	virtual void Compute_Gradient(const real* var,real* grad){}
	virtual void Compute_Constraint(const real* var,real* constraint){if(n_cons==0)return;}
	virtual void Resize_Constraint_Jacobian(int* row_num,int* col_num){if(n_cons==0)return;}
	virtual void Compute_Constraint_Jacobian(const real* var,real* constraint_jacobian){if(n_cons==0)return;}
	virtual void Compute_Hessian(){}
	virtual void Write_Substep(const int frame){}

protected:
	template<class T_VAL> T_VAL* New(const int n){if(n==0)return nullptr;else return (T_VAL*)malloc(sizeof(T_VAL)*n);}
	template<class T_VAL> void Delete(T_VAL* & ptr){if(ptr!=nullptr){free(ptr);ptr=nullptr;}}

	virtual void Allocate_Data()
	{
		////assuming n_cons, n_var already initialized
		Delete(var);var=New<real>(n_var);
		Delete(var_L);if(use_var_L)var_L=New<real>(n_var);
		Delete(var_U);if(use_var_U)var_U=New<real>(n_var);
		Delete(cons_tol);cons_tol=New<real>(n_cons);						
		if(n_cons>0){for(int i=0;i<n_cons;i++)cons_tol[i]=c_tol;}

		if(write_intmed){
			intmed_var.resize(n_var);
			intmed_grad.resize(n_var);
			intmed_cons.resize(n_var);}
	}

	virtual void Delete_Data()
	{
		Delete(var);Delete(var_L);Delete(var_U);Delete(cons_tol);
	}

	virtual void Initialize_NlOpt_Parameters(nlopt_opt& opt){}

	void Optimize_NlOpt()
	{
		nlopt_opt opt=nlopt_create(NLOPT_LD_MMA,n_var);
		if(var_L!=nullptr)nlopt_set_lower_bounds(opt,var_L);
		if(var_U!=nullptr)nlopt_set_upper_bounds(opt,var_U);
		if(n_cons>0)nlopt_add_inequality_mconstraint(opt,n_cons,Nlopt_Constraint,this,cons_tol);
		nlopt_set_min_objective(opt,Nlopt_Objective,this);
		nlopt_set_ftol_rel(opt,tol);
		
		if(nlopt_optimize(opt,var,&obj)<0){std::cout<<"nlopt failed"<<std::endl;}
		else{std::cout<<"nlopt succeed"<<std::endl;}
		if(write_intmed){std::memcpy(intmed_var.data(),var,n_var*sizeof(real));}
		nlopt_destroy(opt);
		Delete(var_L);Delete(var_U);
	}

public:
	static real Nlopt_Objective(unsigned n_var,const real* var,real* grad,void *user_data)
	{
		OptimizerNlOpt* opt=(OptimizerNlOpt*)(user_data);
		if(opt==0){std::cout<<"invalid instance for OptimizerNlOpt"<<std::endl;return (real)0;}
		real obj_value=opt->Compute_Objective(var);
		if(opt->write_intmed){std::memcpy(opt->intmed_var.data(),var,n_var*sizeof(real));opt->intmed_obj=obj_value;}
		if(grad!=0){
			opt->Compute_Gradient(var,grad);
			if(opt->write_intmed)std::memcpy(opt->intmed_grad.data(),grad,n_var*sizeof(real));}
		opt->Write_Substep(opt->iter_count++);
		return obj_value;
	}

	static void Nlopt_Constraint(unsigned n_cons,real* cons,unsigned n_var,const real* var,real* cons_jac,void* user_data)
	{
		OptimizerNlOpt* opt=(OptimizerNlOpt*)(user_data);
		if(opt==0){std::cout<<"invalid instance for OptimizerNlOpt"<<std::endl;return;}	
		if(opt->n_cons>0){
			opt->Compute_Constraint(var,cons);
			opt->Compute_Constraint_Jacobian(var,cons_jac);
		}
	}
};
#endif
#endif
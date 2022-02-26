//#####################################################################
// Optimizer IpOpt interface
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifdef USE_IPOPT
#include <iostream>
#include "Common.h"
#include "AuxFunc.h"
#include "OptimizerIpOpt.h"
#include "IpStdCInterface.h"

IntermediateDataIpOpt::IntermediateDataIpOpt(int _alg_mod,int _iter_count,real _obj_value,real _inf_pr,real _inf_du,real _mu,real _d_norm,real _regularization_size,real _alpha_du,real _alpha_pr,int _ls_trials)
:alg_mod(_alg_mod),iter_count(_iter_count),obj_value(_obj_value),inf_pr(_inf_pr),inf_du(_inf_du),
mu(_mu),d_norm(_d_norm),regularization_size(_regularization_size),alpha_du(_alpha_du),alpha_pr(_alpha_pr),ls_trials(_ls_trials){}

FinalDataIpOpt::FinalDataIpOpt():n_var(0),n_cons(0),obj((real)0){}

void FinalDataIpOpt::Copy_From(const int _n_var,const int _n_cons,const real _obj,const real* _var,const real* _lmda_cons,const real* _lmda_var_L,const real* _lmda_var_U)
{n_var=_n_var;n_cons=_n_cons;obj=_obj;var.resize(n_var);lmda_var_L.resize(n_var);lmda_var_U.resize(n_var);lmda_cons.resize(n_cons);
for(int i=0;i<n_var;i++){var[i]=_var[i];lmda_var_L[i]=_lmda_var_L[i];lmda_var_U[i]=_lmda_var_U[i];}for(int i=0;i<n_cons;i++)lmda_cons[i]=_lmda_cons[i];}

bool FinalDataIpOpt::Copy_To(const int _n_var,const int _n_cons,real* _var,real* _lmda_cons,real* _lmda_var_L,real* _lmda_var_U)
{if(n_var!=_n_var||n_cons!=_n_cons)return false;
for(int i=0;i<n_var;i++){_var[i]=var[i];_lmda_var_L[i]=lmda_var_L[i];_lmda_var_U[i]=lmda_var_U[i];}return true;}

OptimizerIpOpt::OptimizerIpOpt()
:Base(),var(nullptr),var_L(nullptr),var_U(nullptr),cons_L(nullptr),cons_U(nullptr),
lmda_cons(nullptr),lmda_var_L(nullptr),lmda_var_U(nullptr),nnz_cons_jac(0),n_hess(0),obj((real)FLT_MAX),
use_acceptable_tol(false),write_final(true),tol((real)1e-5),max_iter_num(1000),acceptable_iter_num(10),
acceptable_tol_change((real)1e-4), accept_every_trial_step(false), solve_status(Uninited), use_warm_start(false) {}

OptimizerIpOpt::~OptimizerIpOpt()
{Delete_Data();}

void OptimizerIpOpt::Optimize()
{Optimize_IpOpt();}

void OptimizerIpOpt::Allocate_Data()
{
	////assuming n_cons, n_var already initialized
	Delete(var);var=New<real>(n_var);
	
	if (use_elementwise_bound) {
		Delete(var_LB); var_LB = New<real>(n_var);
		Delete(var_UB); var_UB = New<real>(n_var);
	}
	Delete(var_L); var_L = New<real>(n_var);
	Delete(var_U); var_U = New<real>(n_var);
	
	Delete(cons_L);cons_L=New<real>(n_cons);						
	Delete(cons_U);cons_U=New<real>(n_cons);						
	for(int i=0;i<n_cons;i++){cons_L[i]=cons_lb;cons_U[i]=cons_ub;}

	Delete(lmda_cons);lmda_cons=New<real>(n_cons);
	Delete(lmda_var_L);lmda_var_L=New<real>(n_var);
	Delete(lmda_var_U);lmda_var_U=New<real>(n_var);

	if(write_intmed){
		intmed_var.resize(n_var);AuxFunc::Fill(intmed_var,(real)0);
		intmed_grad.resize(n_var);AuxFunc::Fill(intmed_grad,(real)0);
		intmed_cons.resize(n_cons);AuxFunc::Fill(intmed_cons,(real)0);}
}

void OptimizerIpOpt::Set_Boundary() {
	if (use_elementwise_bound) {
		for (int i = 0; i < n_var; i++) { var_L[i] = var_LB[i]; var_U[i] = var_UB[i]; }
	}
	else {
		for (int i = 0; i < n_var; i++) { var_L[i] = var_lb; var_U[i] = var_ub; }
	}
}

void OptimizerIpOpt::Delete_Data()
{
	Delete(var);Delete(var_L);Delete(var_U);Delete(cons_L);Delete(cons_U);
	if (use_elementwise_bound) {
		Delete(var_LB);
		Delete(var_UB);
	}
	Delete(lmda_cons);Delete(lmda_var_L);Delete(lmda_var_U);
}

void OptimizerIpOpt::Initialize_IpOpt_Parameters(IpoptProblemInfo* opt)
{
	std::cout<<"use acc tol: "<<use_acceptable_tol<<", acc iter num: "<<acceptable_iter_num<<
		", acc tol: "<<acceptable_tol_change<<", acc trial: "<<accept_every_trial_step<<", use warm: "<<use_warm_start<<std::endl;

	AddIpoptNumOption(opt,"tol",tol);
	AddIpoptIntOption(opt,"max_iter",max_iter_num);

	if(use_acceptable_tol){
		AddIpoptNumOption(opt,"acceptable_tol",1e10);		////1e-6
		AddIpoptIntOption(opt,"acceptable_iter",acceptable_iter_num);		////15
		AddIpoptNumOption(opt,"acceptable_constr_viol_tol",1e10);	////.01
		AddIpoptNumOption(opt,"acceptable_dual_inf_tol",1e10);	////1e10
		AddIpoptNumOption(opt,"acceptable_compl_inf_tol",1e10);	////.01
		AddIpoptNumOption(opt,"acceptable_obj_change_tol",acceptable_tol_change);}
		
	if(accept_every_trial_step)AddIpoptStrOption(opt,"accept_every_trial_step","yes");
	AddIpoptIntOption(opt,"max_soc",0);
	char* exact_s = "exact";
	char* limited_s = "limited-memory";
	char* ans;
	if(n_hess>0)
		ans = exact_s;
	else
		ans = limited_s; 
	AddIpoptStrOption(opt,"hessian_approximation",ans);
	if(use_warm_start)AddIpoptStrOption(opt,"warm_start_init_point","yes");

	AddIpoptIntOption(opt,"print_level",0);
	AddIpoptStrOption(opt,"print_user_options","yes");
	//AddIpoptStrOption(opt,"print_timing_statistics","yes");	
}

void OptimizerIpOpt::Optimize_IpOpt()
{
	std::cout<<"ipopt, n_var: "<<n_var<<", n_cons: "<<n_cons<<", n_hess: "<<n_hess<<std::endl;
	solve_status=Uninited;

	IpoptProblem opt=CreateIpoptProblem(n_var,var_L,var_U,n_cons,cons_L,cons_U,nnz_cons_jac,n_hess,0,
		&Ipopt_Objective,&Ipopt_Constraint,&Ipopt_Gradient,&Ipopt_Constraint_Jacobian,&Ipopt_Hessian);
	//std::cout<<"var_L: ";for(int i=0;i<n_var;i++)std::cout<<var_L[i]<<", ";std::cout<<std::endl;
	//std::cout<<"var_U: ";for(int i=0;i<n_var;i++)std::cout<<var_U[i]<<", ";std::cout<<std::endl;
	//std::cout<<"cons_L:";for(int i=0;i<n_cons;i++)std::cout<<cons_L[i]<<", ";std::cout<<std::endl;
	//std::cout<<"cons_U:";for(int i=0;i<n_cons;i++)std::cout<<cons_U[i]<<", ";std::cout<<std::endl;
	//std::cout<<"nnz_cons_jac: "<<nnz_cons_jac<<std::endl;
	Delete(var_L);Delete(var_U);Delete(cons_L);Delete(cons_U);
	Initialize_IpOpt_Parameters(opt);
	SetIntermediateCallback(opt,Ipopt_Intermediate); 
		
	//std::cout<<"var: ";for(int i=0;i<n_var;i++)std::cout<<var[i]<<", ";std::cout<<std::endl;
	if(use_warm_start&&final_data.n_var!=0){final_data.Copy_To(n_var,n_cons,var,lmda_cons,lmda_var_L,lmda_var_U);
		for(int i=0;i<n_cons;i++)std::cout<<lmda_cons[i]<<", ";std::cout<<std::endl;}
	enum ApplicationReturnStatus status=IpoptSolve(opt,var,nullptr,&obj,lmda_cons,lmda_var_L,lmda_var_U,this);
	if(use_warm_start||write_final)final_data.Copy_From(n_var,n_cons,obj,var,lmda_cons,lmda_var_L,lmda_var_U);

	if(status==Solve_Succeeded){solve_status=Succeed;std::cout<<"ipopt converged in "<<intmed_data.iter_count<<" iters"<<std::endl;}
	else if(status==Solved_To_Acceptable_Level){solve_status=Acceptable;std::cout<<"ipopt solved to acceptable level in "<<intmed_data.iter_count<<" iters"<<std::endl;}
	else if(status==Maximum_Iterations_Exceeded){solve_status=MaxIter;std::cout<<"ipopt max iterations exceeds in "<<intmed_data.iter_count<< " iters"<<std::endl;}
	else{solve_status=Failed;std::cout<<"ipopt solve failed"<<std::endl;}

	if(write_intmed){std::memcpy(intmed_var.data(),var,n_var*sizeof(real));/*Write_Substep(intmed_data.iter_count+1);*/}
	FreeIpoptProblem(opt);
	Delete(lmda_cons);Delete(lmda_var_L);Delete(lmda_var_U);
}

int OptimizerIpOpt::Ipopt_Objective(int n_var,real* var,int new_var,real* obj_value,void* user_data)
{
	OptimizerIpOpt* opt=(OptimizerIpOpt*)(user_data);
	if(opt==0){std::cout<<"invalid instance for OptimizerIpOpt"<<std::endl;return 0;}
	obj_value[0]=opt->Compute_Objective(var);
	if(opt->write_intmed){
		std::memcpy(opt->intmed_var.data(),var,n_var*sizeof(real));
		opt->intmed_obj=obj_value[0];}
	return 1;	
}

int OptimizerIpOpt::Ipopt_Gradient(int n_var,real* var,int new_var,real* grad,void* user_data)
{
	OptimizerIpOpt* opt=(OptimizerIpOpt*)(user_data);
	if(opt==0){std::cout<<"invalid instance for OptimizerIpOpt"<<std::endl;return 0;}
	opt->Compute_Gradient(var,grad);
	if(opt->write_intmed)std::memcpy(opt->intmed_grad.data(),grad,n_var*sizeof(real));
	return 1;	
}

int OptimizerIpOpt::Ipopt_Constraint(int n_var,real* var,int new_var,int n_cons,real* cons,void* user_data)
{
	if(n_cons==0)return 0;
	OptimizerIpOpt* opt=(OptimizerIpOpt*)(user_data);
	if(opt==0){std::cout<<"invalid instance for OptimizerIpOpt"<<std::endl;return 0;}
	opt->Compute_Constraint(var,cons);
	if(opt->write_intmed)std::memcpy(opt->intmed_cons.data(),cons,n_cons*sizeof(real));
	return 1;
}

int OptimizerIpOpt::Ipopt_Constraint_Jacobian(int n_var,real* var,int new_var,int n_cons,int nele_jac,int* iRow,int* jCol,real* values,void* user_data)
{
	if(n_cons==0)return 0;
	OptimizerIpOpt* opt=(OptimizerIpOpt*)(user_data);
	if(opt==0){std::cout<<"invalid instance for OptimizerIpOpt"<<std::endl;return 0;}
	if (values==nullptr){opt->Resize_Constraint_Jacobian(iRow,jCol);}
	else{
		opt->Compute_Constraint_Jacobian(var,values);
		//std::memcpy(defo->jg.data(),values,n_var*sizeof(real));
	}
	return 1;
}

int OptimizerIpOpt::Ipopt_Hessian(int n_var,real* var,int new_var,real obj_factor,int n_cons,real* lambda,int new_lambda,int n_hess,int* iRow,int* jCol,real* values,void* user_data)
{
	if(n_hess==0)return 0;
	OptimizerIpOpt* opt=(OptimizerIpOpt*)(user_data);
	if(opt==0){std::cout<<"invalid instance for OptimizerIpOpt"<<std::endl;return 0;}

	if(values==0)opt->Allocate_Hessian(iRow,jCol);
	else opt->Compute_Hessian(var,values);
	return 1;
}

int OptimizerIpOpt::Ipopt_Intermediate(int alg_mod,int iter_count,real obj_value,real inf_pr,real inf_du,real mu,real d_norm,real regularization_size,real alpha_du,real alpha_pr,int ls_trials,void* user_data)
{
	//Seperation();std::cout<<"iteration "<<iter_count;Seperation();std::cout<<std::endl;
	OptimizerIpOpt* opt=(OptimizerIpOpt*)(user_data);
	if(opt==0){std::cout<<"invalid instance for OptimizerIpOpt"<<std::endl;return 0;}
	if(opt->write_intmed)opt->intmed_data=
		IntermediateDataIpOpt(alg_mod,iter_count,obj_value,inf_pr,inf_du,mu,d_norm,regularization_size,alpha_du,alpha_pr,ls_trials);
	if(opt->write_opt_substep)opt->Write_Substep(iter_count);
	return 1;	
}
#endif
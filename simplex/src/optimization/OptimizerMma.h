//#####################################################################
// Optimizer MMA
//#####################################################################
#ifndef __OptimizerMma_h__
#define __OptimizerMma_h__
#include "Common.h"
#include "AuxFunc.h"
#include "Optimizer.h"
#include "MMASolver.h"
#include "Timer.h"
#include <cstdlib>
#include <iostream>

class OptimizerMMA : public Optimizer
{typedef Optimizer Base;
public:
	real* var=nullptr;			////variable
	real* var_L=nullptr;		////current lower bound
	real* var_U=nullptr;		////current upper bound

	bool use_elementwise_bound=false;
	real* var_LB=nullptr;		////preset lower bound
	real* var_UB=nullptr;		////preset upper bound

	real* grad=nullptr;
	real* cons=nullptr;
	real* cons_grad=nullptr;

	real tol=(real)1e-3;
	real obj=(real)FLT_MAX;
	int iter_count=0;
	int max_iter_num=1000;
    real movlim=(real).2;

	bool write_intmed=true;
	real intmed_obj=(real)0;

	Array<real> intmed_var;
	Array<real> intmed_grad;
	Array<real> intmed_cons;
	Array<real> intmed_cons_grad;

	bool verbose=false;

	OptimizerMMA():Base(){}
	~OptimizerMMA(){Delete_Data();}

	virtual void Optimize(){Optimize_MMA();}
	virtual real Compute_Objective(const real* var){return (real)0;}
	virtual void Compute_Gradient(const real* var,real* grad){}
	virtual void Compute_Constraint(const real* var,real* constraint){if(n_cons==0)return;}
	virtual void Compute_Constraint_Grad(const real* var,real* constraint_grad){}
	virtual void Write_Substep(const int frame){}

protected:
	template<class T_VAL> T_VAL* New(const int n){if(n==0)return nullptr;else return (T_VAL*)malloc(sizeof(T_VAL)*n);}
	template<class T_VAL> void Delete(T_VAL* & ptr){if(ptr!=nullptr){free(ptr);ptr=nullptr;}}

	virtual void Allocate_Data()
	{
		////assuming n_cons,n_var already initialized
		Delete(var);var=New<real>(n_var);
		Delete(var_L);var_L=New<real>(n_var);
		Delete(var_U);var_U=New<real>(n_var);

		if(use_elementwise_bound){
			Delete(var_LB);var_LB=New<real>(n_var);
			Delete(var_UB);var_UB=New<real>(n_var);}

		Delete(grad);grad=New<real>(n_var);
		Delete(cons);cons=New<real>(n_cons);
		Delete(cons_grad);cons_grad=New<real>(n_var*n_cons);

		intmed_var.resize(n_var);
		if(write_intmed){
			intmed_grad.resize(n_var);
			intmed_cons.resize(n_cons);
			intmed_cons_grad.resize(n_var*n_cons);}

		iter_count=0;
	}

	virtual void Delete_Data()
	{
		Delete(var);
		Delete(var_L);
		Delete(var_U);
		if(use_elementwise_bound){
			Delete(var_LB);
			Delete(var_UB);}
		Delete(grad);
		Delete(cons);
		Delete(cons_grad);
	}

	void Optimize_MMA()
	{
    	MMASolver *mma=new MMASolver(n_var,n_cons);
		mma->ConstraintModification(true);

		////assuming var_LB and var_UB were set outside the class
		real ch=1.0;
        while(ch>tol && iter_count<max_iter_num){
			if(verbose)std::cout<<"============================== MMA loop "<<iter_count<<" start =============================="<<std::endl;
			obj=MMA_Objective(n_var,var,grad);
			MMA_Constraint(n_cons,cons,cons_grad,n_var,var);
			Write_Substep(iter_count);
            //// Set outer move limits
			if(use_elementwise_bound){
				for(int i=0;i<n_var;i++){
					var_U[i]=AuxFunc::Min(var_UB[i],var[i]+movlim);
					var_L[i]=AuxFunc::Max(var_LB[i],var[i]-movlim);}}
			else{
				for(int i=0;i<n_var;i++){
					var_U[i]=AuxFunc::Min(var_ub,var[i]+movlim);
					var_L[i]=AuxFunc::Max(var_lb,var[i]-movlim);}}
			//// Update MMA next step
			MMA_Intermediate(mma);
            //// Compute inf norm on design change
            ch=0.0;for(int i=0;i<n_var;i++){ch=std::max(ch,std::abs(var[i]-intmed_var[i]));}

			if(verbose)std::cout<<"============================== MMA loop "<<iter_count<<" end =============================="<<std::endl;
			iter_count++;
		}

        // Deallocate
		Delete(mma);
	}

	void Set_Var_LB(const real* _var_lb)
	{
		for(int i=0;i<n_var;i++){
			var_LB[i]=_var_lb[i];}
		use_var_L=true;
	}

	void Set_Var_UB(const real* _var_ub)
	{
		for(int i=0;i<n_var;i++){
			var_UB[i]=_var_ub[i];}
		use_var_U=true;
	}

	void Set_Var_LB_And_UB(const real* _var_lb,const real* _var_ub)
	{
		Set_Var_LB(_var_lb);
		Set_Var_UB(_var_ub);
	}

public:
	void MMA_Intermediate(MMASolver* mma)
	{
        // Call the update method
        mma->Update(var,grad,cons,cons_grad,var_L,var_U);
	}

	real MMA_Objective(unsigned n_var,const real* var,real* grad)
	{
		real obj_value=Compute_Objective(var);
		if(write_intmed){
			std::memcpy(&intmed_var[0],var,n_var*sizeof(real));
			intmed_obj=obj_value;}

		if(grad!=0){
			Compute_Gradient(var,grad);
			if(write_intmed)std::memcpy(&intmed_grad[0],grad,n_var*sizeof(real));}

		return obj_value;
	}

	void MMA_Constraint(unsigned n_cons,real* cons,real* cons_grad,unsigned n_var,const real* var)
	{
		if(n_cons>0){
			Compute_Constraint(var,cons);
			Compute_Constraint_Grad(var,cons_grad);
			if(write_intmed){
				std::memcpy(&intmed_cons[0],cons,n_cons*sizeof(real));
				std::memcpy(&intmed_cons_grad[0],cons_grad,(n_cons*n_var)*sizeof(real));}}
	}
};
#endif

//#####################################################################
// Optimizer SnOpt
// Bo Zhu, MIT, 04/17
//#####################################################################
#ifndef __OptimizerSnOpt_h__
#define __OptimizerSnOpt_h__
#ifdef USE_SNOPT
#include "Common.h"
#include "LinearAlgebra.h"
#include "Optimizer.h"
#include "snopt.hh"
#include "snoptProblem.hh"
class SnoptSolver : public snoptProblem
{typedef snoptProblem Base;typedef long long int TLI;
public:
	SnoptSolver():Base(){}

	char*& Cw(){return cw;}
	const char* Cw() const {return cw;}
	TLI LenCw() const {return lencw;}
};

class OptimizerSnOpt : public Optimizer
{typedef Optimizer Base;typedef long long int TLI;
public:
	real* var;
	real* var_L;
	real* var_U;
	real* var_mul;
	TLI* var_state;
	TLI n_F;
	real* F;
	real* F_L;
	real* F_U;
	real* F_mul;
	TLI* F_state;

	int nnz_cons_jac;
	int n_G;
	TLI* G_i;
	TLI* G_j;

	char* var_names;
	char* F_names;

	real tol;
	real c_tol;
	real obj;
	TLI iter_count;

	SnoptSolver snopt;

	OptimizerSnOpt():Base(),var(nullptr),var_L(nullptr),var_U(nullptr),var_mul(nullptr),var_state(nullptr),
		n_F(1),F(nullptr),F_L(nullptr),F_U(nullptr),F_mul(nullptr),F_state(nullptr),nnz_cons_jac(0),n_G(0),G_i(nullptr),G_j(nullptr),
		var_names(nullptr),F_names(nullptr),tol((real)1e-3),c_tol((real)1e-3),obj((real)FLT_MAX),iter_count(0){}
	~OptimizerSnOpt(){Delete_Data();}

	virtual void Optimize(){Optimize_SnOpt();}

	virtual real Compute_Objective(const real* var){return (real)0;}
	virtual void Compute_Gradient(const real* var,real* grad){}
	virtual void Compute_Constraint(const real* var,real* constraint){if(n_cons==0)return;}
	virtual void Resize_Constraint_Jacobian(TLI* row_num,TLI* col_num){if(n_cons==0)return;}	////fill matrix index of cons_jac
	virtual void Resize_Constraint_Jacobian(){if(nnz_cons_jac>0)Resize_Constraint_Jacobian(&G_i[n_var],&G_j[n_var]);}
	virtual void Compute_Constraint_Jacobian(const real* var,real* constraint_jacobian){if(n_cons==0)return;}
	virtual void Compute_Hessian(){}
	virtual void Write_Substep(const int frame){}

protected:
	template<class T_VAL> T_VAL* New(const int n){if(n==0)return nullptr;else return new T_VAL[n];}
	template<class T_VAL> T_VAL* New(const TLI n){if(n==0)return nullptr;else return new T_VAL[n];}
	template<class T_VAL> void Delete(T_VAL* & ptr){if(ptr!=nullptr){delete ptr;ptr=nullptr;}}

	virtual void Allocate_Data()
	{
		////assuming n_cons, n_var are initialized
		Delete(var);var=New<real>(n_var);
		Delete(var_L);var_L=New<real>(n_var);
		Delete(var_U);var_U=New<real>(n_var);
		Delete(var_mul);var_mul=New<real>(n_var);
		Delete(var_state);var_state=New<TLI>(n_var);
		
		n_F=(TLI)(1+n_cons);Delete(F);F=New<real>(n_F);
		Delete(F_L);F_L=New<real>(n_F);
		Delete(F_U);F_U=New<real>(n_F);
		
		Delete(F_mul);F_mul=New<real>(n_F);
		Delete(F_state);F_state=New<TLI>(n_F);

		////assuming nnz_cons_jac is initialized
		n_G=n_var+nnz_cons_jac;
		Delete(G_i);G_i=New<TLI>(n_G);
		Delete(G_j);G_j=New<TLI>(n_G);
		//std::cout<<"n_var: "<<n_var<<", n_cons: "<<n_cons<<", nnz_cons_jac: "<<nnz_cons_jac<<", n_F: "<<n_F<<", n_G: "<<n_G<<std::endl;

		Delete(var_names);var_names=New<char>(8);
		Delete(F_names);F_names=New<char>(8);

		if(write_intmed){
			intmed_var.resize(n_var);Fill(intmed_var,(real)0);
			//ignore intmed_grad and intmed_cons
		}
	}

	virtual void Delete_Data()
	{
		Delete(var);Delete(var_L);Delete(var_U);Delete(var_mul);Delete(var_state);
		Delete(F);Delete(F_L);Delete(F_U);Delete(F_mul);Delete(F_state);
		Delete(G_i);Delete(G_j);Delete(var_names);Delete(F_names);
	}

	virtual void Initialize_Data()
	{
		for(int i=0;i<n_var;i++){var_L[i]=var_lb;var_U[i]=var_ub;}

		if(n_G>=n_var)for(int i=0;i<n_var;i++){G_i[i]=0;G_j[i]=i;}
		for(int i=0;i<n_var;i++){var_L[i]=var_lb;var_U[i]=var_ub;var_state[i]=0;}
		F_L[0]=(real)-FLT_MAX;F_U[0]=(real)FLT_MAX;
		for(int i=0;i<n_cons;i++){F_L[i+1]=cons_lb;F_U[i+1]=cons_ub;}
		for(int i=0;i<n_F;i++){F_state[i]=0;}
	}

	void Optimize_SnOpt()
	{
		OptimizerSnOpt* ptr=this;TLI ptr_size=sizeof(ptr);
		//std::cout<<"ptr_size: "<<ptr_size<<", char: "<<sizeof(char)<<", sizeof(double): "<<sizeof(double)
		//	<<", long int: "<<sizeof(long int)<<", long long int: "<<sizeof(long long int)<<", std::uintptr_t: "<<sizeof(std::uintptr_t)<<std::endl;

		std::uintptr_t dptr=(std::uintptr_t)(void*)(ptr);
		std::memcpy(snopt.Cw(),&dptr,sizeof(dptr));

		//for(int i=0;i<n_G;i++)std::cout<<G_i[i]<<", "<<G_j[i]<<std::endl;

		snopt.setProbName("snopt");
		//snopt.setSpecsFile("snopt_spec");
		snopt.setProblemSize(n_var,n_F);
		snopt.setObjective(0,0);
		snopt.setA(0,nullptr,nullptr,nullptr);
		snopt.setG(n_G,G_i,G_j);
		snopt.setX(var,var_L,var_U,var_mul,var_state);
		snopt.setF(F,F_L,F_U,F_mul,F_state);
		snopt.setXNames(var_names,1);
		snopt.setFNames(F_names,1);
		snopt.setNeA(0);
		snopt.setNeG(n_G);
		snopt.setUserFun(Snopt_Callback);
		snopt.setIntParameter("Derivative option",1);
		//snopt.setIntParameter("Iteration limit",100);
		//snopt.setIntParameter("Minor iterations limit",10);
		//snopt.setIntParameter("Major iterations limit",10);
		//snopt.setRealParameter("Major feasibility tolerance",1);
		//snopt.setRealParameter("Major step limit",1e-2);
		//snopt.setRealParameter("Linesearch tolerance",1e-1);
		//snopt.setIntParameter("Hessian updates",1);
		snopt.setPrintFile("snopt_print");
		std::cout<<"finish setting snopt params"<<std::endl;	////

		//snopt.setIntParameter("Summary file",1);

		TLI cold=0,basis=1,warm=2;
		snopt.solve(cold);

		//std::cout<<"var: ";for(int i=0;i<n_var;i++)std::cout<<var[i]<<", ";std::cout<<std::endl;
	}

public:
	static int Snopt_Callback(TLI* status,TLI* n,real x[],TLI* needF,TLI* neF,real F[],TLI* needG,TLI* neG,real G[],
			 char* cu,TLI* lencu,TLI iu[],TLI* leniu,real ru[],TLI* lenru)
	{
		//static int iter=0;std::cout<<"frame "<<iter++<<", status: "<<*status<<std::endl;
		std::uintptr_t dptr;std::memcpy(&dptr,cu,sizeof(dptr));
		OptimizerSnOpt* ptr=reinterpret_cast<OptimizerSnOpt*>((void*)dptr);

		////compute objective and constraints
		if(*needF>0){
			std::cout<<"obj"<<std::endl;
			F[0]=ptr->Compute_Objective(x);
			if(ptr->n_cons>0)ptr->Compute_Constraint(x,&F[1]);}

		////compute gradient
		if(*needG>0){
			std::cout<<"grad"<<std::endl;
			ptr->Compute_Gradient(x,G);
			if(ptr->n_cons>0){ptr->Compute_Constraint_Jacobian(x,&G[ptr->n_var]);}}

		////write output
		if(*needF>0){
			if(ptr->write_intmed){
				std::memcpy(ptr->intmed_var.data(),x,ptr->n_var*sizeof(real));
				//for(int i=0;i<ptr->n_var;i++)std::cout<<x[i]<<", ";std::cout<<std::endl;
				ptr->intmed_obj=F[0];
				/*std::memcpy(ptr->intmed_grad.data(),G,ptr->n_var*sizeof(real));*/}
			if(ptr->write_opt_substep){static int frame=0;ptr->Write_Substep(frame++);}}

		return 0;
	}
};
#endif
#endif
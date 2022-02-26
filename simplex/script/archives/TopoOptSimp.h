//#####################################################################
// Topology optimization SIMP
// Copyright (c) (2018-), Bo Zhu, boolzhu@gmail.com
// This file is part of SLAX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __TopoOptSimp_h__
#define __TopoOptSimp_h__
#include "AuxFunc.h"
#include "Field.h"
#include "SoftBodyLinearFemGrid.h"
#include "Timer.h"

template<int d,bool use_material_id=false> class TopoOptSimp
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	SoftBodyLinearFemGrid<d>* soft_body=nullptr;
	Grid<d> grid;
	Field<real,d> x;
	real frac=(real).35;
	real p=(real)3;
	int max_iter_num=200;
	real rho_min=(real)1e-4;
	real total_volume=(real)0;
	real r=(real)1/(real)32;

	Field<real,d> dc;
	Field<real,d> fem_variable_coef;	////x^p, output to soft_body
	Field<real,d> cell_vol;				////cell volumes
	Field<real,d> e;
	std::function<void(const int)> Write_Output_Files_Callback;

	bool verbose=true;

	////check material macro
	#define Check_Mat_By_Cell(cell) \
	if constexpr(use_material_id) if(soft_body->material_id(cell)==-1)continue;
	#define Check_Mat_By_Value(mat_id) \
	if constexpr(use_material_id) if(mat_id==-1)continue;

	void Initialize(SoftBodyLinearFemGrid<d>* _soft_body)
	{
		soft_body=_soft_body;
		grid=soft_body->grid;
		x.Resize(grid.cell_counts,(real)0);
		dc.Resize(grid.cell_counts,(real)0);

		fem_variable_coef.Resize(grid.cell_counts);fem_variable_coef.Fill((real)0);
		soft_body->variable_coef=&fem_variable_coef;
		e.Resize(grid.cell_counts);e.Fill((real)0);

		real cell_volume=pow(grid.dx,d);
		total_volume=(real)0;cell_vol.Resize(grid.cell_counts);cell_vol.Fill((real)0);
		iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();Check_Mat_By_Cell(cell);
			cell_vol(cell)=cell_volume;total_volume+=cell_vol(cell);}	
	}

	virtual void Update_FEM()
	{
		static Timer<real> timer;
		////update soft_body variable coef
		int n_cell=grid.Number_Of_Cells();
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i=0;i<n_cell;i++){const VectorDi& cell=grid.Cell_Coord(i);Check_Mat_By_Cell(cell);
			real coef=pow(x(cell),p);
			coef=AuxFunc::Clamp(coef,rho_min,(real)1);	////clamp x^p
			fem_variable_coef(cell)=coef;}
		
		soft_body->u.setZero();
		soft_body->f.setZero();
		SparseFunc::Set_Value(soft_body->K,(real)0);
		if(verbose)timer.Reset();
		soft_body->Update_K_And_f();
		if(verbose){AuxFunc::Seperation_Line();timer.Elapse_And_Output_And_Reset("FEM update K and f");AuxFunc::Seperation_Line();}
		soft_body->Solve();
		if(verbose){AuxFunc::Seperation_Line();timer.Elapse_And_Output_And_Reset("FEM linear solve");AuxFunc::Seperation_Line();}
	}

	real Obj()
	{
		Update_FEM();
		soft_body->Compute_Elastic_Energy(e);
		real sum=(real)0;for(const auto& e0:e.array)sum+=e0;return sum;
	}

	virtual void Fill_X(Field<real,d>& _x,const real frac)
	{
		if constexpr(use_material_id) _x.Fill(frac);
		else{
			iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
				int mat_id=Material_Id(cell);
				if(mat_id==-1)_x(cell)=(real)0;else _x(cell)=frac;}}
	}

	real Total_Volume(const Field<real,d>& x) const
	{real v=(real)0;for(size_type i=0;i<x.array.size();i++)v+=x.array[i]*cell_vol.array[i];return v;}

	virtual void Update_DC()
	{
		//iterate_cell(iter,grid){
		//	VectorDi cell=iter.Coord();int mat_id=Material_Id(cell);Check_Mat_By_Value(mat_id);
		int n_cell=grid.Number_Of_Cells();
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i=0;i<n_cell;i++){const VectorDi& cell=grid.Cell_Coord(i);
			int mat_id=Material_Id(cell);Check_Mat_By_Value(mat_id);
			const MatrixX& K0=soft_body->Ke[mat_id];
			VectorX cell_u;soft_body->Compute_Cell_Displacement(soft_body->u,cell,cell_u);
			dc(cell)=-p*pow(x(cell),p-1)*cell_u.transpose()*K0*cell_u;}
	}

	void Set_Smooth_R(const real ratio=(real)1/(real)32)
	{
		r=std::max(grid.Length()[0]*ratio,grid.dx*(real)1.5);
	}

	virtual void Smooth_DC()
	{
		Field<real,d> dc_smoothed=dc;
		int n=(int)(std::floor(r/grid.dx));
		Grid<d> sub_grid(VectorDi::Ones()*(n*2+1),grid.dx);VectorDi sub_offset=VectorDi::Ones()*n;

		//iterate_cell(iter,grid) {const VectorDi& cell=iter.Coord();Check_Mat_By_Cell(cell);
		int n_cell=grid.Number_Of_Cells();
		#ifdef _OPENMP
		#pragma omp parallel for
		#endif
		for(int i=0;i<n_cell;i++){const VectorDi& cell=grid.Cell_Coord(i);Check_Mat_By_Cell(cell);
			real w=(real)0;real sum=(real)0;
			//real r=grid.dx*(real)1.5;
			//for(int i=0;i<Grid<d>::Number_Of_Nb_R();i++){
				//VectorDi nb_cell=grid.Nb_R(cell,i);
			////mesh independent filter
			iterate_cell(iter_sub,sub_grid){const VectorDi& nb_cell=iter_sub.Coord()+cell-sub_offset;
				if(grid.Valid_Cell(nb_cell)&&Is_Valid_Material(nb_cell)){
					real d0=(grid.Center(cell)-grid.Center(nb_cell)).norm();
					real w0=std::max(r-d0,(real)0);
					sum+=w0*x(nb_cell)*dc(nb_cell);
					w+=w0;}}
			w*=x(cell);if(w!=(real)0)dc_smoothed(cell)=sum/w;}
		dc=dc_smoothed;
	}

	virtual void Initialize_X()
	{
		Fill_X(x,frac);
	}

	virtual void Update_X()
	{
		real l_min=(real)0;real l_max=(real)1e10;real move=(real).1;real tol=(real)1e-4;
		Field<real,d> x_new(x.counts);
		real neg_dc_min=-AuxFunc::Max(dc.array);real neg_dc_max=-AuxFunc::Min(dc.array);
		while(l_max-l_min>tol){real l_mid=(real).5*(l_min+l_max);
			iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();Check_Mat_By_Cell(cell);
				real neg_dc=(real)0;
				if(neg_dc_min<(real)0&&neg_dc_max>0){
					neg_dc=(-dc(cell))/(-neg_dc_min)+(real)1;}
				else{neg_dc=(-dc(cell)-neg_dc_min)/(neg_dc_max-neg_dc_min)*(real)2;}
				real x2=x(cell);real b=sqrt(neg_dc/l_mid);
				x2=std::min(std::min(x(cell)+move,b*x(cell)),(real)1);	////upper bound
				x2=std::max(std::max(x(cell)-move,x2),(real).01);		////lower bound
				x_new(cell)=x2;}
			if(Total_Volume(x_new)-frac*total_volume>0)l_min=l_mid;else l_max=l_mid;}
		x=x_new;
	}

	virtual void Optimize()
	{
		Timer<real> timer;
		Initialize_X();
		real change=(real)1;real tolerance=(real)1e-3;
		int iter=0;
		while(change>tolerance&&iter<max_iter_num){
			if(verbose)timer.Reset();
			Field<real,d> x_old=x;
			Update_FEM();
			if(verbose)timer.Elapse_And_Output_And_Reset("FEM total");
			if(Write_Output_Files_Callback)Write_Output_Files_Callback(iter);
			if(verbose)timer.Elapse_And_Output_And_Reset("IO");
			Update_DC();
			Smooth_DC();
			Update_X();
			change=AuxFunc::Difference_Linf(x.array,x_old.array);
			if(verbose)timer.Elapse_And_Output_And_Reset("Others");
			iter++;}
	}

protected:
	int Material_Id(const VectorDi& cell) const
	{if constexpr(use_material_id) return soft_body->material_id(cell);else return 0;}

	int Material_Id(const int cell_idx) const
	{if constexpr(use_material_id) return soft_body->material_id(cell_idx);else return 0;}

	bool Is_Valid_Material(const VectorDi& cell) const
	{if constexpr(use_material_id) return soft_body->material_id(cell)!=-1;else return true;}

	////For debug purpose
	void Numerical_Derivative()
	{
		Field<real,d> g=x;
		Field<real,d> x0=x;real dx=(real).0001;Field<real,d> e;
		for(auto i=0;i<x.array.size();i++){
			real x1=x.array[i];x.array[i]=x1+dx;
			real e1=Obj();
			x.array[i]=x1-dx;
			real e2=Obj();
			g.array[i]=(e1-e2)/((real)2*dx);
			x.array[i]=x1;}
		Update_FEM();
		Update_DC();
		for(auto i=0;i<g.array.size();i++){std::cout<<"["<<g.array[i]<<", "<<dc.array[i]<<"] ";}std::cout<<std::endl;
	}
};
#endif

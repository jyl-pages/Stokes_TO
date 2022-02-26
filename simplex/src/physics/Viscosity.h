//////////////////////////////////////////////////////////////////////////
// Solve the viscosity equation
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

////TODO: Neumann boundary conditions

#ifndef __Viscosity_h__
#define __Viscosity_h__
#include "MacGrid.h"
#include "FaceField.h"
#include "BoundaryCondition.h"

template<int d> class Viscosity
{Typedef_VectorDii(d);
public:
	MacGrid<d>& mac_grid;
	FaceField<real,d>& velocity;
	FaceField<real,d>& alpha;
	Field<ushort,d>& type;
	BoundaryConditionMacGridViscosity<d>& bc;
	bool use_alpha_for_correction=false;	////TODO: handle alpha in viscosity

	Viscosity(MacGrid<d>& _mac_grid,FaceField<real,d>& _velocity,FaceField<real,d>& _alpha,Field<ushort,d>& _type,BoundaryConditionMacGridViscosity<d>& _bc):
		mac_grid(_mac_grid),velocity(_velocity),alpha(_alpha),type(_type),bc(_bc){}

	////implicit solve
	void Solve(const real dt,const real nu)
	{
		bc.Enforce_Boundary_Conditions(velocity);
		Solve_Velocity_On_Cell_Center(dt,nu);
		bc.Enforce_Boundary_Conditions(velocity);
	}

protected:////Helper functions
	void Solve_Velocity_On_Cell_Center(const real dt,const real nu)
	{
		SparseMatrixT A;VectorX u,b;	////Solve Au=b
		Field<int,d> grid_to_matrix;Array<int> matrix_to_grid;
		grid_to_matrix.Resize(mac_grid.grid.cell_counts,-1);
		
		////Find valid cells
		int n=0;iterate_cell(iter,mac_grid.grid){const VectorDi& cell=iter.Coord();
			if(Is_Fluid_Cell(cell)){grid_to_matrix(cell)=n++;
				matrix_to_grid.push_back(mac_grid.grid.Cell_Index(cell));}}
		
		////Setup A and b
		A.resize(n,n);u.resize(n);u.fill((real)0);b.resize(n);b.fill((real)0);

		real A_dia_coef_inv=pow(mac_grid.grid.dx,2)/(dt*nu);
		Array<TripletT> elements;
		for(auto r=0;r<matrix_to_grid.size();r++){
			const VectorDi& cell=mac_grid.grid.Cell_Coord(matrix_to_grid[r]);
			////off-diagonal elements
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
				VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(Is_Valid_Cell(nb_cell)&&Is_Fluid_Cell(nb_cell)){	
					int c=grid_to_matrix(nb_cell);
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(i,axis,side);
					VectorDi face=cell+VectorDi::Unit(axis)*side;
					real a=alpha(axis,face);
					elements.push_back(TripletT((int)r,(int)c,(real)-a));}}
			////diagonal elements
			real dia_coef=(real)0;
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
				real coef=Dia_Coef(cell,i);dia_coef+=coef;}
			dia_coef+=A_dia_coef_inv;
			elements.push_back(TripletT((int)r,(int)r,dia_coef));}

		A.setFromTriplets(elements.begin(),elements.end());
		A.makeCompressed();		

		Field<VectorD,d> cell_velocity(mac_grid.grid.cell_counts,VectorD::Zero());
		Interpolation<d> intp(mac_grid);
		////map from face to cell
		intp.Interpolate_Faces_To_Cells(velocity,cell_velocity);
		for(int axis=0;axis<d;axis++){
			////update u and rhs
			u.fill((real)0);
			////map from cell to vec
			for(int i=0;i<n;i++){const VectorD& v=cell_velocity.array[matrix_to_grid[i]];b[i]=v[axis]*A_dia_coef_inv;}
			////update bc
			for(auto& iter:bc.psi_D_values){int a=iter.first[0];
				if(a==axis){
					int face_idx=iter.first[1];VectorDi face=mac_grid.Face_Coord(axis,face_idx);real psi_D_value=iter.second;
					for(int i=0;i<2;i++){
						const VectorDi nb_cell=mac_grid.Face_Incident_Cell(axis,face,i);
						if(Is_Valid_Cell(nb_cell)){int i=grid_to_matrix(nb_cell);b[i]+=psi_D_value;}}}}
			
			////solve
			KrylovSolver::ICPCG(A,u,b);
			////map from vec to cell
			for(int i=0;i<n;i++){VectorD& v=cell_velocity.array[matrix_to_grid[i]];v[axis]=u[i];}}
		intp.Interpolate_Cells_To_Faces(cell_velocity,velocity);
	}

	real Dia_Coef(const VectorDi& cell,const int i)
	{
		int axis;VectorDi face;MacGrid<d>::Cell_Incident_Face(cell,i,axis,face);
		/*if(bc.Is_Psi_N(axis,face))return (real)0;else */return alpha(axis,face);	////TODO: handle psi_N
	}

	bool Is_Valid_Cell(const VectorDi& cell) const {return mac_grid.grid.Valid_Cell(cell);}
	bool Is_Fluid_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Fluid;}
};

#endif

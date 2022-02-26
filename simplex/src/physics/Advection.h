//////////////////////////////////////////////////////////////////////////
// Advection schemes (semi-Lagrangian)
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Advection_h__
#define __Advection_h__
#include <Field.h>
#include <FaceField.h>
#include <Interpolation.h>
#include <AuxFunc.h>

namespace Advection
{
////advect a velocity field
template<int d> void Semi_Lagrangian(const real dt,
	const MacGrid<d>& ghost_mac_grid,const FaceField<real,d>& ghost_velocity,
	const MacGrid<d>& mac_grid,FaceField<real,d>& velocity,
	std::function<bool(const int,const Vector<int,d>&)> valid=nullptr)
{	Typedef_VectorDii(d);
	Interpolation<d> ghost_intp(ghost_mac_grid);
	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
		#pragma omp parallel for
		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
			if(valid&&!valid(axis,face))continue;	////skip the invalid faces

			VectorD pos=mac_grid.face_grids[axis].Node(face);
			VectorD vel=ghost_intp.Interpolate_Face_Vectors(ghost_velocity,pos);
			VectorD mid_pos=pos-vel*(real).5*dt;
			vel=ghost_intp.Interpolate_Face_Vectors(ghost_velocity,mid_pos);
			VectorD backtraced_pos=pos-vel*dt;
			real advected_v=ghost_intp.Interpolate_Faces(ghost_velocity,backtraced_pos,axis);
			velocity(axis,face)=advected_v;}}
}

////advect a vector field with referred vector values
////particle position is calculated via ghost_velocity
////old vector value is calculated via ghost_vector
template<int d> void Semi_Lagrangian(const real dt,const FaceField<real,d>& velocity,	/*velocity is on ghost_mac_grid*/
	const MacGrid<d>& ghost_mac_grid,const FaceField<real,d>& ghost_vector,
	const MacGrid<d>& mac_grid,FaceField<real,d>& vector,
	bool use_zero_extrapolation=false,std::function<bool(const int,const Vector<int,d>&)> valid=nullptr)
{	Typedef_VectorDii(d);
	Interpolation<d> intp(ghost_mac_grid);
	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
		#pragma omp parallel for
		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
			if(valid&&!valid(axis,face))continue;	////skip the invalid faces

			VectorD pos=mac_grid.face_grids[axis].Node(face);
			VectorD vel=intp.Interpolate_Face_Vectors(velocity,pos);
			VectorD mid_pos=pos-vel*(real).5*dt;
			vel=intp.Interpolate_Face_Vectors(velocity,mid_pos);
			VectorD backtraced_pos=pos-vel*dt;
			real advected_v;
			if(use_zero_extrapolation&&(!ghost_mac_grid.grid.Inside(backtraced_pos)))advected_v=(real)0;
			else advected_v=intp.Interpolate_Faces(ghost_vector,backtraced_pos,axis);
			vector(axis,face)=advected_v;}}
}

////advect a cell-based field
template<class T,int d> void Semi_Lagrangian(const real dt,const FaceField<real,d>& velocity,
	const MacGrid<d>& ghost_mac_grid,const Field<T,d>& ghost_density,
	const MacGrid<d>& mac_grid,Field<T,d>& density,
	bool use_zero_extrapolation=false,std::function<bool(const Vector<int,d>&)> valid=nullptr)
{	Typedef_VectorDii(d);
	Interpolation<d> intp(ghost_mac_grid);
	int cell_num=mac_grid.grid.Number_Of_Cells();
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){VectorDi cell=mac_grid.grid.Cell_Coord(i);
		if(valid&&!valid(cell))continue;		////skip invalid cells

		VectorD pos=mac_grid.grid.Center(cell);
		VectorD vel=intp.Interpolate_Face_Vectors(velocity,pos);
		VectorD mid_pos=pos-vel*(real).5*dt;
		vel=intp.Interpolate_Face_Vectors(velocity,mid_pos);
		VectorD backtraced_pos=pos-vel*dt;
		T advected_density;
		if(use_zero_extrapolation&&(!ghost_mac_grid.grid.Inside(backtraced_pos)))advected_density*=(real)0;
		else advected_density=intp.Interpolate_Centers(ghost_density,backtraced_pos);
		density(cell)=advected_density;}
}

////this function only works for scalar advection
template<int d> void Semi_Lagrangian_Catmull_Rom(const real dt,const MacGrid<d>& ghost_grid,const FaceField<real,d>& ghost_velocity,const Field<real,d>& ghost_density,
	const MacGrid<d>& mac_grid,Field<real,d>& density)
{	Typedef_VectorDii(d);
	Interpolation<d> intp(ghost_grid);
	int cell_num=mac_grid.grid.Number_Of_Cells();
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){
		VectorDi cell=mac_grid.grid.Cell_Coord(i);
		VectorD pos=mac_grid.grid.Center(cell);
		VectorD vel=intp.Interpolate_Face_Vectors(ghost_velocity,pos);
		VectorD mid_pos=pos-vel*(real).5*dt;
		vel=intp.Interpolate_Face_Vectors(ghost_velocity,mid_pos);
		VectorD backtraced_pos=pos-vel*dt;
		real advected_density=intp.Interpolate_Centers_Catmull_Rom(ghost_density,backtraced_pos);
		density(cell)=advected_density;}
}

////advect a scalar or other type
template<class T,int d> void MacCormack(const real dt,const MacGrid<d>& mac_grid,const FaceField<real,d>& velocity,Field<T,d>& density)
{	Typedef_VectorDii(d);
	Field<T,d> ghost_density=density;
	Field<T,d> intermediate_density=density;
	Semi_Lagrangian<T,d>(dt,velocity,mac_grid,ghost_density,mac_grid,intermediate_density);
	Semi_Lagrangian<T,d>(-dt,velocity,mac_grid,intermediate_density,mac_grid,density);
	int cell_num=mac_grid.grid.Number_Of_Cells();
	#pragma omp parallel for
	for(int i=0;i<cell_num;i++){VectorDi cell=mac_grid.grid.Cell_Coord(i);
		density(cell)=intermediate_density(cell)+(real).5*(ghost_density(cell)-density(cell));}
}

////advect a vector field
template<int d> void MacCormack(const real dt,const MacGrid<d>& mac_grid,const FaceField<real,d>& velocity,FaceField<real,d>& vector,bool use_clamp=false)
{	Typedef_VectorDii(d);
	FaceField<real,d> ghost_vector=vector;
	FaceField<real,d> aux_vector=vector;

	
	Semi_Lagrangian<d>(dt,velocity,mac_grid,vector,mac_grid,aux_vector);
	Semi_Lagrangian<d>(-dt,velocity,mac_grid,aux_vector,mac_grid, vector);
	for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
		#pragma omp parallel for
		for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
			vector(axis,face)= aux_vector(axis,face)+(real).5*(ghost_vector(axis,face)- vector(axis,face));
	}}

	if (use_clamp) {
		Interpolation intp(mac_grid);
		for (int axis = 0; axis < d; axis++) {
			int face_num = mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for (int i = 0; i < face_num; i++) {
				VectorDi face = mac_grid.Face_Coord(axis, i);
				VectorD pos = mac_grid.Face_Center(axis, face);
				VectorD vel = intp.Interpolate_Face_Vectors(vector, pos);
				VectorD backtraced_pos = pos - vel * dt;
				VectorDi backtraced_cell = mac_grid.grid.Cell_Coord(backtraced_pos);
				VectorDi left_face = mac_grid.Cell_Left_Face(axis, backtraced_cell);
				VectorDi right_face = mac_grid.Cell_Right_Face(axis, backtraced_cell);
				real min = std::min(ghost_vector(axis, left_face), ghost_vector(axis, right_face));
				real max = std::max(ghost_vector(axis, left_face), ghost_vector(axis, right_face));
				if (vector(axis, face) < min || vector(axis, face) > max) {//check if overshoot
					vector(axis, face) = intp.Interpolate_Faces(ghost_vector, backtraced_pos, axis);
				}
			}
		}
	}
}
};

#endif

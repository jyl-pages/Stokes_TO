//////////////////////////////////////////////////////////////////////////
// Extrapolation
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __Extrapolation_h__
#define __Extrapolation_h__

#include "LevelSet.h"
#include "MacGrid.h"
#include "Field.h"
#include "FaceField.h"
#include "MacGrid.h"

namespace ExtpFunc
{
	////extrapolate face field with levelset
	//Fan:this has compiling issue on linux, unknown variable bc on line 38 and intp on line 43
	//template<int d> void Extrapolate_1(FaceField<real,d>& field_q,const LevelSet<d>& levelset,
	//	std::function<bool(const int,const Vector<int,d>&)>& Is_Fluid_Face,real narrow_band_width=(real)FLT_MAX)
	//{Typedef_VectorDii(d);
	//	MacGrid<d>& mac_grid=field_q.mac_grid;
	//	const VectorDi& cell_counts=mac_grid.grid.cell_counts;

	//	FaceField<real,d> ghost_field_q=field_q;

	//	////Fill ghost fields with fluid velocities
	//	iterate_face(axis,iter,mac_grid){const VectorDi& face=iter.Coord();
	//		if(Is_Fluid_Face(axis,face))continue;
	//		int n=0;real fluid_v=(real)0;
	//		for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
	//			VectorDi nb_face=Grid<d>::Nb_C(face,i);
	//			if(!mac_grid.face_grids[axis].Valid_Node(nb_face))continue;
	//			if(Is_Fluid_Face(axis,nb_face)){fluid_v+=ghost_field_q(axis,nb_face);n++;}}
	//		if(n>0)ghost_field_q(axis,face)=fluid_v/(real)n;}

	//	iterate_face(axis,iter,mac_grid){const VectorDi& face=iter.Coord();if(bc.Is_Psi_N(axis,face))continue;
	//		VectorD pos=mac_grid.Face_Center(axis,face);real phi=levelset.Phi(pos);
	//		if(phi>(real)0){
	//			if(phi<narrow_band_width){
	//				VectorD interface_pos=levelset.Closest_Point(pos);
	//				field_q(axis,face)=intp.Interpolate_Faces(ghost_field_q,interface_pos,axis);}
	//			else field_q(axis,face)=(real)0;}}
	//}

	////extrapolate field with levelset
	template<class T,int d> void Extrapolate(Field<T,d>& field_q,const LevelSet<d>& levelset,const real narrow_band_width=(real)FLT_MAX)
	{
		////TODO
	}

	////extrapolate face field with cell type
	template<int d> void Extrapolate(FaceField<real,d>& field_q,const Field<ushort,d>& type,const ushort fluid_type,const int narrow_band_width=4)
	{Typedef_VectorDii(d);
		VectorDi cell_counts=type.counts;
		MacGrid<d> mac_grid(cell_counts,(real)1);

		std::function<bool(const int,const VectorDi&)> Is_Fluid_Face=[=](const int axis,const VectorDi& face)->bool
		{
			VectorDi cell_0=face-VectorDi::Unit(axis);VectorDi cell_1=face;
			return (Grid<d>::Valid(cell_0,cell_counts)&&type(cell_0)==fluid_type)||(Grid<d>::Valid(cell_1,cell_counts)&&type(cell_1)==fluid_type);
		}; 

		FaceField<int,d> face_type(cell_counts,0);
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				if(Is_Fluid_Face(axis,face)){face_type(axis,face)=1;}}}

		FaceField<real,d> ghost_q=field_q;
		FaceField<int,d> ghost_face_type=face_type;
		for(int k=0;k<narrow_band_width;k++){
			for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
				//#pragma omp parallel for
				for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
					if(face_type(axis,face)==1){continue;}
					int n=0;real q=(real)0;
					for(int j=0;j<Grid<d>::Number_Of_Nb_C();j++){
						VectorDi nb_face=Grid<d>::Nb_C(face,j);
						if(!mac_grid.Valid_Face(axis,nb_face))continue;
						if(face_type(axis,nb_face)==1){q+=field_q(axis,nb_face);n++;}}
					if(n>0){ghost_q(axis,face)=q/(real)n;ghost_face_type(axis,face)=1;}}}
			field_q=ghost_q;face_type=ghost_face_type;}
	}

	////extrapolate field with cell type
	template<int d> void Extrapolate(Field<real,d>& field_q,Field<ushort,d>& type,const ushort fluid_type,const int narrow_band_width=4)
	{
		////TODO
	}
};


#endif
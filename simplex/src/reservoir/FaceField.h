//////////////////////////////////////////////////////////////////////////
// Face field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FaceField_h__
#define __FaceField_h__
#include <fstream>
#include "MacGrid.h"
#include "Field.h"

template<class T,int d> class FaceField
{Typedef_VectorDi(d);
public:
	MacGrid<d> mac_grid;					////ATTENTION: by default the mac_grid is NOT initialized so you may not read the grid information from it directly!
	ArrayF<Field<T,d>,d> face_fields;

	FaceField(const VectorDi& cell_counts=VectorDi::Zero()){Resize(cell_counts);}
	FaceField(const VectorDi& cell_counts,const T& value){Resize(cell_counts,value);}
	FaceField(const FaceField<T,d>& copy){*this=copy;}
	template<class T1> FaceField(const FaceField<T1, d>& copy) { *this = copy; }
	FaceField<T,d>& operator=(const FaceField<T,d>& copy){mac_grid=copy.mac_grid;face_fields=copy.face_fields;return *this;}
	template<class T1> FaceField<T, d>& operator=(const FaceField<T1, d>& copy) { mac_grid = copy.mac_grid; for (int i = 0; i < d; i++) face_fields[i] = copy.face_fields[i]; return *this; }

	void Resize(const VectorDi& cell_counts);
	void Resize(const VectorDi& cell_counts,const T& value);
	void Fill(const T& value){for(int i=0;i<d;i++)face_fields[i].Fill(value);}
	void Fill(const T& value,const int axis){face_fields[axis].Fill(value);}
	void Fill(const Vector<T,d>& value){for(int i=0;i<d;i++)Fill(value[i],i);}

    inline T& operator()(const int axis,const VectorDi& coord){return face_fields[axis](coord);}
    inline const T& operator()(const int axis,const VectorDi& coord) const {return face_fields[axis](coord);}
    inline T& operator()(const int axis,const int i){return face_fields[axis](i);}
    inline const T& operator()(const int axis,const int i) const {return face_fields[axis](i);}
    inline T& operator()(const int axis,const int i,const int j){return face_fields[axis](i,j);}
    inline const T& operator()(const int axis,const int i,const int j) const {return face_fields[axis](i,j);}
    inline T& operator()(const int axis,const int i,const int j,const int k){return face_fields[axis](i,j,k);}
    inline const T& operator()(const int axis,const int i,const int j,const int k) const {return face_fields[axis](i,j,k);}

	T Max();
	T Min();
	FaceField<T,d>& operator+=(const FaceField<T,d>& f2);
	FaceField<T,d>& operator+=(const Vector<T, d>& df);
	FaceField<T,d>& operator-=(const FaceField<T,d>& f2);
	FaceField<T,d>& operator*=(const FaceField<T,d>& f2);
	FaceField<T,d>& operator*=(const T& c);
	FaceField<T,d>& operator/=(const FaceField<T,d>& f2);

	////Help functions
	T Sum_Of_Incident_Faces(const VectorDi& cell)const;

	////IO
	virtual void Write_Binary(std::ostream& output) const;
	virtual void Write_Binary(const std::string& file_name)const;
	virtual void Read_Binary(std::istream& input);
	virtual void Read_Binary(const std::string& file_name);
	virtual void Write_To_File_3d(const std::string& file_name) const;
};

////Grid matrix mapping
template<int d> void Build_MacGrid_Face_Matrix_Bijective_Mapping(const MacGrid<d>& mac_grid, std::function<bool(const std::pair<int, int>)>& valid_face_func,
	/*rst*/FaceField<int, d>& grid_face_to_matrix,/*rst*/Array<std::pair<int, int>>& matrix_to_grid_face)
{
	Typedef_VectorDii(d);
	grid_face_to_matrix.Resize(mac_grid.grid.cell_counts, -1);
	matrix_to_grid_face.clear(); int c = 0;
	iterate_face(axis, iter, mac_grid) {
		const VectorDi& face = iter.Coord();
		int face_idx = mac_grid.face_grids[axis].Node_Index(face);
		if (valid_face_func(std::pair<int, int>(axis, face_idx))) { grid_face_to_matrix(axis, face) = c++; matrix_to_grid_face.push_back(std::make_pair(axis, face_idx)); }
	}
}

////Aux functions
////Dimension and type conversion
template<class T,int d> void Face_To_Cell_Conversion(const FaceField<T,d>& face_field,Field<Vector<T,d>,d>& cell_field);
template<class T,int d> void Face_To_Node_Conversion(const FaceField<T,d>& face_field,Field<Vector<T,d>,d>& node_field);

////IO functions
template<int d> void Write_Face_Vectors_To_File_3d_Fast(const FaceField<real,d>& face_field,const MacGrid<d>&,const std::string& file_name);
#endif

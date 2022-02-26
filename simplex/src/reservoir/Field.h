//////////////////////////////////////////////////////////////////////////
// Field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Field_h__
#define __Field_h__
#include <fstream>
#include <functional>
#include "Grid.h"

template<class T,int d> class Field
{Typedef_VectorDi(d);
public:
	VectorDi counts;
	Array<T> array;

	////Constructors
	Field(const VectorDi& _counts=VectorDi::Zero());
	Field(const VectorDi& _counts,const T value);
	Field(const VectorDi& _counts,const Array<T>& _array);
    Field(const Field<T,d>& copy);
	template<class T1> Field(const Field<T1, d>& copy) { *this = copy; }
    Field<T,d>& operator=(const Field<T,d>& copy);
	template<class T1> Field<T, d>& operator=(const Field<T1, d>& copy);

	////Resize and access
    void Resize(const VectorDi& _counts);
	void Resize(const VectorDi& _counts,const T value);
	void Resize(const int size);
	void Resize(const int size,const T value);
    
	inline T& operator()(const VectorDi& coord){return array[Index(coord)];} 
    inline const T& operator()(const VectorDi& coord) const {return array[Index(coord)];}
	inline T& operator()(const int i){return array[i];}
	inline const T& operator()(const int i) const {return array[i];}
	inline T& operator()(const int i,const int j){return array[i*counts[1]+j];} 
    inline const T& operator()(const int i,const int j) const {return array[i*counts[1]+j];}
	inline T& operator()(const int i,const int j,const int k){return array[i*counts[1]*counts[2]+j*counts[2]+k];}
	inline const T& operator()(const int i,const int j,const int k) const {return array[i*counts[1]*counts[2]+j*counts[2]+k];}

	////Operators
	Field<T,d>& operator+=(const Field<T,d>& f2);
	Field<T,d>& operator-=(const Field<T,d>& f2);
	Field<T,d>& operator*=(const Field<T,d>& f2);
	Field<T,d>& operator/=(const Field<T,d>& f2);

	////TO PARALLIZE
	template<class T_VAL2> Field<T,d>& operator+=(const T_VAL2& c) {for(auto& a:array)a+=c;return *this;}
	template<class T_VAL2> Field<T,d>& operator-=(const T_VAL2& c) {for(auto& a:array)a-=c;return *this;}
	template<class T_VAL2> Field<T,d>& operator*=(const T_VAL2& c) {for(auto& a:array)a*=c; return *this;}
	template<class T_VAL2> Field<T,d>& operator/=(const T_VAL2& c) {for(auto& a:array)a/=c;return *this;}
	
	void Set_Zero(){memset((void*)&array[0],0,sizeof(T)*array.size());}
    void Fill(const T& value){std::fill(array.begin(),array.end(),value);}
	bool Has(const T& value) const {auto find=std::find(array.begin(),array.end(),value);return find!=array.end();}

	////IO
	virtual void Write_Binary(std::ostream& output) const;
	virtual void Write_Binary(const std::string& file_name)const;
	virtual void Read_Binary(std::istream& input);
	virtual void Read_Binary(const std::string& file_name);
	virtual void Write_To_File_3d(const std::string& file_name) const;	////works for scalar field only if dim conversion is needed

	////Helper functions
protected:
	inline int Index(const Vector1i& c) const {return c[0];}
	inline int Index(const Vector2i& c) const {return c[0]*counts[1]+c[1];}
	inline int Index(const Vector3i& c) const {return c[0]*counts[1]*counts[2]+c[1]*counts[2]+c[2];}
	inline int Index(const int i) const {return i;}
	inline int Index(const int i,const int j) const {return i*counts[1]+j;}
	inline int Index(const int i,const int j,const int k) const {return i*counts[1]*counts[2]+j*counts[2]+k;}
};

template<class T, int d>
template<class T1>
Field<T, d>& Field<T, d>::operator = (const Field<T1, d>& copy) {
	counts = copy.counts; Resize(counts);
#pragma omp parallel for
	for (int i = 0; i < array.size(); i++) array[i] = (T)copy.array[i];
	return *this;
}

////Aux functions
////linear algebra operations
template<class T,int d> void Axpy(T a,const Field<T,d>& x,const Field<T,d>& y,Field<T,d>& axpy);			////axpy=a*x+y
template<class T,int d> T Dot(const Field<T,d>& x,const Field<T,d>& y);										////dot=x dot y
////reduction operations
template<class T,int d> T Sum(const Field<T,d>& x);
template<class T,int d> T Max(const Field<T,d>& x);
template<class T,int d> T Min(const Field<T,d>& x); 	

////Dimension conversion
template<class T,int d1,int d2> void Dim_Conversion(const Field<T,d1>& input,Field<T,d2>& output);
template<class T,int d1,int d2> void VF_Dim_Conversion(const Field<Vector<T,d1>,d1>& input,Field<Vector<T,d2>,d2>& output);
template<class T,int d1,int d2> void TF_Dim_Conversion(const Field<Matrix<T,d1>,d1>& input,Field<Matrix<T,d2>,d2>& output);

////Flood fill
template<int d> void Flood_Fill(Field<int,d>& psi_D,const Grid<d>& grid,const int boundary,const int interior,const int exterior,const Vector<int,d> seed=Vector<int,d>::Ones()*-1);

////Grid matrix mapping
template<int d> void Build_Grid_Node_With_Incident_Cell_Matrix_Bijective_Mapping(const Grid<d>& grid,std::function<bool(const int)>& valid_cell_func,
		/*rst*/Field<int,d>& grid_node_to_matrix,/*rst*/Array<int>& matrix_to_grid_node);
template<int d> void Build_Grid_Node_Matrix_Bijective_Mapping(const Grid<d>& grid,std::function<bool(const int)>& valid_node_func,
		/*rst*/Field<int,d>& grid_node_to_matrix,/*rst*/Array<int>& matrix_to_grid_node);
template<int d> void Build_Grid_Cell_Matrix_Bijective_Mapping(const Grid<d>& grid,std::function<bool(const int)>& valid_cell_func,
		/*rst*/Field<int,d>& grid_cell_to_matrix,/*rst*/Array<int>& matrix_to_grid_cell);

////Type conversion
////Convert a serialized-scalar array (e.g., VectorX,VectorN<real>,Array<real>) to a vector array
template<class T_ARRAY,class T,int d> void Convert_To_Field(const T_ARRAY& input,Field<Vector<T,d>,d>& output);		////Assuming the field is initialized
template<class T_ARRAY,class T,int d> void Convert_To_Field1(const T_ARRAY& input,Field<Vector<T,d>,1>& output);
////Convert a serialized-vector array to a vector array (the vectors of the two arrays are with different dimensions)
template<class T_ARRAY,class T,int d1,int d2> void Convert_To_Field(const T_ARRAY& input,Field<Vector<T,d2>,d2>& output);		////Assuming the field is initialized
template<class T_ARRAY,class T,int d1,int d2> void Convert_To_Field1(const T_ARRAY& input,Field<Vector<T,d2>,1>& output);

////IO helpers
////Write a scalar array that represents a vector (2d or 3d) array into a Field file, 2D vectors will be converted to 3D first
template<class T_ARRAY,class T,int d> void Write_To_Field1_3d(const T_ARRAY& input,const std::string& file_name);
template<class T_ARRAY,class T,int d> void Write_To_Field_3d(const T_ARRAY& input,const Vector<int,d>& counts,const std::string& file_name);

#endif

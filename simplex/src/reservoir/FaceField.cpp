//////////////////////////////////////////////////////////////////////////
// Face Field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <limits>
#include "Interpolation.h"
#include "FaceField.h"
#include "File.h"

template<class T,int d> void FaceField<T,d>::Resize(const VectorDi& cell_counts)
{
	if(mac_grid.grid.cell_counts==cell_counts)return;
	mac_grid.Initialize(cell_counts);
	for(int i=0;i<d;i++)face_fields[i].Resize(mac_grid.face_counts[i]);
}

template<class T,int d> void FaceField<T,d>::Resize(const VectorDi& cell_counts,const T& value)
{
	if(mac_grid.grid.cell_counts==cell_counts)return;
	mac_grid.Initialize(cell_counts);
	for(int i=0;i<d;i++)face_fields[i].Resize(mac_grid.face_counts[i],value);
}

template<class T,int d> T FaceField<T,d>::Max()
{
	////TOACC
	T v_max=-std::numeric_limits<T>::max();
	for(auto i=0;i<d;i++){v_max=std::max(v_max,*std::max_element(face_fields[i].array.begin(),face_fields[i].array.end()));}return v_max;
}

template<class T,int d> T FaceField<T,d>::Min()
{
	////TOACC
	T v_min=std::numeric_limits<T>::max();
	for(auto i=0;i<d;i++){v_min=std::min(v_min,*std::min_element(face_fields[i].array.begin(),face_fields[i].array.end()));}return v_min;
}

template<class T,int d> FaceField<T,d>& FaceField<T,d>::operator+=(const FaceField<T,d>& f2)
{
	for(auto i=0;i<d;i++){face_fields[i]+=f2.face_fields[i];}return *this;
}

template<class T, int d>
FaceField<T, d>& FaceField<T, d>::operator+=(const Vector<T, d>& df)
{
	for (auto i = 0; i < d; i++) face_fields[i] += df[i]; return *this;
}

template<class T,int d> FaceField<T,d>& FaceField<T,d>::operator-=(const FaceField<T,d>& f2)
{
	for(auto i=0;i<d;i++){face_fields[i]-=f2.face_fields[i];}return *this;
}

template<class T,int d> FaceField<T,d>& FaceField<T,d>::operator*=(const FaceField<T,d>& f2)
{
	for(auto i=0;i<d;i++){face_fields[i]*=f2.face_fields[i];}return *this;
}

template<class T,int d> FaceField<T,d>& FaceField<T,d>::operator*=(const T& c)
{
	for(auto i=0;i<d;i++){face_fields[i]*=c;}return *this;
}

template<class T,int d> FaceField<T,d>& FaceField<T,d>::operator/=(const FaceField<T,d>& f2)
{
	for(auto i=0;i<d;i++){face_fields[i]/=f2.face_fields[i];}return *this;
}

template<class T, int d> T FaceField<T, d>::Sum_Of_Incident_Faces(const VectorDi& cell)const {
	int axis; VectorDi face;
	T sum = (T)0;
	for (int i = 0; i < MacGrid<d>::Number_Of_Cell_Incident_Faces(); i++) {
		MacGrid<d>::Cell_Incident_Face(cell, i, axis, face);
		sum += (*this)(axis, face);
	}
	return sum;
}

template<class T,int d> void FaceField<T, d>::Write_Binary(std::ostream& output) const
{
	File::Write_Binary(output,mac_grid.grid.cell_counts);
	for(int i=0;i<d;i++)File::Write_Binary_Array(output,&face_fields[i].array[0],(int)face_fields[i].array.size());
}

template<class T, int d> void FaceField<T, d>::Write_Binary(const std::string& file_name) const
{
	std::ofstream output(file_name, std::ios::binary);
	if (!output) { std::cerr << "FaceField<T, d>::Write_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Write_Binary(output);
	output.close();
}

template<class T,int d> void FaceField<T,d>::Read_Binary(std::istream& input)
{
	VectorDi cell_counts;File::Read_Binary(input,cell_counts);Resize(cell_counts);
    for(int i=0;i<d;i++)File::Read_Binary_Array(input,&face_fields[i].array[0],(int)face_fields[i].array.size());
}

template<class T, int d> void FaceField<T, d>::Read_Binary(const std::string& file_name)
{
	std::ifstream input(file_name, std::ios::binary);
	if (!input) { std::cerr << "FaceField<T, d>::Read_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Read_Binary(input);
	input.close();
}

template<class T,int d> void FaceField<T,d>::Write_To_File_3d(const std::string& file_name) const
{
	if constexpr (d==3){
		Field<Vector<T,d>,d> v;Face_To_Cell_Conversion(*this,v);
		File::Write_Binary_To_File(file_name,v);}
	else{
		Field<Vector<T,d>,d> v;Face_To_Cell_Conversion(*this,v);
		Field<Vector<T,3>,3> v3;VF_Dim_Conversion<T,d,3>(v,v3);
		File::Write_Binary_To_File(file_name,v3);}
}

template class FaceField<ushort,2>;
template class FaceField<ushort,3>;
template class FaceField<int,2>;
template class FaceField<int,3>;
template class FaceField<float,2>;
template class FaceField<float,3>;
template class FaceField<double,2>;
template class FaceField<double,3>;

template<class T,int d> void Face_To_Cell_Conversion(const FaceField<T,d>& face_field,Field<Vector<T,d>,d>& cell_field)
{
	Interpolation<d> intp(face_field.mac_grid);cell_field.Resize(face_field.mac_grid.grid.cell_counts);
	intp.Interpolate_Faces_To_Cells(face_field,cell_field);
}

template<class T,int d> void Face_To_Node_Conversion(const FaceField<T,d>& face_field,Field<Vector<T,d>,d>& node_field)
{
	Interpolation<d> intp(face_field.mac_grid);node_field.Resize(face_field.mac_grid.grid.node_counts);
	intp.Interpolate_Faces_To_Nodes(face_field,node_field);
}

#define Inst_Helper(real,d) \
template void Face_To_Cell_Conversion<real,d>(const FaceField<real,d>&,Field<Vector<real,d>,d>&);\
template void Face_To_Node_Conversion<real,d>(const FaceField<real,d>&,Field<Vector<real,d>,d>&)
Inst_Helper(ushort,2);Inst_Helper(ushort,3);
Inst_Helper(int,2);Inst_Helper(int,3);
Inst_Helper(float,2);Inst_Helper(float,3);
Inst_Helper(double,2);Inst_Helper(double,3);
#undef Inst_Helper


template<int d> void Write_Face_Vectors_To_File_3d_Fast(const FaceField<real,d>& face_field,const MacGrid<d>& mac_grid,const std::string& file_name)
{Typedef_VectorDii(d);
	int n=(int)mac_grid.Number_Of_Faces();
	float* xf=new float[n*8];
	memset((void*)xf,0,sizeof(float)*n*8);

	int s=0;for(int axis=0;axis<d;axis++){
		int face_num=mac_grid.Number_Of_Faces(axis);
		if(axis>0)s+=mac_grid.Number_Of_Faces(axis-1);
		#pragma omp parallel for
		for(int f=0;f<face_num;f++){VectorDi face=mac_grid.Face_Coord(axis,f);int i=s+f;
			VectorD pos=mac_grid.Face_Center(axis,face);
			VectorD pos2=pos+VectorD::Unit(axis)*face_field(axis,face);
			if constexpr (d==2){
				xf[i*8]=(float)pos[0];
				xf[i*8+1]=(float)pos[1];

				xf[i*8+4]=(float)pos2[0];
				xf[i*8+5]=(float)pos2[1];}
			else if constexpr (d==3){
				xf[i*8]=(float)pos[0];
				xf[i*8+1]=(float)pos[1];
				xf[i*8+2]=(float)pos[2];

				xf[i*8+4]=(float)pos2[0];
				xf[i*8+5]=(float)pos2[1];
				xf[i*8+6]=(float)pos2[2];}}}

	std::ofstream output(file_name,std::ios::binary);if(!output)return;
	File::Write_Binary(output,n*8);
	File::Write_Binary_Array(output,xf,n*8);
	delete [] xf;
}

template void Write_Face_Vectors_To_File_3d_Fast<2>(const FaceField<real,2>&,const MacGrid<2>&,const std::string&);
template void Write_Face_Vectors_To_File_3d_Fast<3>(const FaceField<real,3>&,const MacGrid<3>&,const std::string&);

//////////////////////////////////////////////////////////////////////////
// Field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <numeric>
#include <stack>
#include <string>
#include "Field.h"
#include "File.h"
#include "AuxFunc.h"
#include "SparseFunc.h"

template<class T,int d> Field<T,d>::Field(const VectorDi& _counts/*=VectorDi::Zero()*/):counts(_counts){Resize(counts);}

template<class T,int d> Field<T,d>::Field(const VectorDi& _counts,const T value):counts(_counts){Resize(counts,value);}

template<class T,int d> Field<T,d>::Field(const VectorDi& _counts,const Array<T>& _array):counts(_counts)
{Resize(counts);std::memcpy((void*)&array[0],(const void*)&_array[0],sizeof(T)*array.size());}	

template<class T,int d> Field<T,d>::Field(const Field<T,d>& copy){*this=copy;}

template<class T,int d> void Field<T,d>::Resize(const VectorDi& _counts)
{counts=_counts;int n=counts.prod();if(n>0&&n!=array.size())array.resize(n);}

template<class T,int d> void Field<T,d>::Resize(const VectorDi& _counts,const T value)
{counts=_counts;int n=counts.prod();if(n>0){array.resize(n);std::fill(array.begin(),array.end(),value);}}

template<class T,int d> void Field<T,d>::Resize(const int size)
{counts=VectorDi::Ones();counts[0]=size;Resize(counts);}

template<class T,int d> void Field<T,d>::Resize(const int size,const T value)
{counts=VectorDi::Ones();counts[0]=size;Resize(counts,value);}

template<class T,int d> Field<T,d>& Field<T,d>::operator=(const Field<T,d>& copy)
{counts=copy.counts;Resize(counts);std::memcpy((void*)&array[0],(const void*)&copy.array[0],sizeof(T)*array.size());return *this;}

template<class T,int d> Field<T,d>& Field<T,d>::operator+=(const Field<T,d>& f2)
{
	int n=(int)array.size();
	#pragma omp parallel for
	for(int i=0;i<n;i++)array[i]+=f2.array[i];return *this;
}

template<class T,int d> Field<T,d>& Field<T,d>::operator-=(const Field<T,d>& f2)
{
	int n=(int)array.size();
	#pragma omp parallel for
	for(int i=0;i<n;i++)array[i]-=f2.array[i];return *this;
}

template<class T,int d> Field<T,d>& Field<T,d>::operator*=(const Field<T,d>& f2)
{std::cerr<<"Field operator *= not implemented"<<std::endl;return *this;}

template<class T,int d> Field<T,d>& Field<T,d>::operator/=(const Field<T,d>& f2)
{std::cerr<<"Field operator /= not implemented"<<std::endl;return *this;}

#define Inst_Helper(T,d) \
template<> Field<T,d>& Field<T,d>::operator*=(const Field<T,d>& f2) \
{ \
	int n=(int)array.size(); \
	for(int i=0;i<n;i++)array[i]*=f2.array[i];return *this; \
} \
template<> Field<T,d>& Field<T,d>::operator/=(const Field<T,d>& f2) \
{ \
	int n=(int)array.size(); \
	for(int i=0;i<n;i++)array[i]/=f2.array[i];return *this; \
}
Inst_Helper(float,1);Inst_Helper(float,2);Inst_Helper(float,3);
Inst_Helper(double,1);Inst_Helper(double,2);Inst_Helper(double,3);
Inst_Helper(int,1);Inst_Helper(int,2);Inst_Helper(int,3);
Inst_Helper(ushort,1);Inst_Helper(ushort,2);Inst_Helper(ushort,3);
#undef Inst_Helper

template<class T,int d> void Field<T,d>::Write_Binary(std::ostream &output) const
{
    File::Write_Binary(output,counts);
    File::Write_Binary_Array(output,&array[0],(int)array.size());
}

template<class T, int d>
void Field<T, d>::Write_Binary(const std::string& file_name) const
{
	std::ofstream output(file_name, std::ios::binary);
	if (!output) { std::cerr << "Field<T, d>::Write_Binary error: cannot open file " << file_name << "\n "; exit(1); }
	Write_Binary(output);
	output.close();
}

template<class T,int d> void Field<T,d>::Read_Binary(std::istream &input)
{
    File::Read_Binary(input,counts);Resize(counts);
    File::Read_Binary_Array(input,&array[0],(int)array.size());
}

template<class T, int d>
void Field<T, d>::Read_Binary(const std::string& file_name)
{
	std::ifstream input(file_name, std::ios::binary);
	if (!input) { std::cerr << "Field<T, d>::Read_Binary error: cannot open file " << file_name << "\n "; exit(0); }
	Read_Binary(input);
	input.close();
}

template<class T,int d> void Field<T,d>::Write_To_File_3d(const std::string& file_name) const
{	
	if constexpr(d==1||d==3)
		File::Write_Binary_To_File(file_name,*this);
	else{
		////PBG: needs conversion when d==2, works for scalar field only
		Field<T,3> field3;Dim_Conversion<T,d,3>(*this,field3);
		File::Write_Binary_To_File(file_name,field3);}
}

template<> void Field<Vector2,2>::Write_To_File_3d(const std::string& file_name) const
{
	Field<Vector3,3> field3;VF_Dim_Conversion<real,2,3>(*this,field3);
	File::Write_Binary_To_File(file_name,field3);
}

template<> void Field<Vector2i,2>::Write_To_File_3d(const std::string& file_name) const
{
	Field<Vector3i,3> field3;VF_Dim_Conversion<int,2,3>(*this,field3);
	File::Write_Binary_To_File(file_name,field3);
}

template<> void Field<Matrix2,2>::Write_To_File_3d(const std::string& file_name) const
{
	Field<Matrix3,3> field3;TF_Dim_Conversion<real,2,3>(*this,field3);
	File::Write_Binary_To_File(file_name,field3);
}

template<> void Field<Matrix2i,2>::Write_To_File_3d(const std::string& file_name) const
{
	Field<Matrix3i,3> field3;TF_Dim_Conversion<int,2,3>(*this,field3);
	File::Write_Binary_To_File(file_name,field3);
}

template<> void Field<Matrix2,1>::Write_To_File_3d(const std::string& file_name) const
{
	Field<Matrix3,1> field3;field3.Resize(counts);
	AuxFunc::Dim_Conversion_Array<real,2,3>(array,field3.array);
	File::Write_Binary_To_File(file_name,field3);
}

template<> void Field<Matrix2i,1>::Write_To_File_3d(const std::string& file_name) const
{
	Field<Matrix3i,1> field3;field3.Resize(counts);
	AuxFunc::Dim_Conversion_Array<int,2,3>(array,field3.array);
	File::Write_Binary_To_File(file_name,field3);
}

template class Field<double,1>;
template class Field<double,2>;
template class Field<double,3>;
template class Field<float,1>;
template class Field<float,2>;
template class Field<float,3>;
template class Field<int,1>;
template class Field<int,2>;
template class Field<int,3>;
template class Field<short,1>;
template class Field<short,2>;
template class Field<short,3>;
template class Field<ushort,1>;
template class Field<ushort,2>;
template class Field<ushort,3>;

template class Field<Vector2d,1>;
template class Field<Vector3d,1>;

template class Field<Vector2d,2>;
template class Field<Vector3d,2>;
template class Field<Vector4d,2>;
template class Field<Vector5d,2>;
template class Field<Vector6d,2>;

template class Field<Vector2d,3>;
template class Field<Vector3d,3>;
template class Field<Vector4d,3>;
template class Field<Vector5d,3>;
template class Field<Vector6d,3>;

template class Field<Matrix2d,1>;
template class Field<Matrix3d,1>;
template class Field<Matrix2d,2>;
template class Field<Matrix3d,3>;

template class Field<Vector2f,1>;
template class Field<Vector3f,1>;

template class Field<Vector2f,2>;
template class Field<Vector3f,2>;
template class Field<Vector4f,2>;
template class Field<Vector5f,2>;
template class Field<Vector6f,2>;

template class Field<Vector2f,3>;
template class Field<Vector3f,3>;
template class Field<Vector4f,3>;
template class Field<Vector5f,3>;
template class Field<Vector6f,3>;

template class Field<Matrix2f,1>;
template class Field<Matrix3f,1>;
template class Field<Matrix2f,2>;
template class Field<Matrix3f,3>;
template class Field<Vector2i,1>;
template class Field<Vector3i,1>;
template class Field<Vector2i,2>;
template class Field<Vector3i,3>;
template class Field<Matrix2i,1>;
template class Field<Matrix3i,1>;
template class Field<Matrix2i,2>;
template class Field<Matrix3i,3>;
template class Field<Vector<ushort,2>,2>;
template class Field<Vector<ushort,2>,3>;
template class Field<Vector<ushort,3>,2>;
template class Field<Vector<ushort,3>,3>;

template class Field<C,2>;
template class Field<C,3>;
template class Field<Vector<C,2>,2>;
template class Field<Vector<C,2>,3>;
template class Field<Vector<C,3>,2>;
template class Field<Vector<C,3>,3>;

//////////////////////////////////////////////////////////////////////////
////Auxiliary functions: linear algebra

template<class T,int d> void Axpy(T a,const Field<T,d>& x,const Field<T,d>& y,Field<T,d>& axpy)
{
	int n=(int)axpy.array.size();
	#pragma omp parallel for
	for(int i=0;i<n;i++)axpy.array[i]=a*x.array[i]+y.array[i];	
}

template<class T,int d> T Dot(const Field<T,d>& x,const Field<T,d>& y)
{
	int n=(int)x.array.size();
	T dot=(T)0;for(int i=0;i<n;i++)dot+=x.array[i]*y.array[i];return dot;
}

template<class T,int d> T Sum(const Field<T,d>& x)
{
	return std::accumulate(x.array.begin(),x.array.end(),Zero<T>());
}

template<class T,int d> T Max(const Field<T,d>& x)
{
	return *std::max_element(x.array.begin(),x.array.end());
}

template<class T,int d> T Min(const Field<T,d>& x)
{
	return *std::min_element(x.array.begin(),x.array.end());
}

#define Inst_Helper(T,d) \
template void Axpy<T,d>(T,const Field<T,d>&,const Field<T,d>&,Field<T,d>&); \
template T Dot<T,d>(const Field<T,d>&,const Field<T,d>&); \
template T Sum<T,d>(const Field<T,d>&); \
template T Max<T,d>(const Field<T,d>&); \
template T Min<T,d>(const Field<T,d>&);
Inst_Helper(double,2);Inst_Helper(double,3);
Inst_Helper(float,2);Inst_Helper(float,3);
#undef Inst_Helper

//////////////////////////////////////////////////////////////////////////
////Auxiliary functions: dimension conversion
template<class T,int d1,int d2> void Dim_Conversion(const Field<T,d1>& input,Field<T,d2>& output)
{
    Vector<int,d2> c2;AuxFunc::Dim_Conversion<int,d1,d2>(input.counts,c2,1);
    output.Resize(c2);std::copy(input.array.begin(),input.array.end(),output.array.begin());
}

template<class T,int d1,int d2> void VF_Dim_Conversion(const Field<Vector<T,d1>,d1>& input,Field<Vector<T,d2>,d2>& output)
{
	Vector<int,d2> c2;AuxFunc::Dim_Conversion<int,d1,d2>(input.counts,c2,1);
	output.Resize(c2);AuxFunc::Dim_Conversion_Array<T,d1,d2>(input.array,output.array,(T)0);
}

template<class T,int d1,int d2> void TF_Dim_Conversion(const Field<Matrix<T,d1>,d1>& input, Field<Matrix<T,d2>,d2>& output)
{
	Vector<int,d2> c2;AuxFunc::Dim_Conversion<int,d1,d2>(input.counts,c2,1);
	output.Resize(c2);AuxFunc::Dim_Conversion_Array<T,d1,d2>(input.array,output.array);
}

template<int d> void Flood_Fill(Field<int,d>& psi_D,const Grid<d>& grid,const int boundary,
	const int interior,const int exterior,const Vector<int,d> seed/*=Vector<int,d>::Ones()*-1*/)
{	Typedef_VectorDii(d);
	Field<int,d> color=psi_D;
	std::stack<int> cell_index_stack;

	if(seed==VectorDi::Ones()*-1){for(int i=0;i<psi_D.array.size();i++)if(psi_D.array[i]!=boundary){cell_index_stack.push(i);break;}}
	else{cell_index_stack.push(grid.Cell_Index(seed));}

	while(!cell_index_stack.empty()){
		int cell_index=cell_index_stack.top();cell_index_stack.pop();
		VectorDi cell=grid.Cell_Coord(cell_index);color(cell)=exterior;
		for(int i=0;i<grid.Number_Of_Nb_C();i++){VectorDi nb_cell=grid.Nb_C(cell,i);
			if(grid.Valid_Cell(nb_cell)&&color(nb_cell)!=exterior&&color(nb_cell)!=boundary){
				cell_index_stack.push(grid.Cell_Index(nb_cell));}}}

	#pragma omp parallel for
	for(auto i=0;i<color.array.size();i++){
		if(color.array[i]==exterior)psi_D.array[i]=exterior;
		if(color.array[i]!=exterior&&color.array[i]!=boundary)psi_D.array[i]=interior;}
}

template<int d> void Build_Grid_Node_With_Incident_Cell_Matrix_Bijective_Mapping(
	const Grid<d>& grid,std::function<bool(const int)>& valid_cell_func, 
	/*rst*/Field<int, d>& grid_node_to_matrix,/*rst*/Array<int>& matrix_to_grid_node)
{	Typedef_VectorDii(d);
	grid_node_to_matrix.Resize(grid.node_counts);grid_node_to_matrix.Fill(-1);
	matrix_to_grid_node.clear();int c=0;
	iterate_node(iter,grid){const VectorDi& node=iter.Coord();bool valid_node=false;
		for(int i=0;i<Grid<d>::Number_Of_Node_Incident_Cells();i++){
			VectorDi cell=grid.Node_Incident_Cell(node,i);
			if(grid.Valid_Cell(cell)&&valid_cell_func(grid.Cell_Index(cell))){valid_node=true;break;}}
		if(valid_node){grid_node_to_matrix(node)=c++;matrix_to_grid_node.push_back(grid.Node_Index(node));}}
}

template<int d> void Build_Grid_Cell_Matrix_Bijective_Mapping(
	const Grid<d>& grid,std::function<bool(const int)>& valid_cell_func,
	/*rst*/Field<int,d>& grid_cell_to_matrix,/*rst*/Array<int>& matrix_to_grid_cell)
{	Typedef_VectorDii(d);
	grid_cell_to_matrix.Resize(grid.cell_counts,-1);
	matrix_to_grid_cell.clear();int c=0;
	iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
		if(valid_cell_func(grid.Cell_Index(cell))){grid_cell_to_matrix(cell)=c++;matrix_to_grid_cell.push_back(grid.Cell_Index(cell));}}
}

template<int d> void Build_Grid_Node_Matrix_Bijective_Mapping(
	const Grid<d>& grid,std::function<bool(const int)>& valid_node_func, 
	/*rst*/Field<int,d>& grid_node_to_matrix,/*rst*/Array<int>& matrix_to_grid_node)
{	Typedef_VectorDii(d);
	grid_node_to_matrix.Resize(grid.node_counts,-1);
	matrix_to_grid_node.clear();int c=0;
	iterate_node(iter,grid){const VectorDi& node=iter.Coord();
		if(valid_node_func(grid.Node_Index(node))){grid_node_to_matrix(node)=c++;matrix_to_grid_node.push_back(grid.Node_Index(node));}}
}

#define Inst_Helper(T,d1,d2) \
template void Dim_Conversion<T,d1,d2>(const Field<T,d1>&,Field<T,d2>&); \
template void VF_Dim_Conversion<T,d1,d2>(const Field<Vector<T,d1>,d1>&,Field<Vector<T,d2>,d2>&); \
template void TF_Dim_Conversion<T,d1,d2>(const Field<Matrix<T,d1>,d1>&,Field<Matrix<T,d2>,d2>&);
Inst_Helper(double,3,2);Inst_Helper(double,2,3);Inst_Helper(double,2,2);Inst_Helper(double,3,3);
Inst_Helper(float,3,2);Inst_Helper(float,2,3);Inst_Helper(float,2,2);Inst_Helper(float,3,3);
Inst_Helper(int,3,2);Inst_Helper(int,2,3);Inst_Helper(int,2,2);Inst_Helper(int,3,3);
Inst_Helper(ushort,3,2);Inst_Helper(ushort,2,3);Inst_Helper(ushort,2,2);Inst_Helper(ushort,3,3);
#undef Inst_Helper

#define Inst_Helper(d) \
template void Flood_Fill<d>(Field<int,d>&,const Grid<d>&,const int,const int,const int,const Vector<int,d>);
Inst_Helper(2);Inst_Helper(3);
#undef Inst_Helper

#define Inst_Helper(d) \
template void Build_Grid_Node_With_Incident_Cell_Matrix_Bijective_Mapping<d>(const Grid<d>&,std::function<bool(const int)>&,/*rst*/Field<int,d>&,/*rst*/Array<int>&); \
template void Build_Grid_Node_Matrix_Bijective_Mapping<d>(const Grid<d>&,std::function<bool(const int)>&,/*rst*/Field<int,d>&,/*rst*/Array<int>&); \
template void Build_Grid_Cell_Matrix_Bijective_Mapping<d>(const Grid<d>&,std::function<bool(const int)>&,/*rst*/Field<int,d>&,/*rst*/Array<int>&); 
Inst_Helper(2);Inst_Helper(3);
#undef Inst_Helper

template<class T_ARRAY,class T,int d> void Convert_To_Field(const T_ARRAY& u,Field<Vector<T,d>,d>& output)
{
	const int n=(int)output.array.size();
	for(int i=0;i<n;i++){Vector<T,d> v;for(int j=0;j<d;j++)v[j]=u[i*d+j];output.array[i]=v;}
}
#define Inst_Helper(T_ARRAY,T,d) \
template void Convert_To_Field<T_ARRAY,T,d>(const T_ARRAY&,Field<Vector<T,d>,d>&);
Inst_Helper(VectorN<real>,real,2);Inst_Helper(VectorN<real>,real,3);
Inst_Helper(Array<real>,real,2);Inst_Helper(Array<real>,real,3);
#undef Inst_Helper

template<class T_ARRAY,class T,int d> void Convert_To_Field1(const T_ARRAY& u,Field<Vector<T,d>,1>& output)
{
	const int n=(int)u.size()/d;output.Resize(n);
	for(int i=0;i<n;i++){Vector<T,d> v;for(int j=0;j<d;j++)v[j]=u[i*d+j];output.array[i]=v;}
}
#define Inst_Helper(T_ARRAY,T,d) \
template void Convert_To_Field1<T_ARRAY,T,d>(const T_ARRAY&,Field<Vector<T,d>,1>&);
Inst_Helper(VectorN<real>,real,2);Inst_Helper(VectorN<real>,real,3);
Inst_Helper(Array<real>,real,2);Inst_Helper(Array<real>,real,3);
#undef Inst_Helper

template<class T_ARRAY,class T,int d1,int d2> void Convert_To_Field(const T_ARRAY& u,Field<Vector<T,d2>,d2>& output)
{
	const int n=(int)output.array.size();
	for(int i=0;i<n;i++){Vector<T,d2> v=Vector<T,d2>::Zero();for(int j=0;j<d1;j++)v[j]=u[i*d1+j];output.array[i]=v;}
}
#define Inst_Helper(T_ARRAY,T,d1,d2) \
template void Convert_To_Field<T_ARRAY,T,d1,d2>(const T_ARRAY&,Field<Vector<T,d2>,d2>&);
Inst_Helper(VectorN<real>,real,2,3);
Inst_Helper(Array<real>,real,2,3);
#undef Inst_Helper

template<class T_ARRAY,class T,int d1,int d2> void Convert_To_Field1(const T_ARRAY& u,Field<Vector<T,d2>,1>& output)
{
	const int n=(int)u.size()/d1;output.Resize(n);
	for(int i=0;i<n;i++){Vector<T,d2> v=Vector<T,d2>::Zero();for(int j=0;j<d1;j++)v[j]=u[i*d1+j];output.array[i]=v;}
}
#define Inst_Helper(T_ARRAY,T,d1,d2) \
template void Convert_To_Field1<T_ARRAY,T,d1,d2>(const T_ARRAY&,Field<Vector<T,d2>,1>&);
Inst_Helper(VectorN<real>,real,2,3);
Inst_Helper(Array<real>,real,2,3);
#undef Inst_Helper

template<class T_ARRAY,class T,int d> void Write_To_Field1_3d(const T_ARRAY& input,const std::string& file_name)
{
	Field<Vector<T,3>,1> output;
	Convert_To_Field1<T_ARRAY,T,d,3>(input,output);
	File::Write_Binary_To_File(file_name,output);
}
#define Inst_Helper(T_ARRAY,T,d) \
template void Write_To_Field1_3d<T_ARRAY,T,d>(const T_ARRAY&,const std::string& file_name);
Inst_Helper(VectorN<real>,real,2);Inst_Helper(VectorN<real>,real,3);
Inst_Helper(Array<real>,real,2);Inst_Helper(Array<real>,real,3);
#undef Inst_Helper

template<class T_ARRAY,class T,int d> void Write_To_Field_3d(const T_ARRAY& input,const Vector<int,d>& counts,const std::string& file_name)
{
	Vector3i counts_3d;AuxFunc::Dim_Conversion<int,d,3>(counts,counts_3d,1);
	Field<Vector<T,3>,3> output;output.Resize(counts_3d);
	Convert_To_Field<T_ARRAY,real,d,3>(input,output);
	File::Write_Binary_To_File(file_name,output);
}
#define Inst_Helper(T_ARRAY,T,d) \
template void Write_To_Field_3d<T_ARRAY,T,d>(const T_ARRAY&,const Vector<int,d>&,const std::string& file_name);
Inst_Helper(VectorN<real>,real,2);Inst_Helper(VectorN<real>,real,3);
Inst_Helper(Array<real>,real,2);Inst_Helper(Array<real>,real,3);
#undef Inst_Helper



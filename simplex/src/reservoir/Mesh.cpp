//////////////////////////////////////////////////////////////////////////
// Mesh
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Mesh.h"
#include "Hashtable.h"
#include "MeshFunc.h"
#include "File.h"
#include "AuxFunc.h"

////SimplicialMesh
template<int d,int e_d> SimplicialMesh<d,e_d>::SimplicialMesh(const ArrayPtr<VectorD> _vertices)
{if(_vertices==nullptr)vertices=std::make_shared<Array<VectorD> >();else vertices=_vertices;}

template<int d,int e_d> SimplicialMesh<d,e_d>& SimplicialMesh<d,e_d>::operator=(const SimplicialMesh<d,e_d>& copy)	////deep copy
{if(vertices==nullptr)vertices=std::make_shared<Array<VectorD> >();*vertices=*(copy.vertices);elements=copy.elements;return *this;}

template<int d,int e_d> void SimplicialMesh<d,e_d>::Write_Binary(std::ostream&output) const
{
	int vtx_n=(int)(*vertices).size();
	File::Write_Binary(output,vtx_n);
	File::Write_Binary_Array(output,&(*vertices)[0],vtx_n);
	int e_n=(int)elements.size();
	File::Write_Binary(output,e_n);
	if(e_n>0)File::Write_Binary_Array(output,&elements[0],e_n);
}

template<int d,int e_d> void SimplicialMesh<d,e_d>::Read_Binary(std::istream&input)
{
	int vtx_n=0;
	File::Read_Binary(input,vtx_n);
	(*vertices).resize(vtx_n);
	File::Read_Binary_Array(input,&(*vertices)[0],vtx_n);
	int e_n=0;
	File::Read_Binary(input,e_n);
	if(e_n>0){
		elements.resize(e_n);
		File::Read_Binary_Array(input,&elements[0],e_n);}
}

template<int d,int e_d> void SimplicialMesh<d,e_d>::Write_To_File_3d(const std::string& file_name) const
{
	if constexpr (d==3){
		File::Write_Binary_To_File(file_name,*this);}
	else{
		SimplicialMesh<3,e_d> s3;Dim_Conversion(*this,s3);
		File::Write_Binary_To_File(file_name,s3);}
}

template<int d,int e_d> void SimplicialMesh<d,e_d>::Read_Text(std::istream& input)
{
	int vtx_n=0;File::Read_Text(input,vtx_n);
	if(vtx_n>0){
		(*vertices).resize(vtx_n);
		for(int i=0;i<vtx_n;i++)File::Read_Text_Array(input,(*vertices)[i],d);}
	int e_n=0;File::Read_Text(input,e_n);
	if(e_n>0){
		elements.resize(e_n);
		for(int i=0;i<e_n;i++)File::Read_Text_Array(input,elements[i],e_d);}
	int e_t=0;File::Read_Text(input,e_t);
	if(e_t>0){
		if(Textures()==nullptr)textures.reset(new Array<Vector2>());
		(*textures).resize(e_t);
		for(int i=0;i<e_t;i++)File::Read_Text_Array(input,(*textures)[i],2);}
}

template<int d,int e_d> void SimplicialMesh<d,e_d>::Write_Text(std::ostream& output) const
{
	int vtx_n=(int)(*vertices).size();File::Write_Text(output,vtx_n);File::Write_Text(output,'\n');
	if(vtx_n>0){for(int i=0;i<vtx_n;i++){File::Write_Text_Array(output,(*vertices)[i],d,' ');File::Write_Text(output,'\n');}}
	int e_n=(int)elements.size();
	File::Write_Text(output,'\n');File::Write_Text(output,e_n);File::Write_Text(output,'\n');
	if(e_n>0){for(int i=0;i<e_n;i++){File::Write_Text_Array(output,elements[i],e_d,' ');File::Write_Text(output,'\n');}}
}

template class SimplicialMesh<2,2>;
template class SimplicialMesh<2,3>;
template class SimplicialMesh<2,4>;
template class SimplicialMesh<3,2>;
template class SimplicialMesh<3,3>;
template class SimplicialMesh<3,4>;

////TetrahedronMesh
template<int d> TetrahedronMesh<d>::TetrahedronMesh(const ArrayPtr<VectorD> _vertices):Base(_vertices){}

template class TetrahedronMesh<2>;
template class TetrahedronMesh<3>;

////TriangleMesh
template<int d> TriangleMesh<d>::TriangleMesh(const ArrayPtr<VectorD> _vertices):Base(_vertices){}

template class TriangleMesh<2>;
template class TriangleMesh<3>;

////SegmentMesh
template<int d> SegmentMesh<d>::SegmentMesh(const ArrayPtr<VectorD> _vertices):Base(_vertices){}

template class SegmentMesh<2>;
template class SegmentMesh<3>;

////QuadMesh
template<int d> QuadMesh<d>::QuadMesh(const ArrayPtr<VectorD> _vertices):Base(_vertices){}

template class QuadMesh<2>;
template class QuadMesh<3>;

////HexMesh
template<int d> HexMesh<d>::HexMesh(const ArrayPtr<VectorD> _vertices):Base(_vertices){}

template class HexMesh<2>;
template class HexMesh<3>;

////SurfaceQuadMesh
template<int d> SurfaceQuadMesh<d>::SurfaceQuadMesh(const ArrayPtr<VectorD> _vertices):Base(_vertices){}

template class SurfaceQuadMesh<2>;
template class SurfaceQuadMesh<3>;

template<class MESH_T1,class MESH_T2> void Dim_Conversion(const MESH_T1& mesh2,/*rst*/MESH_T2& mesh3)
{
	////TOACC: with openmp
	mesh3.vertices->resize((int)(mesh2.vertices->size()));
	for(auto i=0;i<(*mesh3.vertices).size();i++)
		AuxFunc::Dim_Conversion<real,MESH_T1::Dim(),MESH_T2::Dim()>((*mesh2.vertices)[i],(*mesh3.vertices)[i],(real)0);
	mesh3.elements.resize((int)mesh2.elements.size());
	for(auto i=0;i<mesh2.elements.size();i++)
		AuxFunc::Dim_Conversion<int,MESH_T1::Element_Dim(),MESH_T2::Element_Dim()>(mesh2.elements[i],mesh3.elements[i],(int)-1);
}

#define Inst_Helper(T1,T2) \
template void Dim_Conversion<T1,T2 >(const T1&,T2&);
Inst_Helper(SegmentMesh<2>,SegmentMesh<3>);
Inst_Helper(QuadMesh<2>,QuadMesh<3>);
#undef Inst_Helper

#define Inst_Helper(T1) \
template<> void Dim_Conversion<T1,T1 >(const T1& mesh_input,T1& mesh3){mesh3=mesh_input;}
Inst_Helper(TetrahedronMesh<3>);
Inst_Helper(TriangleMesh<3>);
Inst_Helper(SegmentMesh<3>);
Inst_Helper(QuadMesh<3>);
Inst_Helper(SurfaceQuadMesh<3>);
#undef Inst_Helper
//////////////////////////////////////////////////////////////////////////
// Mesher
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Mesher_h__
#define __Mesher_h__
#include "Common.h"
#include "AuxFunc.h"
#include "Particles.h"
#include "Hashtable.h"

namespace MesherFunc{
//////////////////////////////////////////////////////////////////////////
////Other node(s) in the element
inline Vector<int,1> Other_Nodes_In_Element(const Vector<int,2>& e,const int p)
{return e[0]==p?Vector<int,1>(e[1]):Vector<int,1>(e[0]);}
inline Vector<int,2> Other_Nodes_In_Element(const Vector<int,3>& e,const int p)
{return e[0]==p?Vector<int,2>(e[1],e[2]):(e[1]==p?Vector<int,2>(e[2],e[0]):Vector<int,2>(e[0],e[1]));}
////return (v0,v1,v2) s.t. (p,v0,v1,v2) is a positive tet or (v0,v1,v2) is a face on the tet pointing outward
inline Vector<int,3> Other_Nodes_In_Element(const Vector<int,4>& e,const int p)	
{return e[0]==p?Vector<int,3>(e[1],e[2],e[3]):(e[1]==p?Vector<int,3>(e[2],e[0],e[3]):(e[2]==p?Vector<int,3>(e[3],e[0],e[1]):Vector<int,3>(e[0],e[2],e[1])));}
inline int Other_Node_In_Element(const Vector<int,2>&e,const int p1)
{return e[0]==p1?e[1]:e[0];}
inline int Other_Node_In_Element(const Vector<int,3>& e,const int p1,const int p2)
{return (e[0]!=p1&&e[0]!=p2)?e[0]:(e[1]!=p1&&e[1]!=p2)?e[1]:e[2];}
inline int Other_Node_In_Element(const Vector<int,4>& e,const int p1,const int p2,const int p3)
{for(int i=0;i<4;i++)if(e[i]!=p1&&e[i]!=p2&&e[i]!=p3)return e[i];return 0;}
////return (b,d), s.t. (a,b,c,d) is a positive tet
inline Vector<int,2> Other_Nodes_In_Element(const Vector<int,4>& e,const int a,const int c)	
{Vector4i tet=Reordered(e,a);Vector3i tri=Reordered(Vector3i(tet[1],tet[2],tet[3]),c);int d=tri[1];int b=tri[2];return Vector2i(b,d);}

////Replace vertex in element
template<int dim> inline Vector<int,dim> Replaced(const Vector<int,dim>& e,const int v0,const int v1) ////replace v0 with v1
{Vector<int,dim> replaced_e=e;for(int i=0;i<dim;i++)if(replaced_e[i]==v0)replaced_e[i]=v1;return replaced_e;}
template<int dim> inline void Replace_Vertex(Vector<int,dim>& e,const int v0,const int v1) ////replace v0 with v1
{for(int i=0;i<dim;i++)if(e[i]==v0)e[i]=v1;}
template<int dim> inline Vector<int,dim> Replaced_Two(const Vector<int,dim>& e,const int v0,const int v1,const int v2) ////replace v0 and v1 with v2
{Vector<int,dim> replaced_e=e;for(int i=0;i<dim;i++)if(replaced_e(i)==v0||replaced_e(i)==v1)replaced_e(i)=v2;return replaced_e;} 
template<int dim> inline bool Replaced_Two(const Vector<int,dim>& e,Vector<int,dim>& replaced_e,const int v0,const int v1,const int v2) ////replace v0 and v1 with v2
{replaced_e=e;bool replaced=false;for(int i=0;i<dim;i++)if(replaced_e[i]==v0||replaced_e[i]==v1){replaced_e[i]=v2;replaced=true;}return replaced;} 

////Convert element dimension
template<int d2> inline Vector<int,d2> To_V(const Vector4i& input){Vector<int,d2> output=Vector<int,d2>::Zero();for(int i=0;i<d2;i++)output[i]=input[i];return output;}
template<> inline Vector2i To_V<2>(const Vector4i& input){return Vector2i(input[0],input[1]);}
template<> inline Vector3i To_V<3>(const Vector4i& input){return Vector3i(input[0],input[1],input[2]);}
inline Vector2i To_V2(const Vector4i& input){return Vector2i(input[0],input[1]);}
inline Vector3i To_V3(const Vector4i& input){return Vector3i(input[0],input[1],input[2]);}
inline Vector4i To_V4(const Vector2i& input){return Vector4i(input[0],input[1],-1,-1);}
inline Vector4i To_V4(const Vector3i& input){return Vector4i(input[0],input[1],input[2],-1);}
inline Vector4i To_V4(int i,int j){return Vector4i(i,j,-1,-1);}
inline Vector4i To_V4(int i,int j,int k){return Vector4i(i,j,k,-1);}

////Check whether vertex exists
inline bool Element_Has(const Vector<int,3>& e,const Vector<int,2>& f){return AuxFunc::Has(e,f[0])&&AuxFunc::Has(e,f[1]);}
inline bool Element_Has(const Vector<int,4>& e,const Vector<int,3>& f){return AuxFunc::Has(e,f[0])&&AuxFunc::Has(e,f[1])&&AuxFunc::Has(e,f[2]);}
inline bool Element_Has_Not(const Vector<int,2>& e,const int f){return e[0]!=f||e[1]!=f;}
inline bool Element_Has_Not(const Vector<int,3>& e,const int f){return e[0]!=f||e[1]!=f||e[2]!=f;}
inline bool Element_Has_Not(const Vector<int,4>& e,const int f){return e[0]!=f||e[1]!=f||e[2]!=f||e[3]!=f;}

////Check shared vertices
inline int Shared_Vertex(const Vector<int,3>& e1,const Vector<int,3>& e2,const Vector<int,3>& e3)
{return (AuxFunc::Has(e2,e1[0])&&AuxFunc::Has(e3,e1[0]))?e1[0]:((AuxFunc::Has(e2,e1[1])&&AuxFunc::Has(e3,e1[1]))?e1[1]:(AuxFunc::Has(e2,e1[2])&&AuxFunc::Has(e2,e1[2])?e1[2]:-1));}
inline int Shared_Vertex(const Vector<int,2>& e1,const Vector<int,2>& e2)
{return AuxFunc::Has(e2,e1[0])?e1[0]:(AuxFunc::Has(e2,e1[1])?e1[1]:-1);}
inline Vector<int,2> Shared_Edge(const Vector<int,3>& e1,const Vector<int,3>& e2)
{if(!AuxFunc::Has(e2,e1[0]))return Vector<int,2>(e1[1],e1[2]);else if(!AuxFunc::Has(e2,e1[1]))return Vector<int,2>(e1[2],e1[0]);else return Vector<int,2>(e1[0],e1[1]);}
};

template<int d,class T_MESH> class MeshAux
{Typedef_VectorDi(d);
public:

	const T_MESH& mesh;
	HashtableMultiValue<int,int> vertex_element_hashtable;
	HashtableMultiValue<Vector2i,int> edge_element_hashtable;

	MeshAux(const T_MESH& _mesh):mesh(_mesh){}

	virtual void Initialize(){}
	virtual void Update_Vertex_Element_Hashtable(){}
	virtual void Update_Edge_Element_Hashtable(){}
};

template<int d,class T_MESH,class T_PARTICLES> class Mesher
{Typedef_VectorD(d);
public:
	T_MESH& mesh;
	T_PARTICLES& particles;

	Mesher(T_MESH& _mesh,T_PARTICLES& _particles):mesh(_mesh),particles(_particles){}

	virtual void Initialize(){}
	virtual void Initialize_Parameters(){}
	virtual void Update(){}
	virtual void Dbg_Output() const {}

protected:
	Heap<int,std::greater<int> > invalid_particle_heap;
	Hashset<int> invalid_particle_hashset;
	Heap<int,std::greater<int> > invalid_element_heap;
	Hashset<int> invalid_element_hashset;

	////Basic operations
	real Length(const Vector2i& edge) const {return (particles.X(edge[0])-particles.X(edge[1])).norm();}

	//////////////////////////////////////////////////////////////////////////
	////Element operations
	int Set_Element(const int e_idx,const Vector4i& vtx){mesh.elements[e_idx]=vtx;return e_idx;}
	void Remove_Element(const int e_idx){invalid_element_heap.push(e_idx);invalid_element_hashset.insert(e_idx);}
	int Append_Element(const Vector4i& vtx){mesh.elements.push_back(vtx);return (int)mesh.elements.size()-1;}
public: 
	bool Valid_Element(const int e_idx) const {return invalid_element_hashset.find(e_idx)==invalid_element_hashset.end();}

	int Add_Element(const Vector4i& vtx)
	{
		if(!invalid_element_heap.empty()){
			int e_idx=invalid_element_heap.top();invalid_element_heap.pop();invalid_element_hashset.erase(e_idx);
			Set_Element(e_idx,vtx);return e_idx;}
		else return Append_Element(vtx);
	}

	//////////////////////////////////////////////////////////////////////////
	////Particle operations
public:
	bool Valid_Particle(const int i) const {return invalid_particle_hashset.find(i)==invalid_particle_hashset.end();}
	virtual int Add_Particle()
	{if(!invalid_particle_heap.empty()){int p=invalid_particle_heap.top();invalid_particle_heap.pop();invalid_particle_hashset.erase(p);return p;}
	else return Append_Particle();}
protected:
	void Remove_Particle(const int i){invalid_particle_heap.push(i);invalid_particle_hashset.insert(i);}
	int Append_Particle(){particles.Add_Element();return particles.Size()-1;}

	int Split_Particles(const int a,const int c)
	{
		int v=Add_Particle();
		particles.X(v)=(real).5*(particles.X(a)+particles.X(c));
		particles.V(v)=(real).5*(particles.V(a)+particles.V(c));
		particles.M(v)=particles.M(a)+particles.M(c);
		return v;
	}

	int Collapse_Particles(const int a,const int c)
	{
		Remove_Particle(a);Remove_Particle(c);int v=Add_Particle();
		particles.X(v)=(real).5*(particles.X(a)+particles.X(c));
		particles.V(v)=(real).5*(particles.V(a)+particles.V(c));
		particles.M(v)=particles.M(a)+particles.M(c);
		return v;
	}

	//////////////////////////////////////////////////////////////////////////
	////Get vertices from element
	void Element_Vertices(const int e_idx,int& v0,int& v1) const {v0=mesh.elements[e_idx][0];v1=mesh.elements[e_idx][1];}
	void Element_Vertices(const int e_idx,int& v0,int& v1,int& v2) const {v0=mesh.elements[e_idx][0];v1=mesh.elements[e_idx][1];v2=mesh.elements[e_idx][2];}
	Vector2i Segment_Vertices(const int e_idx) const {return MesherFunc::To_V2(mesh.elements[e_idx]);}
	Vector3i Triangle_Vertices(const int e_idx) const {return MesherFunc::To_V3(mesh.elements[e_idx]);}
	template<class T_ARRAY> int Element_Vertices_Edges(const Vector4i& v,const int ds,const int v0,T_ARRAY& edges) const
	{int c=0;for(int i=0;i<ds;i++){if(v[i]!=v0)edges[c++]=Vector2i(v0,v[i]);}return c;}
	template<int dim> int Element_Vertex_Index(const Vector<int,dim>& v,const int v0) const
	{for(int i=0;i<dim;i++)if(v[i]==v0)return i;return -1;}
};
#endif
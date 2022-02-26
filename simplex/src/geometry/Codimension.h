//////////////////////////////////////////////////////////////////////////
// Codimension
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __CodimMesh_h__
#define __CodimMesh_h__
#include <algorithm>
#include "Mesh.h"
#include "File.h"
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "Mesher.h"
#include "Constants.h"
#include "SparseFunc.h"
#include "Particles.h"

namespace Codimension{

template<int d> using CodimParticles=Particles<d>;

template<int d> class CodimMesh : public SimplicialMesh<d,4>
{Typedef_VectorDii(d);Typedef_VectorEi(4);typedef SimplicialMesh<d,4> Base;
public:
	using Base::vertices;using Base::elements;
	using Base::Vertices;using Base::Elements;

	CodimMesh(const ArrayPtr<VectorD> _vertices=nullptr):Base(_vertices){}

	void Add_Element(const int v0,const int v1){elements.push_back(v0,v1,-1,-1);}
	void Add_Element(const int v0,const int v1,const int v2){elements.push_back(v0,v1,v2,-1);}
	void Add_Element(const int v0,const int v1,const int v2,const int v3){elements.push_back(v0,v1,v2,v3);}
	void Add_Element(const Vector2i& vtx){elements.push_back(Vector4i(vtx[0],vtx[1],-1,-1));}
	void Add_Element(const Vector3i& vtx){elements.push_back(Vector4i(vtx[0],vtx[1],vtx[2],-1));}
	void Add_Element(const Vector4i& vtx){elements.push_back(vtx);}
	int Dim(const int i) const {return Dim(elements[i]);}
	static int Dim(const Vector4i& e){return e[3]!=-1?4:(e[2]!=-1?3:2);}

	void Element_Vertices(const int idx,ArrayF<VectorD,d+1>& vtx) const {for(int i=0;i<d+1;i++)vtx[i]=Vertices()[Elements()[idx][i]];}	////assuming the element is with d+1 dimension
	void Element_Vertices(const int idx,Array<VectorD>& vtx) const {int ds=Dim(idx);vtx.resize(ds);for(int i=0;i<ds;i++)vtx[i]=Vertices()[Elements()[idx][i]];}
};

////Auxiliary structures
template<int d> class CodimMeshAux : public MeshAux<d,CodimMesh<d> >
{Typedef_VectorDi(d);typedef MeshAux<d,CodimMesh<d> > Base;
public:
	using Base::vertex_element_hashtable;using Base::edge_element_hashtable;
	using Base::mesh;
	
	HashtableMultiValue<Vector3i,int> face_element_hashtable;
	HashtableMultiValue<int,int> vertex_ancestor_hashtable;

	CodimMeshAux(const CodimMesh<d>& _mesh):Base(_mesh){}
	virtual void Initialize();

	////hashtable operations
	virtual void Update_Vertex_Element_Hashtable();
	virtual void Update_Edge_Element_Hashtable();
	void Update_Face_Element_Hashtable();
	void Add_New_Element_To_Vertex_Element_Hashtable(const Vector4i& vtx,const int i,const int ds);
	void Add_New_Element_To_Edge_Element_Hashtable(const Vector4i& vtx,const int i,const int ds);
	void Add_New_Element_To_Face_Element_Hashtable(const Vector4i& vtx,const int i,const int ds);
	void Add_New_Element_To_Hashtables(const Vector4i& vtx,const int i,const int ds);
	void Vertex_Incident_Vertices(const int v,Hashset<int>& incident_vertex_hashset) const;
	void Remove_Element_From_Hashtables(const int i,const int ds);

	////boundary
	int Boundary_Type(const int v) const;
	bool Is_On_Boundary(const int v) const;
	bool On_Dimension_Boundary_2(const int v) const;
	bool On_Dimension_Boundary_3(const int v) const;
	bool On_Dimension_Boundary_4(const int v) const;
	bool On_Dimension_Boundary_3(const Vector2i& edge) const;
};

////Mesher
template<int d> class CodimMesher : public Mesher<d,CodimMesh<d>,CodimParticles<d> >
{Typedef_VectorD(d);typedef Mesher<d,CodimMesh<d>,CodimParticles<d> > Base;
public:
	using Base::particles;using Base::mesh;using Base::invalid_particle_hashset;using Base::invalid_particle_heap;

	using Base::Valid_Particle;using Base::Valid_Element;using Base::Element_Vertices_Edges;using Base::Split_Particles;using Base::Element_Vertices;
	using Base::Collapse_Particles;using Base::Length;using Base::Triangle_Vertices;using Base::Remove_Particle;using Base::Append_Particle;
	using Base::Set_Element;using Base::Remove_Element;using Base::Add_Element;

	CodimMeshAux<d> aux;

	////meshing parameters
	real length_scale=(real)1;
	real edge_split_ratio=(real)1.5;
	real edge_collapse_ratio=(real).5;
	real small_angle=pi*(real).1;
	real edge_snap_ratio=(real).1;
	real edge_flip_dihedral_angle=pi*(real).75;

	////meshing options
	bool transition=true;
	bool verbose=false;
	bool dbg=false;
	bool use_triangle_edge_split=true;
	bool use_triangle_edge_collapse=true;
	bool use_triangle_edge_snap=true;	////not the basic operation
	bool use_boundary_aware_collapse=true;
	bool track_vertex_ancestor=true;

	CodimMesher(CodimMesh<d>& _mesh,CodimParticles<d>& _particles):Base(_mesh,_particles),aux(_mesh){}
	virtual void Initialize();
	virtual void Update();	////local topological repairs

public:
	////Checking operations
	bool Check_Split_Segment_Edge(const Vector2i& edge) const;
	bool Check_Collapse_Segment_Edge(const Vector2i& edge) const;
	bool Check_Split_Triangle_Edge(const Vector2i& edge) const; 
	bool Check_Collapse_Triangle_Edge(const Vector2i& edge,const int e) const; 
	bool Check_Snap_Triangle_Edge(const Vector2i& edge,const int e) const;

	////Atomic meshing operations
	////Atomic segment
	void Split_Segment_Edge(const int i);
	void Collapse_Segment_Edge(const int i);
	////Atomic triangle
	int Split_Triangle_And_Tetrahedron_Edge(const Vector2i& edge);
	bool Collapse_Triangle_Edge(const Vector2i& ac,bool enforce_op=false);
	bool Flip_Triangle_Edge(const Vector2i& ac,bool enforce_op=false);	////assuming ac is ordered
	bool Snap_Triangle_Edge(const Vector2i& ac,const int e);
	////Atomic add/remove
	int Add_Triangle(const ArrayF<VectorD,3>& vtx_pos);
	bool Remove_Triangle(const int e);
	int Add_Segment(const ArrayF<VectorD,2>& vtx_pos);
	bool Remove_Segment(const int e);
	////Atomic aux functions
	bool Has_Inverted_Triangle(const int a,const VectorD& new_pos,const List<int>& a_hat_elements) const;
	bool Flip_Improve_Triangle_Quality(const Vector2i& ac,const int b,const int g,int& e_i,int& e_j) const;
	virtual int Add_Particle();
	////Prune operations
	void Pruned_Copy(CodimMesh<d>& pruned_mesh,CodimParticles<d>& pruned_particles,std::shared_ptr<Array<int> >* remap=nullptr);
	void Prune(std::shared_ptr<Array<int> >* remap=nullptr);
	////Length scale
	void Set_Default_Length_Scale();
	////Hierarchy
	void Update_Vertex_Ancestor_Hashtable();
	void Coarsen(const real coarsen_factor);
	void Build_Hierarchy(const int levels,const real coarsen_factor,Array<CodimMesh<d> >& hier_meshes/*,Array<SparseMatrixT>& mappings*/);

protected:
	////State buffer
	Array<int> buffer_elements=Array<int>(16,0);
	void Clear_Buffer_Elements(){buffer_elements.clear();}
	void Add_Buffer_Element(const int i){buffer_elements.push_back(i);}
	////Data access
	inline const VectorD& X(const int i) const {return particles.X(i);}
	inline VectorD& X(const int i){return particles.X(i);}
	inline const Vector4i& E(const int i) const {return mesh.elements[i];}
	inline Vector4i& E(const int i){return mesh.elements[i];}

protected:
	////Helper functions
	void Vertex_Incident_Edges(const int v,const List<int>& incident_elements,Hashset<Vector2i>& edge_hashset);
	void Vertex_Incident_Edges(const int v,Hashset<Vector2i>& edge_hashset);
	bool Edge_Incident_Triangle_Pair(const Vector2i& ac,int& b,int& dd,int& e_i,int& e_j) const;
	int Edge_Type(const List<int>& incident_elements) const;
	int Edge_Type(const Vector2i& edge) const;
	real Measure_Triangle_Quality(const Array<Vector3i>& tris_v) const;
	real Default_Length_Scale() const;
	
public:
	////Dbg
	void Dbg_Output() const;
	bool Dbg_Check(const std::string& info) const;
	bool Validate() const;
	bool Validate_Triangle_Order() const;
	bool Validate_Vertex_Element_Hashtable() const;
	bool Validate_Edge_Element_Hashtable() const;
	bool Validate_Dihedral_Angles() const;
	bool Validate_Normal_2D() const;	////2D only
	bool Validate_Tetrahedron_Order() const;
	bool Validate_Coarsening_Mapping() const;
};

////Element iterator
template<int d> class CodimMeshElementIterator
{Typedef_VectorDii(d);
	const CodimMesh<d>* mesh=nullptr;
	const CodimMesher<d>* mesher=nullptr;
	const Array<Vector4i>* elements=nullptr;
	int index=-1;
public:
	CodimMeshElementIterator(const CodimMesh<d>& _mesh,const CodimMesher<d>& _mesher)
		:mesh(&_mesh),mesher(&_mesher),elements(&mesh->Elements()){Next();}
	CodimMeshElementIterator(const CodimMesh<d>& _mesh,const CodimMesher<d>& _mesher,const Array<Vector4i>& _elements)
		:mesh(&_mesh),mesher(&_mesher),elements(&_elements){Next();}

	bool Valid() const {return index<elements->size();}
	void Next();
	int Index() const {return index;}

	int Dim() const {return mesh->Dim(index);}
	const Vector4i& Element() const {return mesh->Elements()[index];}
	void Element_And_Ds(int& ds,Vector4i& e) const {ds=mesh->Dim(index);e=mesh->Elements()[index];}
	void Element_Vertices(ArrayF<VectorD,d+1>& vtx) const {mesh->Element_Vertices(index,vtx);}
	void Element_Vertices(Array<VectorD>& vtx) const {mesh->Element_Vertices(index,vtx);}
};

////Macros for using GridIterator
#define iterate_codim_mesh_ele(iter,mesh,mesher) \
for(CodimMeshElementIterator<d> iter(mesh,mesher);iter.Valid();iter.Next())
#define iterate_codim_ele_array(iter,mesh,mesher,ele) \
for(CodimMeshElementIterator<d> iter(mesh,mesher,ele);iter.Valid();iter.Next())

////Helper functions
void Split_Dimension(const CodimMesh<3>& mesh,SegmentMesh<3>& seg_mesh,TriangleMesh<3>& tri_mesh,TetrahedronMesh<3>& tet_mesh);
void Split_Dimension(const CodimMesh<2>& mesh,SegmentMesh<2>& seg_mesh,TriangleMesh<2>& tri_mesh);

void Merge_Dimension(CodimMesh<3>& mesh,const SegmentMesh<3>* seg_mesh=nullptr,const TriangleMesh<3>* tri_mesh=nullptr,const TetrahedronMesh<3>* tet_mesh=nullptr);
void Merge_Dimension(CodimMesh<2>& mesh,const SegmentMesh<2>* seg_mesh=nullptr,const TriangleMesh<2>* tri_mesh=nullptr,const TetrahedronMesh<2>* tet_mesh=nullptr);
}
#endif
//////////////////////////////////////////////////////////////////////////
// Codimension
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include <stack>
#include <queue>
#include "Codimension.h"

namespace Codimension{
using namespace AuxFunc;
using namespace MeshFunc;
using namespace MesherFunc;

inline bool Has_Duplicate_Vertices(const Vector<int,3>& e){return e[0]==e[1]||e[0]==e[2]||e[1]==e[2];}

inline Vector<int,3> Simplex_Opposite_Face(const Vector<int,4>& e,const int v)
{
    if(v==e[0])return Vector<int,3>(e[1],e[2],e[3]);
    else if(v==e[1])return Vector<int,3>(e[0],e[3],e[2]);
    else if(v==e[2])return Vector<int,3>(e[0],e[1],e[3]);
    else return Vector<int,3>(e[0],e[2],e[1]);
}

////split tet abcd to vbcd and vabd
inline void Split_Tetrahedron_Vertices(const Vector<int,4>& tet_vtx,const int a,const int c,const int v,Vector<int,4>& vbcd,Vector<int,4>& vabd)
{   
    Vector<int,2> b_and_d=Other_Nodes_In_Element(tet_vtx,a,c);int b=b_and_d[0];int d=b_and_d[1];
    if(Unique_Ordered(tet_vtx)!=Unique_Ordered(Vector<int,4>(a,b,c,d))){int tmp=b;b=d;d=tmp;}
    vbcd=Vector<int,4>(v,b,c,d);vabd=Vector<int,4>(v,a,b,d);
}

////split tri acd to avd and vcd (or adv and vdc)
inline void Split_Triangle_Vertices(const Vector<int,3>& tri_vtx,const int a,const int c,const int v,Vector<int,3>& avd,Vector<int,3>& vcd)
{
    int d=Other_Node_In_Element(tri_vtx,a,c);bool flip=(Unique_Ordered(tri_vtx)!=Unique_Ordered(Vector<int,3>(a,c,d)));
    if(!flip){avd=Vector<int,3>(a,v,d);vcd=Vector<int,3>(v,c,d);}else{avd=Vector<int,3>(a,d,v);vcd=Vector<int,3>(v,d,c);}
}

////split tet abcd 
inline void Split_Tetrahedron_Vertices(const Vector<int,4>& tet,const Vector<int,3>& face,const int v,Vector<int,4>& adbv,Vector<int,4>& bdcv,Vector<int,4>& cdav)
{
    int a=face[0];int b=face[1];int c=face[2];int d=Other_Node_In_Element(tet,a,b,c);if(Unique_Ordered(tet)!=Unique_Ordered(Vector<int,4>(a,b,c,d))){int tmp=b;b=c;c=tmp;}
    adbv=Vector<int,4>(a,d,b,v);bdcv=Vector<int,4>(b,d,c,v);cdav=Vector<int,4>(c,d,a,v);
}

real Volume_Length_Measure(const Vector<real,3>& a,const Vector<real,3>& b,const Vector<real,3>& c,const Vector<real,3>& d)
{
    real edge_rms=sqrt(((a-b).squaredNorm()+(a-c).squaredNorm()+(a-d).squaredNorm()+(b-c).squaredNorm()+(b-d).squaredNorm()+(c-d).squaredNorm())/(real)6);
    return Tetrahedron_Volume(a,b,c,d)/pow(edge_rms,3)*(real)6*sqrt((real)2);
}

template<int d> real Area_Length_Measure(const Vector<real,d>& a,const Vector<real,d>& b,const Vector<real,d>& c)
{
    real edge_rms=sqrt(((a-b).squaredNorm()+(a-c).squaredNorm()+(b-c).squaredNorm())*(real)one_third);
    return Triangle_Area(a,b,c)/pow(edge_rms,2)*(real)4/(real)3*sqrt((real)3);
}
template real Area_Length_Measure<2>(const Vector2&,const Vector2&,const Vector2&);
template real Area_Length_Measure<3>(const Vector3&,const Vector3&,const Vector3&);

inline real Min_Cosine_Angle(const Vector<real,3>& a,const Vector<real,3>& b,const Vector<real,3>& c)
{
    Vector<real,3> ab=(b-a).normalized();Vector<real,3> bc=(c-b).normalized();Vector<real,3> ca=(a-c).normalized();
    return (std::min)(-ab.dot(ca),(std::min)(-ab.dot(bc),-ca.dot(bc)));
}

template<class real> inline real Max_Cosine_Angle(const Vector<real,3>& a,const Vector<real,3>& b,const Vector<real,3>& c)
{
    Vector<real,3> ab=(b-a).normalized();Vector<real,3> bc=(c-b).normalized();Vector<real,3> ca=(a-c).normalized();
    return (std::max)(-ab.dot(ca),(std::max)(-ab.dot(bc),-ca.dot(bc)));
}

////types
enum class BoundaryTypeName : int {
	ParticleOnSegmentMesh=0x1,
	ParticleOnSegmentMeshBoundary=0x2,
	ParticleOnTriangleMesh=0x4,
	ParticleOnTriangleMeshBoundary=0x8,
	ParticleOnTetrahedronMesh=0x10,
	ParticleOnTetrahedronMeshBoundary=0x20};

inline void Set_On_Segment_Mesh(int& type){type|=(int)BoundaryTypeName::ParticleOnSegmentMesh;}
inline void Set_On_Triangle_Mesh(int& type){type|=(int)BoundaryTypeName::ParticleOnTriangleMesh;}
inline void Set_On_Tetrahedron_Mesh(int& type){type|=(int)BoundaryTypeName::ParticleOnTetrahedronMesh;}
inline void Set_On_Segment_Mesh_Boundary(int& type){type|=(int)BoundaryTypeName::ParticleOnSegmentMeshBoundary;}
inline void Set_On_Triangle_Mesh_Boundary(int& type){type|=(int)BoundaryTypeName::ParticleOnTriangleMeshBoundary;}
inline void Set_On_Tetrahedron_Mesh_Boundary(int& type){type|=(int)BoundaryTypeName::ParticleOnTetrahedronMeshBoundary;}
inline bool Is_On_Segment_Mesh(const int type){return (type & (int)BoundaryTypeName::ParticleOnSegmentMesh)!=0;}
inline bool Is_On_Triangle_Mesh(const int type){return (type & (int)BoundaryTypeName::ParticleOnTriangleMesh)!=0;}
inline bool Is_On_Tetrahedron_Mesh(const int type){return (type & (int)BoundaryTypeName::ParticleOnTetrahedronMesh)!=0;}
inline bool Is_On_Segment_Mesh_Boundary(const int type){return (type & (int)BoundaryTypeName::ParticleOnSegmentMeshBoundary)!=0;}
inline bool Is_On_Triangle_Mesh_Boundary(const int type){return (type & (int)BoundaryTypeName::ParticleOnTriangleMeshBoundary)!=0;}
inline bool Is_On_Tetrahedron_Mesh_Boundary(const int type){return (type & (int)BoundaryTypeName::ParticleOnTetrahedronMeshBoundary)!=0;}
inline bool Is_On_Segment_Mesh_Interior(const int type){return Is_On_Segment_Mesh(type) && !Is_On_Segment_Mesh_Boundary(type);}
inline bool Is_On_Triangle_Mesh_Interior(const int type){return Is_On_Triangle_Mesh(type) && !Is_On_Triangle_Mesh_Boundary(type);}
inline bool Is_On_Tetrahedron_Mesh_Interior(const int type){return Is_On_Tetrahedron_Mesh(type) && !Is_On_Tetrahedron_Mesh_Boundary(type);}
inline bool Is_Isolated(const int type){return type==0;}
inline bool Is_On_Codimensional_Mesh_Free_Surface(const int type){return Is_On_Segment_Mesh(type)||Is_On_Triangle_Mesh(type)||Is_On_Tetrahedron_Mesh_Boundary(type);}
inline bool Is_On_Codimensional_Mesh_Boundary(const int type){
    return (Is_On_Segment_Mesh_Boundary(type)&&!Is_On_Triangle_Mesh_Interior(type)&&!Is_On_Tetrahedron_Mesh_Interior(type))||
    (Is_On_Triangle_Mesh_Boundary(type)&&!Is_On_Tetrahedron_Mesh_Interior(type)&&!Is_On_Segment_Mesh_Interior(type))||
    (Is_On_Segment_Mesh_Boundary(type)&&!Is_On_Triangle_Mesh_Interior(type)&&!Is_On_Tetrahedron_Mesh_Interior(type));}

//////////////////////////////////////////////////////////////////////////
////CodimMeshAux
template<int d> void CodimMeshAux<d>::Initialize()
{
	Update_Vertex_Element_Hashtable();
	Update_Edge_Element_Hashtable();
}

template<int d> void CodimMeshAux<d>::Update_Vertex_Element_Hashtable()
{
	vertex_element_hashtable.clear();
	for(int i=0;i<(int)mesh.elements.size();i++){const Vector4i& vtx=mesh.elements[i];const int ds=mesh.Dim(i);
		Add_New_Element_To_Vertex_Element_Hashtable(vtx,i,ds);}
}

template<int d> void CodimMeshAux<d>::Update_Edge_Element_Hashtable()
{
	edge_element_hashtable.clear();
	for(int i=0;i<(int)mesh.elements.size();i++){const Vector4i& vtx=mesh.elements[i];const int ds=mesh.Dim(i);
		Add_New_Element_To_Edge_Element_Hashtable(vtx,i,ds);}
}

template<int d> void CodimMeshAux<d>::Update_Face_Element_Hashtable()
{
	face_element_hashtable.clear();
	for(int i=0;i<(int)mesh.elements.size();i++){const Vector4i& vtx=mesh.elements[i];
		const int ds=mesh.Dim(i);Add_New_Element_To_Face_Element_Hashtable(vtx,i,ds);}
}

template<int d> void CodimMeshAux<d>::Add_New_Element_To_Vertex_Element_Hashtable(const Vector4i& vtx,const int i,const int ds)
{for(int j=0;j<ds;j++){Add(vertex_element_hashtable,vtx[j],i);}}

template<int d> void CodimMeshAux<d>::Add_New_Element_To_Edge_Element_Hashtable(const Vector4i& vtx,const int i,const int ds)
{ArrayF<Vector2i,6> edges;int edge_n=Element_Edges(vtx,ds,edges);
for(int j=0;j<edge_n;j++){Add(edge_element_hashtable,Unique_Ordered(edges[j]),i);}}

template<int d> void CodimMeshAux<d>::Add_New_Element_To_Face_Element_Hashtable(const Vector4i& vtx,const int i,const int ds)
{ArrayF<Vector3i,4> faces;int face_n=Element_Faces(vtx,ds,faces);
for(int j=0;j<face_n;j++){Add(face_element_hashtable,Sorted(faces[j]),i);}}

template<int d> void CodimMeshAux<d>::Add_New_Element_To_Hashtables(const Vector4i& vtx,const int i,const int ds)
{
	Add_New_Element_To_Vertex_Element_Hashtable(vtx,i,ds);
	Add_New_Element_To_Edge_Element_Hashtable(vtx,i,ds);
}

template<int d> void CodimMeshAux<d>::Vertex_Incident_Vertices(const int v,Hashset<int>& incident_vertex_hashset) const
{
	List<int> vertex_elements;Value_List(vertex_element_hashtable,v,vertex_elements);
	for(const int& i:vertex_elements){const Vector4i& vtx=mesh.elements[i];const int ds=mesh.Dim(i);
		for(int j=0;j<ds;j++){if(vtx[j]!=v)incident_vertex_hashset.insert(vtx[j]);}}
}

template<int d> void CodimMeshAux<d>::Remove_Element_From_Hashtables(const int i,const int ds)
{
	const Vector4i& vtx=mesh.elements[i];
	for(int j=0;j<ds;j++){Remove(vertex_element_hashtable,vtx[j],i);}
	ArrayF<Vector2i,6> edges;int edge_n=Element_Edges(vtx,ds,edges);
	for(int j=0;j<edge_n;j++)Remove(edge_element_hashtable,Unique_Ordered(edges[j]),i);
}

template<int d> int CodimMeshAux<d>::Boundary_Type(const int v) const
{
	int type=0;
	List<int> vertex_elements;Value_List(vertex_element_hashtable,v,vertex_elements);
	for(const int& i:vertex_elements){const int ds=mesh.Dim(i);
		switch(ds){
		case 2:Set_On_Segment_Mesh(type);break;
		case 3:Set_On_Triangle_Mesh(type);break;
		case 4:Set_On_Tetrahedron_Mesh(type);break;}}
	if(Is_On_Segment_Mesh(type)){bool is_on_seg_boundary=On_Dimension_Boundary_2(v);if(is_on_seg_boundary)Set_On_Segment_Mesh_Boundary(type);}
	if(Is_On_Triangle_Mesh(type)){bool is_on_tri_boundary=On_Dimension_Boundary_3(v);if(is_on_tri_boundary)Set_On_Triangle_Mesh_Boundary(type);}
	if(Is_On_Tetrahedron_Mesh(type)){bool is_on_tet_boundary=On_Dimension_Boundary_4(v);if(is_on_tet_boundary)Set_On_Tetrahedron_Mesh_Boundary(type);}
	return type;
}

template<int d> bool CodimMeshAux<d>::Is_On_Boundary(const int v) const
{
	int b=Boundary_Type(v);
	return Is_On_Segment_Mesh_Boundary(b)||Is_On_Triangle_Mesh_Boundary(b)||Is_On_Tetrahedron_Mesh_Boundary(b);
}

template<int d> bool CodimMeshAux<d>::On_Dimension_Boundary_2(const int v) const
{
	List<int> vertex_elements;Value_List(vertex_element_hashtable,v,vertex_elements);
	int seg_n=0;for(const int& i:vertex_elements){const int ds=mesh.Dim(i);if(ds==2)seg_n++;}return seg_n==1;
}

template<int d> bool CodimMeshAux<d>::On_Dimension_Boundary_3(const int v) const
{
	Hashset<int> v_incident_vertices;Vertex_Incident_Vertices(v,v_incident_vertices);
	for(const int& i:v_incident_vertices){if(On_Dimension_Boundary_3(Unique_Ordered(v,i)))return true;}return false;
}

template<int d> bool CodimMeshAux<d>::On_Dimension_Boundary_4(const int v) const
{
	std::cerr<<"On_Dimension_Boundary_4 TOIMPL"<<std::endl;return false;
}

template<int d> bool CodimMeshAux<d>::On_Dimension_Boundary_3(const Vector2i& edge) const
{
	List<int> edge_elements;Value_List(edge_element_hashtable,edge,edge_elements);
	int tri_n=0;for(const int& i:edge_elements){const int ds=mesh.Dim(i);if(ds==3)tri_n++;}return tri_n==1;		
}

template class CodimMeshAux<2>;
template class CodimMeshAux<3>;

template<int d> void CodimMesher<d>::Initialize()
{
	aux.Initialize();

	////initialize ancestor hashtable
	if(track_vertex_ancestor)Update_Vertex_Ancestor_Hashtable();
}

template<int d> void CodimMesher<d>::Update_Vertex_Ancestor_Hashtable()
{
	aux.vertex_ancestor_hashtable.clear();
		for(int i=0;i<particles.Size();i++){if(!Valid_Particle(i))continue;
			Add(aux.vertex_ancestor_hashtable,i,i);}	
}

template<int d> void CodimMesher<d>::Set_Default_Length_Scale()
{
	length_scale=Default_Length_Scale();
}

////Calculate the averaged edge length
template<int d> real CodimMesher<d>::Default_Length_Scale() const
{
	real length=(real)0;int n=0;
	for(auto i=0;i<mesh.elements.size();i++){
		if(!Valid_Element(i))continue;const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		switch(ds){
			case 2:{////segment operations
				const Vector2i v=To_V2(vtx);
				length+=(X(v[0])-X(v[1])).norm();n++;
			}break;
			case 3:{
				const Vector3i v=To_V3(vtx);
				length+=((X(v[0])-X(v[1])).norm()+(X(v[1])-X(v[2])).norm()+(X(v[2])-X(v[0])).norm());n+=3;			
			}break;
			case 4:break;}}
	if(n>0)return length/(real)n;return length;
}

template<int d> void CodimMesher<d>::Update()
{
	int op_num=0;int pre_op_num=0;int max_op_num=std::numeric_limits<int>::max();int min_op_num=3;int max_iter_num=10;

	for(int iter=0;iter<max_iter_num;iter++){op_num=0;
		for(int i=0;i<mesh.elements.size();i++){if(!Valid_Element(i))continue;
			const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
			ArrayF<Vector2i,6> edges;int edge_n=Element_Edges(vtx,ds,edges);

			switch(ds){
			case 2:{	////segment operations
				const Vector2i edge=Unique_Ordered(To_V2(vtx));
				if(Check_Split_Segment_Edge(edge)){Split_Segment_Edge(i);}
				//else if(Check_Collapse_Segment_Edge(edge)){Collapse_Segment_Edge(i);}
			}break;
			case 3:{
				////triangle edge split
				if(use_triangle_edge_split){
					Vector2i edge_to_split=Vector2i(-1,-1);
					real max_edge_split_length=std::numeric_limits<real>::min();
					for(int j=0;j<edge_n;j++){
						const Vector2i edge=Unique_Ordered(edges[j]);
						if(Check_Split_Triangle_Edge(edge)){real length=Length(edge);
							if(length>max_edge_split_length){edge_to_split=edge;max_edge_split_length=length;}}}
					////split the edge with the longest length among all candidates 
					if(edge_to_split!=Vector2i(-1,-1)){
						Split_Triangle_And_Tetrahedron_Edge(edge_to_split);op_num++;
						if(dbg)Dbg_Check("edge split");continue;}}

				////triangle edge collapse
				if(use_triangle_edge_collapse){
					Vector2i edge_to_collapse=Vector2i(-1,-1);
					real min_edge_collapse_length=std::numeric_limits<real>::max();
					for(int j=0;j<edge_n;j++){
						const Vector2i edge=Unique_Ordered(edges[j]);
						if(Check_Collapse_Triangle_Edge(edge,i)){real length=Length(edge);
							if(length<min_edge_collapse_length){edge_to_collapse=edge;min_edge_collapse_length=length;}}}
					////collapse the edge with the shortest length among all candidates 
					if(edge_to_collapse!=Vector2i(-1,-1)){
						if(Collapse_Triangle_Edge(edge_to_collapse)){op_num++;
							if(dbg)Dbg_Check("edge collapse");continue;}}}

				////triangle edge snap
				if(use_triangle_edge_snap){
					Vector2i edge_to_snap=Vector2i(-1,-1);
					real max_edge_snap_length=std::numeric_limits<real>::min();
					for(int j=0;j<edge_n;j++){const Vector2i edge=Unique_Ordered(edges[j]);
						if(Check_Snap_Triangle_Edge(edge,i)){real length=Length(edge);
							if(length>max_edge_snap_length){edge_to_snap=edge;max_edge_snap_length=length;}}}
					////snap the edge with the longest length among all candidates 
					if(edge_to_snap!=Vector2i(-1,-1)){
						if(Snap_Triangle_Edge(edge_to_snap,i)){op_num++;
							if(dbg)Dbg_Check("edge snap");continue;}}}
			}break;}
		}

		for(int i=0;i<mesh.elements.size();i++){if(!Valid_Element(i))continue;
			const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);if(ds!=3)continue;
			ArrayF<Vector2i,6> edges;int edge_n=Element_Edges(vtx,ds,edges);

			////triangle edge flip
			for(int j=0;j<edge_n;j++){
				const Vector2i edge=Unique_Ordered(edges[j]);
				if(Flip_Triangle_Edge(edge)){op_num++;
					if(dbg)Dbg_Check("edge flip");continue;}}}

		////check if meshing operations converge
		if(abs(op_num-pre_op_num)<min_op_num)break;
		else{pre_op_num=op_num;op_num=0;std::cout<<"op num: "<<pre_op_num<<std::endl;}}

	////Dbg
	//Validate_Edge_Element_Hashtable();
	//Validate_Vertex_Element_Hashtable();
	//Validate_Triangle_Order();
	//Validate_Dihedral_Angles();
}

//////////////////////////////////////////////////////////////////////////
////check meshing operation conditions

template<int d> bool CodimMesher<d>::Check_Split_Segment_Edge(const Vector2i& edge) const
{
	return Length(edge)>length_scale*edge_split_ratio;
}

template<int d> bool CodimMesher<d>::Check_Collapse_Segment_Edge(const Vector2i& edge) const
{return Length(edge)<length_scale*edge_collapse_ratio;}

template<int d> bool CodimMesher<d>::Check_Split_Triangle_Edge(const Vector2i& edge) const
{
	////check angle
	List<int> incident_elements;Value_List(aux.edge_element_hashtable,edge,incident_elements);
	for(auto& e:incident_elements){
		Vector3i vtx=To_V3(E(e));int p=Other_Node_In_Element(vtx,edge[0],edge[1]);
		real angle=MeshFunc::Three_Point_Angle<d>(X(edge[0]),X(p),X(edge[1]));
		if(angle<small_angle)return false;}

	////check length
	return Length(edge)>length_scale*edge_split_ratio;
}

template<int d> bool CodimMesher<d>::Check_Collapse_Triangle_Edge(const Vector2i& edge,const int e) const
{
	////check transition, does not allow an edge vertex connecting to elements that are not triangles
	if(transition){for(int i=0;i<2;i++){
		List<int> vtx_incident_elements;Value_List(aux.vertex_element_hashtable,edge[i],vtx_incident_elements);
		for(const auto& i:vtx_incident_elements)if(mesh.Dim(i)!=3){return false;}}}

	////check length
	bool invalid_edge_length=Length(edge)<length_scale*edge_collapse_ratio;
	if(invalid_edge_length)return true;	

	////check angle
	Vector3i vtx=Triangle_Vertices(e);ArrayF<VectorD,3> vtx_pos;for(int i=0;i<3;i++)vtx_pos[i]=X(vtx[i]);
	ArrayF<real,3> angles;MeshFunc::Triangle_Angles<d>(vtx_pos,angles);
	int p=0;for(int i=0;i<3;i++){if(vtx[i]!=edge[0]&&vtx[i]!=edge[1]){p=i;break;}}
	Vector3i idx(0,1,2);idx=Reordered(idx,p);
	bool invalid_angle=angles[idx[0]]<small_angle&&angles[idx[1]]>small_angle&&angles[idx[2]]>small_angle;

	if(invalid_angle){return true;}
	return false;
}

template<int d> bool CodimMesher<d>::Check_Snap_Triangle_Edge(const Vector2i& edge,const int e) const
{
	////check angle
	Vector3i vtx=To_V3(E(e));int p=Other_Node_In_Element(vtx,edge[0],edge[1]);
	VectorD v0=(X(edge[0])-X(p)).normalized();
	VectorD v1=(X(edge[1])-X(p)).normalized();
	real angle=Angle_Between(v0,v1);return angle>(pi-(real)2*small_angle);
	
	////check length
	//real t;real dis=Distance_From_Point_To_Segment<d>(mesh.Vertices()[p],mesh.Vertices()[edge[0]],mesh.Vertices()[edge[1]],t);
	//return dis<length_scale*edge_snap_ratio;
}

//////////////////////////////////////////////////////////////////////////
////meshing operations
template<int d> void CodimMesher<d>::Split_Segment_Edge(const int i)
{
	int a,c;Element_Vertices(i,a,c);
	int v=Split_Particles(a,c);
	if(verbose)std::cout<<"[Split segment edge]: "<<i<<": "<<a<<", "<<c<<std::endl;

	Remove_Element(i);
	int j=Add_Element(To_V4(a,v));Add_Buffer_Element(j);
	int k=Add_Element(To_V4(v,c));Add_Buffer_Element(k);

	Replace(aux.vertex_element_hashtable,a,i,j);
	Replace(aux.vertex_element_hashtable,c,i,k);
	Add(aux.vertex_element_hashtable,v,j);
	Add(aux.vertex_element_hashtable,v,k);

	Remove(aux.edge_element_hashtable,Unique_Ordered(a,c));
	Add(aux.edge_element_hashtable,Unique_Ordered(a,v),j);
	Add(aux.edge_element_hashtable,Unique_Ordered(v,c),k);
}

template<int d> void CodimMesher<d>::Collapse_Segment_Edge(const int i)
{
	int a,c;Element_Vertices(i,a,c);
	if(verbose)std::cout<<"[Collapse segment edge]: "<<i<<": "<<a<<", "<<c<<std::endl;
	List<int> a_incident_elements;Value_List(aux.vertex_element_hashtable,a,a_incident_elements);a_incident_elements.remove(i);
	List<int> c_incident_elements;Value_List(aux.vertex_element_hashtable,c,c_incident_elements);c_incident_elements.remove(i);
	Hashset<Vector2i> a_incident_edges;Vertex_Incident_Edges(a,a_incident_elements,a_incident_edges);
	Hashset<Vector2i> c_incident_edges;Vertex_Incident_Edges(c,c_incident_elements,c_incident_edges);

	int v=Collapse_Particles(a,c);

	Remove_Element(i);
	for(const auto& i:a_incident_elements){Replace_Vertex<4>(E(i),a,v);Add_Buffer_Element(i);}
	for(const auto& i:c_incident_elements){Replace_Vertex<4>(E(i),c,v);Add_Buffer_Element(i);}

	Remove(aux.vertex_element_hashtable,a,i);
	Remove(aux.vertex_element_hashtable,c,i);
	Merge_Keys(aux.vertex_element_hashtable,a,c,v);

	Remove(aux.edge_element_hashtable,Unique_Ordered(a,c));
	for(const auto& edge:a_incident_edges){Replace_Key(aux.edge_element_hashtable,edge,Unique_Ordered(Replaced<2>(edge,a,v)));}
	for(const auto& edge:c_incident_edges){Replace_Key(aux.edge_element_hashtable,edge,Unique_Ordered(Replaced<2>(edge,c,v)));}
}

template<int d> int CodimMesher<d>::Split_Triangle_And_Tetrahedron_Edge(const Vector2i& edge)
{
	const int a=edge[0];const int c=edge[1];
	int v=Split_Particles(a,c);
	if(verbose)std::cout<<"[Split triangle edge]: "<<a<<", "<<c<<"->"<<v<<std::endl;

	List<int> incident_elements;Value_List(aux.edge_element_hashtable,edge,incident_elements);
	for(const auto& i:incident_elements){const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		switch(ds){
		case 3:{
			const Vector3i& tri_vtx=To_V3(vtx);
			int dd=Other_Node_In_Element(tri_vtx,a,c);
			Vector3i avd,vcd;Split_Triangle_Vertices(tri_vtx,a,c,v,avd,vcd);
			Remove_Element(i);
			int j=Add_Element(To_V4(avd));Add_Buffer_Element(j);
			int k=Add_Element(To_V4(vcd));Add_Buffer_Element(k);

			Replace(aux.vertex_element_hashtable,a,i,j);
			Replace(aux.vertex_element_hashtable,c,i,k);	
			Replace(aux.vertex_element_hashtable,dd,i,j);
			Add(aux.vertex_element_hashtable,dd,k);
			Add(aux.vertex_element_hashtable,v,Array<int>{j,k});

			Replace(aux.edge_element_hashtable,Unique_Ordered(a,dd),i,j);
			Replace(aux.edge_element_hashtable,Unique_Ordered(c,dd),i,k);
			Add(aux.edge_element_hashtable,Unique_Ordered(a,v),j);
			Add(aux.edge_element_hashtable,Unique_Ordered(v,c),k);
			Add(aux.edge_element_hashtable,Unique_Ordered(v,dd),Array<int>{j,k});		
		}break;
		case 4:{
			Vector2i edge2=Other_Nodes_In_Element(vtx,a,c);int b=edge2[0];int dd=edge2[1];
			Remove_Element(i);
			int j=Add_Element(Vector4i(a,b,v,dd));Add_Buffer_Element(j);
			int k=Add_Element(Vector4i(b,c,v,dd));Add_Buffer_Element(k);

			//std::cout<<"abcdv: "<<a<<", "<<b<<", "<<c<<", "<<dd<<", "<<v<<std::endl;
			//std::cout<<"jk: "<<j<<", "<<k<<std::endl;

			Replace(aux.vertex_element_hashtable,a,i,j);
			Replace(aux.vertex_element_hashtable,c,i,k);
			Remove(aux.vertex_element_hashtable,b,i);
			Add(aux.vertex_element_hashtable,b,{j,k});
			Remove(aux.vertex_element_hashtable,dd,i);
			Add(aux.vertex_element_hashtable,dd,{j,k});
			Add(aux.vertex_element_hashtable,v,{j,k});

			Replace(aux.edge_element_hashtable,Unique_Ordered(a,b),i,j);
			Replace(aux.edge_element_hashtable,Unique_Ordered(a,dd),i,j);
			Replace(aux.edge_element_hashtable,Unique_Ordered(b,c),i,k);
			Replace(aux.edge_element_hashtable,Unique_Ordered(dd,c),i,k);
			Remove(aux.edge_element_hashtable,Unique_Ordered(a,c));

			Add(aux.edge_element_hashtable,Unique_Ordered(b,v),{j,k});
			Add(aux.edge_element_hashtable,Unique_Ordered(dd,v),{j,k});
			Add(aux.edge_element_hashtable,Unique_Ordered(a,v),j);
			Add(aux.edge_element_hashtable,Unique_Ordered(c,v),k);
			Remove(aux.edge_element_hashtable,Unique_Ordered(b,dd),i);
			Add(aux.edge_element_hashtable,Unique_Ordered(b,dd),{j,k});
		}break;}
	}
	Remove(aux.edge_element_hashtable,Unique_Ordered(a,c));
	return v;
}

constexpr int ND(int d){return d==2?1:3;}
template<int d> bool CodimMesher<d>::Has_Inverted_Triangle(const int a,const VectorD& new_pos,const List<int>& a_hat_elements) const
{
	bool has_inverted_triangle=false;
	bool has_large_angle=false;
	real cos_epsilon=(real)1e-10;
	real max_cos=(real)1-cos_epsilon;
	for(const int& i:a_hat_elements){
		const Vector3i& vtx=Triangle_Vertices(i);
		ArrayF<VectorD,3> vtx_pos;for(int j=0;j<3;j++)vtx_pos[j]=X(vtx[j]);
		Vector<real,ND(d)> old_normal=Normal(vtx_pos[0],vtx_pos[1],vtx_pos[2]);
		for(int j=0;j<3;j++){if(vtx[j]==a)vtx_pos[j]=new_pos;}
		Vector<real,ND(d)> new_normal=Normal(vtx_pos[0],vtx_pos[1],vtx_pos[2]);
		if(old_normal.dot(new_normal)<cos_epsilon){has_inverted_triangle=true;break;}
	}
	return has_inverted_triangle;
}

template<int d> bool CodimMesher<d>::Collapse_Triangle_Edge(const Vector2i& ac,bool enforce_op/*=false*/)
{
	const int a=ac[0];const int c=ac[1];
	List<int> ac_incident_elements;Value_List(aux.edge_element_hashtable,ac,ac_incident_elements);
	List<int> a_incident_elements;Value_List(aux.vertex_element_hashtable,a,a_incident_elements);
	List<int> c_incident_elements;Value_List(aux.vertex_element_hashtable,c,c_incident_elements);
	List<int> a_hat_elements;Set_Difference(a_incident_elements,ac_incident_elements,a_hat_elements);
	List<int> c_hat_elements;Set_Difference(c_incident_elements,ac_incident_elements,c_hat_elements);

	////merge aware of boundary vertices
	VectorD pos_v=(real).5*(X(a)+X(c));
	if(use_boundary_aware_collapse){
		bool boundary_a=aux.Is_On_Boundary(a);bool boundary_c=aux.Is_On_Boundary(c);
		if(boundary_a!=boundary_c){
			if(boundary_a) pos_v=X(a);
			else pos_v=X(c);}}

	////check inverted elements
	int merged_pos_state=0;	////0-mid, 1-a, 2-c
	real merge_coef=(real).1;
	bool has_inverted_element=Has_Inverted_Triangle(a,pos_v,a_hat_elements)||Has_Inverted_Triangle(c,pos_v,c_hat_elements);
	if(has_inverted_element){
		pos_v=((real)1.-merge_coef)*X(a)+merge_coef*X(c);merged_pos_state=1;
		has_inverted_element=Has_Inverted_Triangle(a,pos_v,a_hat_elements)||Has_Inverted_Triangle(c,pos_v,c_hat_elements);
		if(has_inverted_element){
			pos_v=((real)1.-merge_coef)*X(c)+merge_coef*X(a);merged_pos_state=2;
			has_inverted_element=Has_Inverted_Triangle(a,pos_v,a_hat_elements)||Has_Inverted_Triangle(c,pos_v,c_hat_elements);}}

	////output a warning regardless of verbose or dbg state
	if(enforce_op&&has_inverted_element){std::cerr<<"Error [Codimension] Collapse_Triangle_Edge enforce_op but has_inverted_element"<<std::endl;}
	if(!enforce_op&&has_inverted_element)return false;

	Hashset<Vector2i> a_hat_edges;Vertex_Incident_Edges(a,a_hat_elements,a_hat_edges);
	Hashset<Vector2i> c_hat_edges;Vertex_Incident_Edges(c,c_hat_elements,c_hat_edges);

	int v=Collapse_Particles(a,c);
	if(track_vertex_ancestor)Merge_Keys(aux.vertex_ancestor_hashtable,a,c,v);
	X(v)=pos_v;

	if(verbose)std::cout<<"[Collapse triangle edge]: "<<a<<", "<<c<<"->"<<v<<", state: "<<merged_pos_state<<std::endl;

	Remove(aux.edge_element_hashtable,ac);
	ArrayF<Vector2i,6> edges;
	for(const auto& i:ac_incident_elements){const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		for(int j=0;j<ds;j++)Remove(aux.vertex_element_hashtable,vtx[j],i);
		int edge_n=Element_Edges(vtx,ds,edges);
		for(int j=0;j<edge_n;j++){Vector2i e=Unique_Ordered(edges[j]);Remove(aux.edge_element_hashtable,e,i);}}

	Merge_Keys(aux.vertex_element_hashtable,a,c,v);

	for(const auto& edge:a_hat_edges){Replace_Key(aux.edge_element_hashtable,edge,Unique_Ordered(Replaced<2>(edge,a,v)));}
	for(const auto& edge:c_hat_edges){Replace_Key(aux.edge_element_hashtable,edge,Unique_Ordered(Replaced<2>(edge,c,v)));}

	Hashset<Vector2i> degenerated_segments;
	if(transition){
		for(const auto& i:ac_incident_elements){const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
			switch(ds){
			case 3:{
				Vector3i tri_vtx=To_V3(vtx);
				int b=Other_Node_In_Element(tri_vtx,a,c);
				degenerated_segments.insert(Unique_Ordered(Vector2i(b,v)));
			}break;}}}

	for(const auto& i:ac_incident_elements)Remove_Element(i);
	for(const auto& i:a_hat_elements){Replace_Vertex<4>(E(i),a,v);Add_Buffer_Element(i);}
	for(const auto& i:c_hat_elements){Replace_Vertex<4>(E(i),c,v);Add_Buffer_Element(i);}

	////triangles degenerate to segments
	if(transition){
		for(const auto& edge:degenerated_segments){
			if(Value_Count(aux.edge_element_hashtable,edge)==0){
				Vector4i vtx=To_V4(edge);
				int i=Add_Element(vtx);
				aux.Add_New_Element_To_Hashtables(vtx,i,2);}}}

	return true;
}

template<int d> bool CodimMesher<d>::Flip_Improve_Triangle_Quality(const Vector2i& ac,const int b,const int g,int& e_i,int& e_j) const
{
	const int a=ac[0];const int c=ac[1];

	bool flip=false;
	Array<Vector3i> old_tris;old_tris.push_back(Vector3i(a,c,g));old_tris.push_back(Vector3i(a,b,c));
	Array<Vector3i> new_tris;new_tris.push_back(Vector3i(a,b,g));new_tris.push_back(Vector3i(c,g,b));
	real old_measure=Measure_Triangle_Quality(old_tris);real new_measure=Measure_Triangle_Quality(new_tris);
	real measure_coef=(real).98;if(old_measure<measure_coef*new_measure){flip=true;}
	return flip;
}

template<int d> bool CodimMesher<d>::Flip_Triangle_Edge(const Vector2i& ac,bool enforce_op/*=false*/)
{
	const int a=ac[0];const int c=ac[1];
	int b,dd,e_i,e_j;bool find_tri_pair=Edge_Incident_Triangle_Pair(ac,b,dd,e_i,e_j);
	if(!find_tri_pair)return false;

	if(!enforce_op){
		if(!Flip_Improve_Triangle_Quality(ac,b,dd,e_i,e_j))return false;
		real dihedral_angle=Dihedral_Angle(mesh.Vertices()[b],mesh.Vertices()[dd],mesh.Vertices()[c],mesh.Vertices()[a]);
		if(dihedral_angle<edge_flip_dihedral_angle||dihedral_angle>two_pi-edge_flip_dihedral_angle)return false;}

	if(verbose)std::cout<<"[Flip triangle edge]: "<<a<<", "<<c<<", "<<b<<", "<<dd<<": "<<e_i<<", "<<e_j<<std::endl;

	Set_Element(e_i,To_V4(a,b,dd));Add_Buffer_Element(e_i);
	Set_Element(e_j,To_V4(c,dd,b));Add_Buffer_Element(e_j);

	Remove(aux.vertex_element_hashtable,a,e_j);
	Remove(aux.vertex_element_hashtable,c,e_i);
	Add(aux.vertex_element_hashtable,dd,e_j);
	Add(aux.vertex_element_hashtable,b,e_i);

	Remove(aux.edge_element_hashtable,ac);
	Add(aux.edge_element_hashtable,Unique_Ordered(b,dd),Array<int>{e_i,e_j});
	Replace(aux.edge_element_hashtable,Unique_Ordered(a,b),e_j,e_i);
	Replace(aux.edge_element_hashtable,Unique_Ordered(c,dd),e_i,e_j);

	return true;
}

template<int d> bool CodimMesher<d>::Snap_Triangle_Edge(const Vector2i& ac,const int e)
{
	Vector3i vtx=To_V3(E(e));int p=Other_Node_In_Element(vtx,ac[0],ac[1]);
	int v=Split_Triangle_And_Tetrahedron_Edge(ac);
	if(verbose)std::cout<<"[Snap triangle edge]: "<<ac.transpose()<<", "<<v<<", "<<p<<std::endl;
	return Collapse_Triangle_Edge(Unique_Ordered(Vector2i(p,v))/*,true*/);
}

template<int d> int CodimMesher<d>::Add_Triangle(const ArrayF<VectorD,3>& vtx_pos)
{
	Vector4i vtx={-1,-1,-1,-1};for(int i=0;i<3;i++){int p=Add_Particle();vtx[i]=p;particles.X(p)=vtx_pos[i];}
	int i=Add_Element(vtx);aux.Add_New_Element_To_Hashtables(vtx,i,3);return i;
}

template<int d> bool CodimMesher<d>::Remove_Triangle(const int e)
{
	const Vector4i vtx=E(e);Remove_Element(e);aux.Remove_Element_From_Hashtables(e,3);
	for(int i=0;i<3;i++)if(Value_Count(aux.vertex_element_hashtable,vtx[i])==0)Remove_Particle(vtx[i]);return true;
}

template<int d> int CodimMesher<d>::Add_Segment(const ArrayF<VectorD,2>& vtx_pos)
{
	Vector4i vtx={-1,-1,-1,-1};for(int i=0;i<2;i++){int p=Add_Particle();vtx[i]=p;particles.X(p)=vtx_pos[i];}
	int i=Add_Element(vtx);aux.Add_New_Element_To_Hashtables(vtx,i,2);return i;	
}

template<int d> bool CodimMesher<d>::Remove_Segment(const int e)
{
	const Vector4i vtx=E(e);Remove_Element(e);aux.Remove_Element_From_Hashtables(e,2);
	for(int i=0;i<2;i++)if(Value_Count(aux.vertex_element_hashtable,vtx[i])==0)Remove_Particle(vtx[i]);return true;
}

//////////////////////////////////////////////////////////////////////////
////helper functions
template<int d> int CodimMesher<d>::Add_Particle()
{
	if(!invalid_particle_heap.empty()){
		int p=invalid_particle_heap.top();invalid_particle_heap.pop();
		invalid_particle_hashset.erase(p);
		//particles.Id(p)=particles.Id_Generator();
		return p;}
	else return Append_Particle();
}

template<int d> void CodimMesher<d>::Pruned_Copy(CodimMesh<d>& pruned_mesh,CodimParticles<d>& pruned_particles,std::shared_ptr<Array<int> >* _remap)
{
	pruned_mesh.elements.clear();
	for(int i=0;i<(int)mesh.elements.size();i++)
		if(Valid_Element(i)){pruned_mesh.elements.push_back(E(i));}
	
	int pn=particles.Size()-(int)invalid_particle_hashset.size();pruned_particles.Resize(pn);
	std::shared_ptr<Array<int> > remap=std::make_shared<Array<int> >(particles.Size(),-1);
	int j=0;for(int i=0;i<particles.Size();i++)
		if(Valid_Particle(i)){pruned_particles.Copy_Attribute_From(j,particles,i);(*remap)[i]=j;j++;}

	for(auto& e:pruned_mesh.elements){const int ds=CodimMesh<d>::Dim(e);
		for(int i=0;i<ds;i++)e[i]=(*remap)[e[i]];for(int i=ds;i<4;i++)e[i]=-1;}

	if(_remap!=nullptr)(*_remap)=remap;
}

template<int d> void CodimMesher<d>::Prune(std::shared_ptr<Array<int> >* _remap)
{
	CodimMesh<d> pruned_mesh;CodimParticles<d> pruned_particles;
	Pruned_Copy(pruned_mesh,pruned_particles,_remap);
	mesh=pruned_mesh;particles=pruned_particles;
}

template<int d> void CodimMesher<d>::Vertex_Incident_Edges(const int v,const List<int>& incident_elements,Hashset<Vector2i>& edge_hashset)
{
	ArrayF<Vector2i,6> edges;
	for(const int& i:incident_elements){const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		int edge_n=Element_Vertices_Edges(vtx,ds,v,edges);
		for(int j=0;j<edge_n;j++)edge_hashset.insert(Unique_Ordered(edges[j]));}
}

template<int d> void CodimMesher<d>::Vertex_Incident_Edges(const int v,Hashset<Vector2i>& edge_hashset)
{
	List<int> v_incident_elements;Value_List(aux.vertex_element_hashtable,v,v_incident_elements);
	Vertex_Incident_Edges(v,v_incident_elements,edge_hashset);
}

template<int d> bool CodimMesher<d>::Edge_Incident_Triangle_Pair(const Vector2i& ac,int& b,int& g,int& e_i,int& e_j) const
{
	const int a=ac[0];const int c=ac[1];Array<int> ac_incident_elements;
	Value_Array(aux.edge_element_hashtable,ac,ac_incident_elements);
	if(ac_incident_elements.size()!=2)return false;
	e_i=ac_incident_elements[0];e_j=ac_incident_elements[1];
	if(CodimMesh<2>::Dim(E(e_i))!=3||CodimMesh<2>::Dim(E(e_j))!=3)return false;

	{Vector3i vtx=Reordered(To_V3(E(e_i)),a);if(vtx[2]==c){int tmp=e_i;e_i=e_j;e_j=tmp;}}	////swap e_i and e_j if e_i is not acd
	const Vector3i acd=To_V3(E(e_i));g=Other_Node_In_Element(acd,a,c);
	const Vector3i abc=To_V3(E(e_j));b=Other_Node_In_Element(abc,a,c);
	return true;
}

template<int d> int CodimMesher<d>::Edge_Type(const List<int>& incident_elements) const
{
	int type=0;
	for(const auto& i:incident_elements){const int ds=mesh.Dim(i);
		switch(ds){
		case 2:Set_On_Segment_Mesh(type);
		case 3:Set_On_Triangle_Mesh(type);
		case 4:Set_On_Tetrahedron_Mesh(type);}}
	return type;
}

template<int d> int CodimMesher<d>::Edge_Type(const Vector2i& edge) const
{
	List<int> incident_elements;Value_List(aux.edge_element_hashtable,edge,incident_elements);
	return Edge_Type(incident_elements);
}

template<int d> real CodimMesher<d>::Measure_Triangle_Quality(const Array<Vector3i>& tris_v) const
{
	real measure=std::numeric_limits<real>::max();
	for(const auto& v:tris_v) measure=(std::min)(measure,Area_Length_Measure<d>(X(v[0]),X(v[1]),X(v[2])));
	return measure;
}

template<int d> void CodimMesher<d>::Coarsen(const real coarsen_factor)
{
	length_scale*=coarsen_factor;
	track_vertex_ancestor=true;
	Update_Vertex_Ancestor_Hashtable();
	Update();
	Validate_Coarsening_Mapping();
}

template<int d> void CodimMesher<d>::Build_Hierarchy(const int levels,const real coarsen_factor,Array<CodimMesh<d> >& meshes/*,Array<SparseMatrixT>& mappings*/)
{
	meshes.push_back(mesh);
	for(int i=0;i<levels;i++){
		CodimMesh<d> fine_pruned_mesh;CodimParticles<d> fine_pruned_particles;
		std::shared_ptr<Array<int> > fine_remap=std::make_shared<Array<int> >();
		Pruned_Copy(fine_pruned_mesh,fine_pruned_particles,&fine_remap);
		Coarsen(coarsen_factor);
		CodimMesh<d> coarse_pruned_mesh;CodimParticles<d> coarse_pruned_particles;
		std::shared_ptr<Array<int> > coarse_remap=std::make_shared<Array<int> >();
		Pruned_Copy(coarse_pruned_mesh,coarse_pruned_particles,&coarse_remap);
		//Array<TripletT> elements;
		for(int i=0;i<(*coarse_remap).size();i++){if((*coarse_remap)[i]==-1)continue;
			List<int> anc;Value_List(aux.vertex_ancestor_hashtable,i,anc);
			std::cout<<(*coarse_remap)[i]<<": ";
			for(auto& a:anc){std::cout<<(*fine_remap)[a]<<", ";}std::cout<<std::endl;}}	
}

template<int d> bool CodimMesher<d>::Validate_Coarsening_Mapping() const
{
	for(int i=0;i<particles.Size();i++){if(!Valid_Particle(i))continue;
		List<int> anc;Value_List(aux.vertex_ancestor_hashtable,i,anc);
		std::cout<<i<<": ";for(auto& a:anc)std::cout<<a<<", ";std::cout<<std::endl;}
	return true;
}

void Split_Dimension(const CodimMesh<3>& mesh,SegmentMesh<3>& seg_mesh,TriangleMesh<3>& tri_mesh,TetrahedronMesh<3>& tet_mesh)
{
	seg_mesh.vertices=mesh.vertices;
	tri_mesh.vertices=mesh.vertices;
	tet_mesh.vertices=mesh.vertices;
	for(int i=0;i<(int)mesh.elements.size();i++){
		const auto& vtx=mesh.elements[i];int ds=CodimMesh<3>::Dim(vtx);
		switch(ds){
			case 2:seg_mesh.elements.push_back(To_V<2>(vtx));break;
			case 3:tri_mesh.elements.push_back(To_V<3>(vtx));break;
			case 4:tet_mesh.elements.push_back(vtx);break;}}
}

void Split_Dimension(const CodimMesh<2>& mesh,SegmentMesh<2>& seg_mesh,TriangleMesh<2>& tri_mesh)
{
	seg_mesh.vertices=mesh.vertices;
	tri_mesh.vertices=mesh.vertices;
	for(int i=0;i<(int)mesh.elements.size();i++){
		const auto& vtx=mesh.elements[i];int ds=CodimMesh<2>::Dim(vtx);
		switch(ds){
			case 2:seg_mesh.elements.push_back(To_V<2>(vtx));break;
			case 3:tri_mesh.elements.push_back(To_V<3>(vtx));break;}}
}

void Merge_Dimension(CodimMesh<3>& mesh,const SegmentMesh<3>* seg_mesh/*=nullptr*/,const TriangleMesh<3>* tri_mesh/*=nullptr*/,const TetrahedronMesh<3>* tet_mesh/*=nullptr*/)
{
	mesh.Vertices()=(seg_mesh!=nullptr?seg_mesh->Vertices():(tri_mesh!=nullptr?tri_mesh->Vertices():tet_mesh->Vertices()));	////assuming simplex meshes share vertices
	if(seg_mesh!=nullptr)for(int i=0;i<(int)seg_mesh->elements.size();i++){
		mesh.elements.push_back(To_V4(seg_mesh->elements[i]));}
	if(tri_mesh!=nullptr)for(int i=0;i<(int)tri_mesh->elements.size();i++){
		mesh.elements.push_back(To_V4(tri_mesh->elements[i]));}
	if(tet_mesh!=nullptr)for(int i=0;i<(int)tet_mesh->elements.size();i++){
		mesh.elements.push_back(tet_mesh->elements[i]);}
}

void Merge_Dimension(CodimMesh<2>& mesh,const SegmentMesh<2>* seg_mesh/*=nullptr*/,const TriangleMesh<2>* tri_mesh/*=nullptr*/,const TetrahedronMesh<2>* tet_mesh/*=nullptr*/)
{
	mesh.Vertices()=(seg_mesh!=nullptr?seg_mesh->Vertices():tri_mesh->Vertices());	////assuming simplex meshes share vertices
	if(seg_mesh!=nullptr)for(int i=0;i<(int)seg_mesh->elements.size();i++){
		mesh.elements.push_back(To_V4(seg_mesh->elements[i]));}
	if(tri_mesh!=nullptr)for(int i=0;i<(int)tri_mesh->elements.size();i++){
		mesh.elements.push_back(To_V4(tri_mesh->elements[i]));}
}

////TOFIX
////CodimMeshAux
//template<int d2> bool On_Dimension(const int v) const;
//{
//	List<int> vertex_elements;Value_List(vertex_element_hashtable,v,vertex_elements);
//	for(const int& i:vertex_elements){const int ds=mesh.Dim(i);if(ds==d2)return true;}return false;
//}

//template<int d2> bool On_Dimension(const Vector2i& edge) const	////assuming sorted
//{
//	List<int> edge_elements;Value_List(edge_element_hashtable,edge,edge_elements);
//	for(const int& i:edge_elements){const int ds=mesh.Dim(i);if(ds==d2)return true;}return false;
//}

//template<int ds,int dg_ds> bool Degenerate_By_Edge(const Vector4i& vtx,const Vector2i& ac,const int v,Vector4i& dg_vtx){return false;}
//template<> bool Degenerate_By_Edge<3,2>(const Vector4i& vtx,const Vector2i& ac,const int v,Vector4i& dg_vtx)	////tri->seg
//{
//	const int a=ac[0];const int c=ac[1];
//	int d=Other_Node_In_Element(To_V3(vtx),a,c);
//	const Vector2i ad=Unique_Ordered(a,d);
//	const Vector2i cd=Unique_Ordered(c,d);
//	if(aux.On_Dimension_Boundary_3(ad)&&aux.On_Dimension_Boundary_3(cd)){dg_vtx=To_V4(d,v);return true;}
//	else{return false;}
//}
//template<> bool Degenerate_By_Edge<4,3>(const Vector4i& vtx,const Vector2i& ac,const int v,Vector4i& dg_vtx)	////tet->tri
//{return false;}	////TODO

//////////////////////////////////////////////////////////////////////////
////Debug functions
template<int d> void CodimMesher<d>::Dbg_Output() const
{
	//for(int i=0;i<particles.Size();i++)particles.Boundary(i)=aux.Boundary_Type(i);
	std::cout<<"particles: "<<particles.Size()<<std::endl;  
	//for(int i=0;i<particles.Size();i++)std::cout<<i<<", X: "<<X(i).transpose()<<", bdry: "<<particles.Boundary(i)<<", valid: "<<Valid_Particle(i)<<std::endl;
	std::cout<<"elements: "<<mesh.elements.size()<<std::endl;
	for(int i=0;i<(int)mesh.elements.size();i++)std::cout<<i<<", vtx: "<<E(i).transpose()<<", valid: "<<Valid_Element(i)<<std::endl;

	std::cout<<"vtx_ele_hashtable: "<<aux.vertex_element_hashtable.size()<<std::endl;
	for(const auto iter:aux.vertex_element_hashtable)std::cout<<iter.first<<": "<<iter.second<<std::endl;
	std::cout<<"edge_ele_hashtable: "<<aux.edge_element_hashtable.size()<<std::endl;
	for(const auto iter:aux.edge_element_hashtable)std::cout<<iter.first.transpose()<<": "<<iter.second<<std::endl;
}

template<int d> bool CodimMesher<d>::Dbg_Check(const std::string& info) const
{
	bool valid_dihedral_angle=Validate_Dihedral_Angles();
	bool valid_normal_2d=(d==2?true:Validate_Normal_2D());
	if(!valid_dihedral_angle||!valid_normal_2d){
		std::cerr<<"-------------------- Invalid case for "<<info<<" --------------------"<<std::endl;return false;}
	return true;
}

template<int d> bool CodimMesher<d>::Validate() const
{
	bool valid=true;
	valid=valid&&Validate_Edge_Element_Hashtable();
	valid=valid&&Validate_Vertex_Element_Hashtable();
	if(d==2)valid=valid&&Validate_Triangle_Order();
	else valid=valid&&Validate_Tetrahedron_Order();
	valid=valid&&Validate_Dihedral_Angles();
	if(!valid){std::cout<<"Error occurs in CodimMesher::Validate()"<<std::endl;}
	return valid;
}

template<int d> bool CodimMesher<d>::Validate_Triangle_Order() const
{
	bool valid=true;
	for(int i=0;i<mesh.elements.size();i++){if(!Valid_Element(i))continue;
		const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		ArrayF<Vector2i,6> edges;int edge_n=Element_Edges(vtx,ds,edges);
		for(int j=0;j<edge_n;j++){Vector2i ac=edges[j];int a=ac[0];int c=ac[1];	
			int b,dd,e_i,e_j;bool find_tri_pair=Edge_Incident_Triangle_Pair(ac,b,dd,e_i,e_j);
			if(!find_tri_pair)continue;
			Vector3i vtx_i=To_V3(E(e_i));vtx_i=Reordered(vtx_i,a);
			Vector3i vtx_j=To_V3(E(e_j));vtx_j=Reordered(vtx_j,a);
			if(!(vtx_i[1]==c&&vtx_j[1]!=c||vtx_i[1]!=c&&vtx_j[1]==c)){valid=false;
				std::cerr<<"Error: [Codimension] Invalid edge-triangle order: "<<a<<", "<<c<<": "<<b<<", "<<d<<": "<<vtx_i.transpose()<<", "<<vtx_j.transpose()<<std::endl;}}}
	return valid;
}

template<int d> bool CodimMesher<d>::Validate_Tetrahedron_Order() const
{
	bool valid=true;
	for(int i=0;i<mesh.elements.size();i++){if(!Valid_Element(i))continue;
		const Vector4i& v=E(i);
		real vol=MeshFunc::Tetrahedron_Volume(X(v[0]),X(v[1]),X(v[2]),X(v[3]));
		if(vol<(real)0){valid=false;
			std::cerr<<"Error: [Codimension] Invalid tetrahedron: "<<i<<", "<<v.transpose()<<std::endl;}}
	return valid;
}

template<int d> bool CodimMesher<d>::Validate_Vertex_Element_Hashtable() const
{
	bool valid=true;
	for(auto& iter:aux.vertex_element_hashtable){const int a=iter.first;
		List<int> incident_elements;Value_List(aux.vertex_element_hashtable,a,incident_elements);
		for(auto& i:incident_elements){
			const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
			switch(ds){
			case 2:{
				Vector2i seg_vtx=To_V2(vtx);
				if(!AuxFunc::Has(seg_vtx,a)){valid=false;
					std::cerr<<"Error: [Codimension] Invalid vertex-segment in vertex_element_hashtable: "<<a<<": "<<vtx.transpose()<<std::endl;}
			}break;
			case 3:{
				Vector3i tri_vtx=To_V3(vtx);
				if(!AuxFunc::Has(tri_vtx,a)){valid=false;
					std::cerr<<"Error: [Codimension] Invalid vertex-triangle in vertex_element_hashtable: "<<a<<": "<<vtx.transpose()<<std::endl;}
			}break;}}}
	return valid;
}

template<int d> bool CodimMesher<d>::Validate_Edge_Element_Hashtable() const
{
	bool valid=true;
	for(auto& iter:aux.edge_element_hashtable){Vector2i edge=iter.first;
		List<int> incident_elements;Value_List(aux.edge_element_hashtable,edge,incident_elements);
		for(auto& i:incident_elements){
			const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
			switch(ds){
			case 2:{
				Vector2i seg_vtx=To_V2(vtx);
				if(!AuxFunc::Has(seg_vtx,edge[0])||!AuxFunc::Has(seg_vtx,edge[1])){valid=false;
					std::cerr<<"Error: [Codimension] Invalid edge-segment in edge_element_hashtable: "<<edge.transpose()<<": "<<vtx.transpose()<<std::endl;}
			}break;
			case 3:{
				Vector3i tri_vtx=To_V3(vtx);
				if(!AuxFunc::Has(tri_vtx,edge[0])||!AuxFunc::Has(tri_vtx,edge[1])){valid=false;
					std::cerr<<"Error: [Codimension] Invalid edge-triangle in edge_element_hashtable: "<<edge.transpose()<<": "<<vtx.transpose()<<std::endl;}
			}break;}}}
	return valid;
}

template<int d> bool CodimMesher<d>::Validate_Dihedral_Angles() const
{
	real threshold=(real).1*pi;bool valid=true;
	for(int i=0;i<mesh.elements.size();i++){if(!Valid_Element(i))continue;
		const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		ArrayF<Vector2i,6> edges;int edge_n=Element_Edges(vtx,ds,edges);
		for(int j=0;j<edge_n;j++){Vector2i ac=edges[j];int a=ac[0];int c=ac[1];	
			int b,dd,e_i,e_j;bool find_tri_pair=Edge_Incident_Triangle_Pair(ac,b,dd,e_i,e_j);
			if(!find_tri_pair)continue;
			real angle=Dihedral_Angle(X(a),X(c),X(b),X(d));
			if(angle<threshold||angle>two_pi-threshold){valid=false;
				std::cerr<<"Error: [Codimension] Invalid dihedral angle: "<<a<<", "<<c<<", "<<b<<", "<<d<<": "<<angle<<std::endl;}}}
	return valid;
}

////for 2D only
template<int d> bool CodimMesher<d>::Validate_Normal_2D() const
{
	bool valid=true;
	for(int i=0;i<mesh.elements.size();i++){if(!Valid_Element(i))continue;
		const Vector4i& vtx=E(i);const int ds=mesh.Dim(i);
		switch(ds){
		case 3:{
			Vector<real,ND(d)> normal=Normal(X(vtx[0]),X(vtx[1]),X(vtx[2]));
			if(normal[0]<0.){valid=false;std::cerr<<"Error: [Codimension] Invalid normal for element: "<<vtx.transpose()<<", "<<normal[0]<<std::endl;}
		}break;}}
	return valid;
}

template class CodimMesher<2>;
template class CodimMesher<3>;

template<int d> void CodimMeshElementIterator<d>::Next()
{
	do{index++;}while(!mesher->Valid_Element(index)&&Valid());
}

template class CodimMeshElementIterator<2>;
template class CodimMeshElementIterator<3>;

};
//////////////////////////////////////////////////////////////////////////
// Mesh functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <numeric>
#include <iostream>
#include "GeometryPrimitives.h"
#include "Constants.h"
#include "Hashtable.h"
#include "Grid.h"
#include "Mesh.h"
#include "MeshFunc.h"
#ifdef USE_TRI2D
#include "Triangulation2D.h"	
#endif

namespace MeshFunc{
    using namespace AuxFunc;

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Normals
    Vector2 Normal(const Vector2& v1,const Vector2& v2)
	{Vector2 v12=(v2-v1).normalized();return Vector2(v12[1],-v12[0]);}
    
	Vector1 Normal(const Vector2& p1,const Vector2& p2,const Vector2& p3)
	{Vector2 p12=p2-p1;Vector2 p13=p3-p1;return Cross(p12,p13).normalized();}
    
	Vector3 Normal(const Vector3& p1,const Vector3& p2,const Vector3& p3)
	{return (p2-p1).cross(p3-p1).normalized();}

    template<> Vector<real,2> Normal<2>(const ArrayF<Vector<real,2>,2>& v){return Normal(v[0],v[1]);}
    template<> Vector<real,3> Normal<3>(const ArrayF<Vector<real,3>,3>& v){return Normal(v[0],v[1],v[2]);}

    template<> void Update_Normals<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,Array<Vector2>& normals)
    {
        normals.resize(vertices.size(),Vector2::Zero());
        for(const auto& v:elements){Vector2 n=Normal(vertices[v[0]],vertices[v[1]]);for(int j=0;j<2;j++){normals[v[j]]+=n;}}
        for(auto& n:normals){n.normalize();}
    }
    template<> void Update_Normals<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,Array<Vector3>& normals)
    {
        normals.resize(vertices.size(),Vector3::Zero());
        for(const auto& v:elements){Vector3 n=Normal(vertices[v[0]],vertices[v[1]],vertices[v[2]]);for(int j=0;j<3;j++){normals[v[j]]+=n;}}
        for(auto& n:normals){n.normalize();}
    }

    template<> Vector2 Element_Normal<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,const int i)
    {return Normal(vertices[elements[i][0]],vertices[elements[i][1]]);}
    
	template<> Vector3 Element_Normal<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,const int i)
    {return Normal(vertices[elements[i][0]],vertices[elements[i][1]],vertices[elements[i][2]]);}

    template<> Vector2 Element_Area_Weighted_Normal<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,const int i)
    {return Area_Weighted_Normal(vertices[elements[i][0]],vertices[elements[i][1]]);}
    
	template<> Vector3 Element_Area_Weighted_Normal<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,const int i)
    {return Area_Weighted_Normal(vertices[elements[i][0]],vertices[elements[i][1]],vertices[elements[i][2]]);}

	Vector2 Area_Weighted_Normal(const Vector2& v1,const Vector2& v2)
	{Vector2 v12=(v2-v1).normalized();return Vector2(v12[1],-v12[0]);}

	Vector1 Area_Weighted_Normal(const Vector2& p1,const Vector2& p2,const Vector2& p3)
	{Vector2 p12=p2-p1;Vector2 p13=p3-p1;return Cross(p12,p13);}

	Vector3 Area_Weighted_Normal(const Vector3& p1,const Vector3& p2,const Vector3& p3)
	{return (p2-p1).cross(p3-p1);}

	////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Element size and center
    real Tetrahedron_Volume(const Vector3& v0,const Vector3& v1,const Vector3& v2,const Vector3& v3){return (real)one_sixth*abs((v1-v0).cross(v2-v0).dot(v3-v0));}
    real Tetrahedron_Volume(const ArrayF<Vector3,4>& tet){return Tetrahedron_Volume(tet[0],tet[1],tet[2],tet[3]);}
	real Tetrahedron_Volume(const Vector2&,const Vector2&,const Vector2&,const Vector2&){return (real)0;}

	real Triangle_Area(const Vector2& v0,const Vector2& v1,const Vector2& v2){Vector2 v01=v1-v0;Vector2 v02=v2-v0;return (real).5*abs(Cross(v01,v02)[0]);}
    real Triangle_Area(const ArrayF<Vector2,3>& tri){return Triangle_Area(tri[0],tri[1],tri[2]);}
    real Triangle_Area(const Vector3& v0,const Vector3& v1,const Vector3& v2){return (real).5*((v1-v0).cross(v2-v0)).norm();}
    real Triangle_Area(const ArrayF<Vector3,3>& tri){return Triangle_Area(tri[0],tri[1],tri[2]);}
	
	real Simplex_Size(const ArrayF<Vector2,2>& seg){return (seg[0]-seg[1]).norm();}
	real Simplex_Size(const ArrayF<Vector3,2>& seg){return (seg[0]-seg[1]).norm();}
	real Simplex_Size(const ArrayF<Vector2,3>& tri){return Triangle_Area(tri);}
	real Simplex_Size(const ArrayF<Vector3,4>& tet){return Tetrahedron_Volume(tet);}

    template<> real Element_Size<2>(Array<Vector2>& v){return (v[1]-v[0]).norm();}
    template<> real Element_Size<3>(Array<Vector3>& v){return (real).5*((v[1]-v[0]).cross(v[2]-v[0])).norm();}

    template<> real Element_Size<2,2>(const Array<Vector2>& vertices,const Array<Vector2i>& elements,const int i)
    {return (vertices[elements[i][1]]-vertices[elements[i][0]]).norm();}
    template<> real Element_Size<3,3>(const Array<Vector3>& vertices,const Array<Vector3i>& elements,const int i)
    {return (real).5*((vertices[elements[i][1]]-vertices[elements[i][0]]).cross(vertices[elements[i][2]]-vertices[elements[i][0]])).norm();}

	template<int d> Vector<real,d> Simplex_Center(const ArrayF<Vector<real,d>,d+1>& v)
	{
		Vector<real,d> center=Vector<real,d>::Zero();
		for(int i=0;i<v.size();i++){center+=v[i];}center/=(real)(v.size());return center;		
	}
	template Vector2 Simplex_Center<2>(const ArrayF<Vector2,3>&);
	template Vector3 Simplex_Center<3>(const ArrayF<Vector3,4>&);

    template<> Vector<real,2> Element_Center<2>(const Array<Vector2>& v){return (real).5*(v[0]+v[1]);}
    template<> Vector<real,3> Element_Center<3>(const Array<Vector3>& v){return (real)one_third*(v[0]+v[1]+v[2]);}

    template<int d,int e_d> Vector<real,d> Element_Center(const Array<Vector<real,d> >& vertices,const Array<Vector<int,e_d> >& elements,const int i)
    {Vector<real,d> center=Vector<real,d>::Zero();for(int j=0;j<e_d;j++)center+=vertices[elements[i][j]];center/=(real)e_d;return center;}

#define comma ,
#define Inst_Helper(d,e_d) \
template Vector<real,d> Element_Center<d,e_d>(const Array<Vector<real,d> >&,const Array<Vector<int,e_d> >&,const int)
Inst_Helper(2,2);Inst_Helper(2,3);Inst_Helper(3,3);Inst_Helper(3,4);
#undef Inst_Helper
#undef comma

	template<int d> void Triangle_Angles(const ArrayF<Vector<real,d>,3>& tri,ArrayF<real,3>& angles)
	{
		Vector<real,d> v01=(tri[1]-tri[0]).normalized();Vector<real,d> v02=(tri[2]-tri[0]).normalized();angles[0]=Angle_Between(v01,v02);
		Vector<real,d> v12=(tri[2]-tri[1]).normalized();angles[1]=Angle_Between(v12,-v01);
		angles[2]=pi-angles[0]-angles[1];	
	}
	template void Triangle_Angles<2>(const ArrayF<Vector2,3>&,ArrayF<real,3>&);
	template void Triangle_Angles<3>(const ArrayF<Vector3,3>&,ArrayF<real,3>&);
	
	template<int d> real Three_Point_Angle(const Vector<real,d>& a,const Vector<real,d>& b,const Vector<real,d>& c)
	{
		Vector<real,d> ba=(a-b).normalized();Vector<real,d> bc=(c-b).normalized();return Angle_Between(bc,ba);
	}
	template real Three_Point_Angle<2>(const Vector2&,const Vector2&,const Vector2&);
	template real Three_Point_Angle<3>(const Vector3&,const Vector3&,const Vector3&);

    ////////////////////////////////////////////////////////////////////////////////////////////////////
    ////Element edges
    template<class T_ARRAY> int Element_Edges(const Vector2i& v,T_ARRAY& edges)
    {edges[0]=v;return 1;}
    template<class T_ARRAY> int Element_Edges(const Vector3i& v,T_ARRAY& edges)
    {edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[1],v[2]);edges[2]=Vector2i(v[2],v[0]);return 3;}
    template<class T_ARRAY> int Element_Edges(const Vector4i& v,T_ARRAY& edges)
    {edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[0],v[2]);edges[2]=Vector2i(v[0],v[3]);edges[3]=Vector2i(v[1],v[2]);edges[4]=Vector2i(v[2],v[3]);edges[5]=Vector2i(v[3],v[1]);return 6;}
    template<class T_ARRAY> int Element_Edges(const Vector4i& v,const int ds,T_ARRAY& edges)
    {switch(ds){case 2:return Element_Edges(Vector2i(v[0],v[1]),edges);case 3:return Element_Edges(Vector3i(v[0],v[1],v[2]),edges);case 4:return Element_Edges(v,edges);default:return 0;}}

    template<class T_ARRAY> int Quad_Edges(const Vector4i& v,T_ARRAY& edges)
    {edges[0]=Vector2i(v[0],v[1]);edges[1]=Vector2i(v[1],v[2]);edges[2]=Vector2i(v[2],v[3]);edges[3]=Vector2i(v[3],v[0]);return 4;}

    template<class T_ARRAY> int Element_Faces(const Vector3i& v,T_ARRAY& faces)
    {faces[0]=v;return 1;}
    template<class T_ARRAY> int Element_Faces(const Vector4i& v,T_ARRAY& faces)
    {faces[0]=Vector3i(v[0],v[1],v[3]);faces[1]=Vector3i(v[1],v[2],v[3]);faces[2]=Vector3i(v[2],v[0],v[3]);faces[3]=Vector3i(v[1],v[0],v[2]);return 4;}
    template<class T_ARRAY> int Element_Faces(const Vector4i& v,const int ds,T_ARRAY& faces)
    {switch(ds){case 3:return Element_Faces(Vector3i(v[0],v[1],v[2]),faces);case 4:return Element_Faces(v,faces);default:return 0;}}

#define comma ,
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Edges<T_ARRAY >(const T_VEC& v,T_ARRAY& edges)
    Inst_Helper(Vector2i,Array<Vector2i>);Inst_Helper(Vector2i,ArrayF<Vector2i comma 1>);Inst_Helper(Vector2i,ArrayF<Vector2i comma 6>);
    Inst_Helper(Vector3i,Array<Vector2i>);Inst_Helper(Vector3i,ArrayF<Vector2i comma 3>);Inst_Helper(Vector3i,ArrayF<Vector2i comma 6>);
    Inst_Helper(Vector4i,Array<Vector2i>);Inst_Helper(Vector4i,ArrayF<Vector2i comma 6>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Edges<T_ARRAY >(const T_VEC& v,const int ds,T_ARRAY& edges)
    Inst_Helper(Vector4i,Array<Vector2i>);Inst_Helper(Vector4i,ArrayF<Vector2i comma 6>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Quad_Edges<T_ARRAY >(const T_VEC& v,T_ARRAY& edges)
    Inst_Helper(Vector4i,Array<Vector2i>);Inst_Helper(Vector4i,ArrayF<Vector2i comma 4>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Faces<T_ARRAY >(const T_VEC& v,T_ARRAY& faces)
    Inst_Helper(Vector3i,Array<Vector3i>);Inst_Helper(Vector3i,ArrayF<Vector3i comma 1>);Inst_Helper(Vector3i,ArrayF<Vector3i comma 4>);
    Inst_Helper(Vector4i,Array<Vector3i>);Inst_Helper(Vector4i,ArrayF<Vector3i comma 4>);
#undef Inst_Helper
#define Inst_Helper(T_VEC,T_ARRAY) \
template int Element_Faces<T_ARRAY >(const T_VEC& v,const int ds,T_ARRAY& faces)
    Inst_Helper(Vector4i,Array<Vector3i>);Inst_Helper(Vector4i,ArrayF<Vector3i comma 4>);
#undef Inst_Helper
#undef comma

	//////////////////////////////////////////////////////////////////////////
	////Mesh edges
	template<int d,int e_d> void Get_Edges(const SimplicialMesh<d,e_d>& mesh,Array<Vector2i>& edges)
	{
		Hashset<Vector2i> edge_hashset;ArrayF<Vector2i,6> element_edges;
		for(const auto& vtx:mesh.elements){
			int n=MeshFunc::Element_Edges(vtx,element_edges);
			for(int i=0;i<n;i++)edge_hashset.insert(Sorted(element_edges[i]));}
		for(const auto& edge:edge_hashset)edges.push_back(edge);
	}	

	template<int d,int e_d> void Get_Vertices(const SimplicialMesh<d,e_d>& mesh,Array<int>& vertices)
	{
		Hashset<int> vtx_hashset;
		for(const auto& vtx:mesh.elements){
			for(int i=0;i<e_d;i++)vtx_hashset.insert(vtx[i]);}
		for(const auto& v:vtx_hashset)vertices.push_back(v);
	}

#define comma ,
#define Inst_Helper(d,e_d) \
template void Get_Edges(const SimplicialMesh<d,e_d>&,Array<Vector2i>&); \
template void Get_Vertices(const SimplicialMesh<d,e_d>&,Array<int>&);
Inst_Helper(2,2);Inst_Helper(2,3);Inst_Helper(3,2);Inst_Helper(3,3);Inst_Helper(3,4);
#undef Inst_Helper
#undef comma

	template<int d> real Average_Edge_Length(const Array<Vector<real,d> >& X,const Array<Vector<int,3> >& tri)
	{
		real length=(real)0;int n=0;
		for(auto i=0;i<tri.size();i++){const Vector3i& v=tri[i];
			length+=((X[v[0]]-X[v[1]]).norm()+(X[v[1]]-X[v[2]]).norm()+(X[v[2]]-X[v[0]]).norm());n+=3;}		
		if(n>0)return length/(real)n;return length;
	}
	template real Average_Edge_Length<2>(const Array<Vector<real,2> >&,const Array<Vector3i>&);
	template real Average_Edge_Length<3>(const Array<Vector<real,3> >&,const Array<Vector3i>&);
	
	template<int d> real Average_Edge_Length(const Array<Vector<real,d> >& X,const Array<Vector<int,2> >& seg)
	{
		real length=(real)0;int n=(int)seg.size();
		for(auto i=0;i<n;i++){const Vector2i& v=seg[i];
			length+=((X[v[0]]-X[v[1]]).norm());}		
		if(n>0)return length/(real)n;return length;
	}
	template real Average_Edge_Length<2>(const Array<Vector<real,2> >&,const Array<Vector2i>&);
	template real Average_Edge_Length<3>(const Array<Vector<real,3> >&,const Array<Vector2i>&);
	
    //////////////////////////////////////////////////////////////////////////
    ////Mesh converter
    void Volumetric_Mesh_Surface(const Array<Vector4i>& vol_elements,Array<Vector3i>& surf_elements)
    {
        ArrayF<Vector3i,4> faces;Hashtable<Vector3i,int> boundary_face_hashtable;
        for(int i=0;i<(int)vol_elements.size();i++){const Vector4i& vtx=vol_elements[i];
            Element_Faces(vtx,faces);for(int i=0;i<4;i++){const Vector3i& key=Sorted(faces[i]);
                auto find=boundary_face_hashtable.find(key);
                if(find!=boundary_face_hashtable.end()){boundary_face_hashtable.erase(key);}
                else{int sign=Is_Same_Order(key,faces[i])?1:-1;boundary_face_hashtable.insert(std::make_pair(key,sign));}}}
        for(const auto& iter:boundary_face_hashtable){Vector3i tri=iter.first;
            if(iter.second<0)tri=Reversed(iter.first);surf_elements.push_back(tri);}
    }

	void Volumetric_Mesh_Surface(const Array<Vector3i>& vol_elements,Array<Vector2i>& surf_elements)
	{
		ArrayF<Vector2i,3> edges;Hashtable<Vector2i,int> boundary_edge_hashtable;
		for(int i=0;i<(int)vol_elements.size();i++){const Vector3i& vtx=vol_elements[i];Element_Edges(vtx,edges);
			for(int i=0;i<3;i++){const Vector2i& key=Sorted(edges[i]);
				auto find=boundary_edge_hashtable.find(key); 
				if(find!=boundary_edge_hashtable.end()){boundary_edge_hashtable.erase(key);}
				else{int sign=(key==edges[i])?1:-1;boundary_edge_hashtable.insert(std::make_pair(key,sign));}}}
		for(const auto& iter:boundary_edge_hashtable){Vector2i edge=iter.first;
			if(iter.second<0)edge=Reversed(edge);surf_elements.push_back(edge);}
	}

	Vector2 Barycentric_Coord(const Vector2& p0,const Vector2& p1,const Vector2& p2,const Vector2& p)
	{
		Matrix2 m;m<<p0[0]-p[0],p1[0]-p2[0],p0[1]-p2[1],p1[1]-p2[1];return m.inverse()*(p-p2);
	}

	Vector3 Barycentric_Coord(const Vector3& p0,const Vector3& p1,const Vector3& p2,const Vector3& p3,const Vector3& p)
	{
		////TOIMPL
		return Vector3::Zero();
	}

	template<> Vector2 Barycentric_Coord<2>(const ArrayF<Vector2,3>& vtx,const Vector2& p)
	{
		return Barycentric_Coord(vtx[0],vtx[1],vtx[2],p);
	}

	template<> Vector3 Barycentric_Coord<3>(const ArrayF<Vector3,4>& vtx,const Vector3& p)
	{
		return Barycentric_Coord(vtx[0],vtx[1],vtx[2],vtx[3],p);
	}

	template<class T,int d> T Barycentric_Interpolation(const ArrayF<T,d+1>& values,const Vector<real,d>& coord)
	{
		T v=Zero<T>();real c=(real)1;
		for(int i=0;i<d;i++){c-=coord[i];v+=values[i]*coord[i];}v+=values[d]*c;return v;
	}
	#define Inst_Helper(T,d)	\
	template T Barycentric_Interpolation<T,d>(const ArrayF<T,d+1>&,const Vector<real,d>&);
	Inst_Helper(real,2);Inst_Helper(real,3);Inst_Helper(Vector2,2);Inst_Helper(Vector3,3);
	#undef Inst_Helper
	
	template<int d> bool Inside(const ArrayF<Vector<real,d>,d+1>& vtx,const Vector<real,d>& p)
	{
		Vector<real,d> bc=Barycentric_Coord<d>(vtx,p);real c=1;for(int i=0;i<d;i++){if(bc[i]<0||bc[i]>1)return false;c-=bc[i];}if(c<0||c>1)return false;return true;
	}
	#define Inst_Helper(d)	\
	template bool Inside<d>(const ArrayF<Vector<real,d>,d+1>&,const Vector<real,d>&);
	Inst_Helper(2);Inst_Helper(3);
	#undef Inst_Helper
	
	//////////////////////////////////////////////////////////////////////////
	////Mesh transformation
	template<int d> Box<d> Bounding_Box(const Array<Vector<real,d> >& vertices)
	{
		Box<d> box=Box<d>::Infi_Min();
		for(auto& v:vertices){box.min_corner=box.min_corner.cwiseMin(v);box.max_corner=box.max_corner.cwiseMax(v);}return box;
	}
	template Box<2> Bounding_Box<2>(const Array<Vector2>&);
	template Box<3> Bounding_Box<3>(const Array<Vector3>&);

	template<int d> Box<d> Bounding_Box(const Vector<real,d>* vertices,int vn)
	{
		Box<d> box=Box<d>::Infi_Min();
		for(int i=0;i<vn;i++){box.min_corner=box.min_corner.cwiseMin(vertices[i]);box.max_corner=box.max_corner.cwiseMax(vertices[i]);}return box;
	}
	template Box<2> Bounding_Box<2>(const Vector2*,int);
	template Box<3> Bounding_Box<3>(const Vector3*,int);
		
	template<int d> void Rescale(Array<Vector<real,d> >& vertices,const real longest_length)
	{
		Box<d> box=Bounding_Box<d>(vertices);Vector<real,d> length=box.Edge_Lengths();int axis=AuxFunc::Max_Index(length);
		real rescale=(length[axis]>(real)0)?longest_length/length[axis]:(real)1;for(auto& v:vertices)v*=rescale;
	}
	template void Rescale<2>(Array<Vector2>&,const real);
	template void Rescale<3>(Array<Vector3>&,const real);

	template<int d> Vector<real,d> Center(const Array<Vector<real,d> >& vertices)
	{
		Vector<real,d> sum=Vector<real,d>::Zero();for(auto& v:vertices)sum+=v;
		return sum/=(real)vertices.size();
	}
	
	template Vector2 Center<2>(const Array<Vector2>&);
	template Vector3 Center<3>(const Array<Vector3>&);

	template<int d> void Translate(Array<Vector<real,d> >& vertices,const Vector<real,d>& trans)
	{for(auto& v:vertices)v+=trans;}
	template void Translate<2>(Array<Vector2>&,const Vector2&);
	template void Translate<3>(Array<Vector3>&,const Vector3&);

	template<int d> void Translate_Center_To(Array<Vector<real,d> >& vertices,const Vector<real,d>& target)
	{Vector<real,d> center=Center<d>(vertices);Vector<real,d> trans=target-center;Translate<d>(vertices,trans);}
	template void Translate_Center_To<2>(Array<Vector2>&,const Vector2&);
	template void Translate_Center_To<3>(Array<Vector3>&,const Vector3&);

	void Rotate(Array<Vector2>& vertices,const real a)
	{for(auto& v:vertices){Vector2 r;r[0]=cos(a)*v[0]-sin(a)*v[1];r[1]=sin(a)*v[0]+cos(a)*v[1];v=r;}}
	void Rotate(Array<Vector3>& vertices,const AngleAxis& rot)
	{for(auto& v:vertices)v=rot*v;}

	void Rotate(Array<Vector2>& vertices, const Vector2& rot) {
		const real a = -atan2(rot[0], rot[1]);
		for (auto& v : vertices) { Vector2 r; r[0] = cos(a)*v[0] - sin(a)*v[1]; r[1] = sin(a)*v[0] + cos(a)*v[1]; v = r; }}

	void Rotate(Array<Vector3>& vertices, const Vector3& rot) {
		Quaternion quat = Quaternion::FromTwoVectors(rot.normalized(), Vector3::Unit(1));
		{for (auto& v : vertices)v = quat * v; }}

	void RotateBack(Array<Vector2>& vertices, const Vector2& rot) {
		const real a = atan2(rot[0], rot[1]);
		for (auto& v : vertices) { Vector2 r; r[0] = cos(a)*v[0] - sin(a)*v[1]; r[1] = sin(a)*v[0] + cos(a)*v[1]; v = r; }}
	void RotateBack(Array<Vector3>& vertices, const Vector3& rot) {
		Quaternion quat = Quaternion::FromTwoVectors(Vector3::Unit(1), rot.normalized());
		{for (auto& v : vertices)v = quat * v; }}

	//////////////////////////////////////////////////////////////////////////
	////Mesh initialization
	
	////Five tetrahedrons per cube, PhysBAM implementation TETRAHEDRON_MESH.cpp
	template<class T_MESH> void Initialize_Lattice_Mesh_Implementation(const Vector3i& counts,const real dx,T_MESH* mesh,const Array<int>* flag)
	{
		const int m=counts[0]+1,n=counts[1]+1,p=counts[2]+1;Grid<3> grid(counts,dx);
		for(int i=1;i<=m-1;i++)for(int j=1;j<=n-1;j++)for(int k=1;k<=p-1;k++){
			if(flag!=nullptr&&(*flag)[grid.Cell_Index(Vector3i(i-1,j-1,k-1))]==-1)continue;
			if((i+j+k)%2==0){
				mesh->elements.push_back(Vector4i(i+m*(j-1)+m*n*(k-1),i+1+m*(j-1)+m*n*(k-1),i+m*j+m*n*(k-1),i+m*(j-1)+m*n*k));
				mesh->elements.push_back(Vector4i(i+1+m*(j-1)+m*n*(k-1),i+1+m*(j-1)+m*n*k,i+1+m*j+m*n*k,i+m*(j-1)+m*n*k));
				mesh->elements.push_back(Vector4i(i+m*j+m*n*(k-1),i+1+m*j+m*n*(k-1),i+1+m*j+m*n*k,i+1+m*(j-1)+m*n*(k-1)));
				mesh->elements.push_back(Vector4i(i+m*j+m*n*k,i+1+m*j+m*n*k,i+m*(j-1)+m*n*k,i+m*j+m*n*(k-1)));
				mesh->elements.push_back(Vector4i(i+1+m*(j-1)+m*n*(k-1),i+m*(j-1)+m*n*k,i+1+m*j+m*n*k,i+m*j+m*n*(k-1)));}
			else{
				mesh->elements.push_back(Vector4i(i+m*(j-1)+m*n*(k-1),i+1+m*(j-1)+m*n*(k-1),i+1+m*j+m*n*(k-1),i+1+m*(j-1)+m*n*k));
				mesh->elements.push_back(Vector4i(i+m*(j-1)+m*n*(k-1),i+m*j+m*n*(k-1),i+m*j+m*n*k,i+1+m*j+m*n*(k-1)));
				mesh->elements.push_back(Vector4i(i+m*j+m*n*k,i+1+m*(j-1)+m*n*k,i+m*(j-1)+m*n*k,i+m*(j-1)+m*n*(k-1)));
				mesh->elements.push_back(Vector4i(i+m*j+m*n*k,i+1+m*j+m*n*k,i+1+m*(j-1)+m*n*k,i+1+m*j+m*n*(k-1)));
				mesh->elements.push_back(Vector4i(i+m*j+m*n*k,i+m*(j-1)+m*n*(k-1),i+1+m*j+m*n*(k-1),i+1+m*(j-1)+m*n*k));}}
		for(auto& e:mesh->elements){e-=Vector4i::Ones();}
		for(int k=0;k<p;k++)for(int j=0;j<n;j++)for(int i=0;i<m;i++)mesh->Vertices().push_back(Vector3((real)i,(real)j,(real)k)*dx);
		if(flag)Prune_Unused_Vertices(mesh);
	}

 	void Initialize_Lattice_Mesh(const Vector3i& counts,const real dx,TetrahedronMesh<3>* mesh,const Array<int>* flag)
	{
		Initialize_Lattice_Mesh_Implementation<TetrahedronMesh<3> >(counts,dx,mesh,flag);
	}

	void Initialize_Lattice_Mesh(const Vector2i& counts,const real dx,TriangleMesh<2>* mesh,const Array<int>* flag/*=nullptr*/)
	{
		Initialize_Herring_Bone_Mesh(counts[0]+1,counts[1]+1,dx,mesh);	////TODO: implement flag
	}
		
	template<int d> void Initialize_Lattice_Mesh(const Vector<int,d>& counts,const real dx,HexMesh<d>* mesh,const Array<int>* flag/*=nullptr*/)
	{
		Grid<d> grid(counts,dx);
		iterate_node(iter,grid){const Vector<int,d>& node=iter.Coord();const Vector<real,d> pos=grid.Node(node);mesh->Vertices().push_back(pos);}
		iterate_cell(iter,grid){const Vector<int,d>& cell=iter.Coord();const int idx=grid.Cell_Index(cell);
			if(flag!=nullptr&&(*flag)[idx]==-1)continue;
			ArrayF<int,Pow(2,d)> v;grid.Cell_Incident_Node_Indices(cell,v);
			Vector<int,Pow(2,d)> vtx;for(int i=0;i<Pow(2,d);i++)vtx[i]=v[i];
			mesh->Elements().push_back(vtx);}
		if(flag)Prune_Unused_Vertices(mesh);
	}

	template void Initialize_Lattice_Mesh<2>(const Vector2i&,const real,HexMesh<2>*,const Array<int>*);
	template void Initialize_Lattice_Mesh<3>(const Vector3i&,const real,HexMesh<3>*,const Array<int>*);

	////Cone mesh, n+1 vertices, n triangles
	void Initialize_Cone_Mesh(const real r,const real h,const int n,TriangleMesh<3>* mesh,const int axis)	
	{
		int a0,a1;if(axis==0){a0=1;a1=2;}else if(axis==1){a0=2;a1=0;}else{a0=0;a1=1;}
		mesh->Vertices().push_back(Vector3::Unit(axis)*h);
		real delta=(real)two_pi/(real)n;
		for(int i=0;i<n;i++){
			Vector3 v=Vector3::Unit(a0)*r*cos(delta*(real)i)+Vector3::Unit(a1)*r*sin(delta*(real)i);
			mesh->Vertices().push_back(v);}
		for(int i=1;i<n;i++)mesh->elements.push_back(Vector3i(0,i,i+1));
		mesh->elements.push_back(Vector3i(0,n,1));
	}

	////Cylinder mesh, 2n vertices, 2n triangles
	void Initialize_Cylinder_Mesh(const real r,const real h,const int n,TriangleMesh<3>* mesh,const int axis)	
	{
		int a0,a1;if(axis==0){a0=1;a1=2;}else if(axis==1){a0=2;a1=0;}else{a0=0;a1=1;}
		real delta=(real)two_pi/(real)n;
		for(int i=0;i<n;i++){
			Vector3 p0=Vector3::Unit(a0)*(r*cos(delta*(real)i))+Vector3::Unit(a1)*(r*sin(delta*(real)i));
			Vector3 p1=p0+Vector3::Unit(axis)*h;
			mesh->Vertices().push_back(p0);
			mesh->Vertices().push_back(p1);}
		for(int j=0;j<n-1;j++){int i=2*j;
			mesh->elements.push_back(Vector3i(i,i+2,i+1));
			mesh->elements.push_back(Vector3i(i+1,i+2,i+3));}
		mesh->elements.push_back(Vector3i(2*n-2,0,2*n-1));
		mesh->elements.push_back(Vector3i(2*n-1,0,1));
	}

	void Initialize_Circle_Mesh(const real r,const int n,SegmentMesh<2>* mesh)
	{
		real angle=two_pi/(real)n;
		for(int i=0;i<n;i++){
			real theta=angle*(real)i;
			Vector2 pos(r*cos(theta),r*sin(theta));
			mesh->Vertices().push_back(pos);
			mesh->Elements().push_back(Vector2i(i,(i+1)%n));}
	}

	void Initialize_Oval_Mesh(const real a,const real b,const int n,SegmentMesh<2>* mesh)
	{
		real angle=two_pi/(real)n;
		for(int i=0;i<n;i++){
			real theta=angle*((real)i+(real).5);
			Vector2 pos(a*cos(theta),b*sin(theta));
			mesh->Vertices().push_back(pos);
			mesh->Elements().push_back(Vector2i(i,(i+1)%n));}
	}

	template<int d> real Initialize_Segment_Mesh(const Vector<real, d>& v0, const Vector<real, d>& v1, int seg_num, SegmentMesh<d>* mesh, bool include_v0, bool include_v1) {
		int N = seg_num - 1 + (int)include_v0 + (int)include_v1;
		real h = (v1 - v0).norm() / seg_num;
		Vector<real, d> u = (v1 - v0).normalized();
		if (include_v0) mesh->Vertices().push_back(v0);
		for (int i = 0; i < seg_num; i++) {
			mesh->Vertices().push_back(v0 + i * h * u);
		}
		if (include_v1) mesh->Vertices().push_back(v1);
		for (int i = 0; i < mesh->Vertices().size() - 1; i++) {
			mesh->Elements().push_back(Vector2i(i, i + 1));
		}
		return h;
	}
	template real Initialize_Segment_Mesh<2>(const Vector<real, 2>& v0, const Vector<real, 2>& v1, int seg_num, SegmentMesh<2>* mesh, bool include_v0, bool include_v1);
	template real Initialize_Segment_Mesh<3>(const Vector<real, 3>& v0, const Vector<real, 3>& v1, int seg_num, SegmentMesh<3>* mesh, bool include_v0, bool include_v1);

	template<class T_MESH> void Initialize_Herring_Bone_Mesh(const int m, const int n, const real dx, T_MESH* mesh, int axis_0, int axis_1)
	{
		using VectorD=Vector<real,mesh->Dim()>;
		mesh->elements.resize(2*(m-1)*(n-1));int t=0;
		for(int i=1;i<=m-1;i++)for(int j=1;j<=n-1;j++){ // counterclockwise node ordering
			if(i%2){mesh->elements[t++]=Vector3i(i+m*(j-1),i+1+m*(j-1),i+m*j);mesh->elements[t++]=Vector3i(i+1+m*(j-1),i+1+m*j,i+m*j);}
			else{mesh->elements[t++]=Vector3i(i+m*(j-1),i+1+m*(j-1),i+1+m*j);mesh->elements[t++]=Vector3i(i+m*(j-1),i+1+m*j,i+m*j);}}
		for(size_type i=0;i<mesh->elements.size();i++){mesh->elements[i]-=Vector3i::Ones();}
			//int tmp=mesh->elements[i][1];mesh->elements[i][1]=mesh->elements[i][2];mesh->elements[i][2]=tmp;}
		for(int j=0;j<n;j++)for(int i=0;i<m;i++){VectorD pos=VectorD::Zero();pos[axis_0]=(real)i*dx;pos[axis_1]=(real)j*dx;mesh->Vertices().push_back(pos);}
	}

	#define Inst_Helper(T1) \
	template void Initialize_Herring_Bone_Mesh<T1 >(const int,const int,const real,T1*,int,int);
	Inst_Helper(TriangleMesh<2>);
	Inst_Helper(TriangleMesh<3>);
	#undef Inst_Helper

	template<int d> void Subdivide(const TriangleMesh<d>* input,TriangleMesh<d>* output)
	{
		Array<Vector2i> edges;Get_Edges(*input,edges);
		output->Vertices()=input->Vertices();
		Hashtable<Vector2i,int> edge_vtx_hashtable;
		for(const auto& e:edges){
			Vector<real,d> pos=(real).5*(input->Vertices()[e[0]]+input->Vertices()[e[1]]);
			output->Vertices().push_back(pos);
			int i=(int)output->Vertices().size()-1;
			edge_vtx_hashtable.insert(std::make_pair(e,i));}

		for(const auto& v:input->elements){int v3,v4,v5;
			{auto search=edge_vtx_hashtable.find(Sorted(Vector2i(v[0],v[1])));if(search==edge_vtx_hashtable.end())continue;v3=search->second;}
			{auto search=edge_vtx_hashtable.find(Sorted(Vector2i(v[1],v[2])));if(search==edge_vtx_hashtable.end())continue;v4=search->second;}
			{auto search=edge_vtx_hashtable.find(Sorted(Vector2i(v[2],v[0])));if(search==edge_vtx_hashtable.end())continue;v5=search->second;}
			output->elements.push_back(Vector3i(v[0],v3,v5));
			output->elements.push_back(Vector3i(v3,v[1],v4));
			output->elements.push_back(Vector3i(v5,v4,v[2]));
			output->elements.push_back(Vector3i(v3,v4,v5));}
	}

	template<int d> void Subdivide(TriangleMesh<d>* mesh)
	{
		Array<Vector2i> edges;Get_Edges(*mesh,edges);
		Hashtable<Vector2i,int> edge_vtx_hashtable;
		for(const auto& e:edges){
			Vector<real,d> pos=(real).5*(mesh->Vertices()[e[0]]+mesh->Vertices()[e[1]]);
			mesh->Vertices().push_back(pos);
			int i=(int)mesh->Vertices().size()-1;
			edge_vtx_hashtable.insert(std::make_pair(e,i));}

		auto n=mesh->elements.size();
		for(auto i=0;i<n;i++){const Vector3i v=mesh->elements[i];int v3,v4,v5;
			{auto search=edge_vtx_hashtable.find(Sorted(Vector2i(v[0],v[1])));if(search==edge_vtx_hashtable.end())continue;v3=search->second;}
			{auto search=edge_vtx_hashtable.find(Sorted(Vector2i(v[1],v[2])));if(search==edge_vtx_hashtable.end())continue;v4=search->second;}
			{auto search=edge_vtx_hashtable.find(Sorted(Vector2i(v[2],v[0])));if(search==edge_vtx_hashtable.end())continue;v5=search->second;}
			mesh->elements.push_back(Vector3i(v[0],v3,v5));
			mesh->elements.push_back(Vector3i(v3,v[1],v4));
			mesh->elements.push_back(Vector3i(v5,v4,v[2]));
			mesh->elements[i]=Vector3i(v3,v4,v5);}
	}

	#define Inst_Helper(d) \
	template void Subdivide<d>(const TriangleMesh<d>*,TriangleMesh<d>*);\
	template void Subdivide<d>(TriangleMesh<d>*);
	Inst_Helper(2);Inst_Helper(3);
	#undef Inst_Helper

	void Initialize_Icosahedron_Mesh(const real scale,/*rst*/TriangleMesh<3>* mesh)
	{
		////http://donhavey.com/blog/tutorials/tutorial-3-the-icosahedron-sphere/
		const real tao=1.61803399f;
		real vtx_pos[12][3]={{1,tao,0},{-1,tao,0},{1,-tao,0},{-1,-tao,0},{0,1,tao},{0,-1,tao},{0,1,-tao},{0,-1,-tao},{tao,0,1},{-tao,0,1},{tao,0,-1},{-tao,0,-1}};
		int ele[20][3]={{0,1,4},{1,9,4},{4,9,5},{5,9,3},{2,3,7},{3,2,5},{7,10,2},{0,8,10},{0,4,8},{8,2,10},{8,4,5},{8,5,2},{1,0,6},{11,1,6},{3,9,11},{6,10,7},{3,11,7},{11,6,7},{6,0,10},{9,1,11}};		

		mesh->Clear();
		int vtx_num=12;mesh->Vertices().resize(vtx_num);for(int i=0;i<vtx_num;i++){mesh->Vertices()[i]=Vector3(vtx_pos[i][0],vtx_pos[i][1],vtx_pos[i][2])*scale;}
		int ele_num=20;mesh->elements.resize(ele_num);for(int i=0;i<ele_num;i++)mesh->elements[i]=Vector3i(ele[i][0],ele[i][1],ele[i][2]);
	}

	void Initialize_Sphere_Mesh(const real r,/*rst*/TriangleMesh<3>* mesh,const int sub/*=2*/)
	{
		Initialize_Icosahedron_Mesh(r,mesh);for(int i=0;i<sub;i++)Subdivide(mesh);
		for(auto& v:mesh->Vertices()){real length=v.norm();real rs=r/length;v*=rs;}
	}

	void Initialize_Ellipsoid_Mesh(const real r,/*rst*/TriangleMesh<3>* mesh,const real a, const real b, const real c, const int sub/*=2*/)
	{
		Initialize_Icosahedron_Mesh(r, mesh); for (int i = 0; i < sub; i++)Subdivide(mesh);
		Matrix3 trans = Eigen::Scaling(a, b, c);
		for (auto& v : mesh->Vertices()) { real length = v.norm(); real rs = r / length; v *= rs; v = trans*v; }
	}

	//////////////////////////////////////////////////////////////////////////
	////Mesh operations
	template<class T_MESH> void Merge(const Array<T_MESH*>& meshes, T_MESH* merged_mesh)
	{
		merged_mesh->vertices->clear();merged_mesh->elements.clear();
		int n=0;for(auto& m:meshes){
			for(auto& p:m->Vertices()){
				merged_mesh->vertices->push_back(p);}
			for(auto& e:m->elements){
				auto e1=e;for(int i=0;i<m->Element_Dim();i++)e1[i]+=n;
				merged_mesh->elements.push_back(e1);}
			n+=(int)m->vertices->size();}
	}

	#define Inst_Helper(T_MESH) \
	template void Merge<T_MESH>(const Array<T_MESH*>&, T_MESH*);
	Inst_Helper(TriangleMesh<2>);Inst_Helper(TriangleMesh<3>);Inst_Helper(TetrahedronMesh<3>);
	#undef Inst_Helper

	template<class T_MESH> void Prune_Unused_Vertices(T_MESH* mesh)
	{
		const int d=mesh->Dim();const int e_d=mesh->Element_Dim();
		Hashset<int> flag;for(int i=0;i<mesh->Elements().size();i++)for(int j=0;j<e_d;j++){flag.insert(mesh->Elements()[i][j]);}
		Array<int> mapping(mesh->Vertices().size(),-1);for(int i=0,c=0;i<mesh->Vertices().size();i++){if(flag.find(i)!=flag.end())mapping[i]=c++;}
		Array<Vector<real,d> > pruned_vertices;for(int i=0;i<mesh->Vertices().size();i++)if(mapping[i]!=-1)pruned_vertices.push_back(mesh->Vertices()[i]);
		mesh->Vertices()=pruned_vertices;
		for(int i=0;i<mesh->Elements().size();i++){for(int j=0;j<e_d;j++){int k=mesh->Elements()[i][j];mesh->Elements()[i][j]=mapping[k];}}
	}

	#define Inst_Helper(T_MESH) \
	template void Prune_Unused_Vertices<T_MESH>(T_MESH*);
	Inst_Helper(TriangleMesh<2>);Inst_Helper(TriangleMesh<3>);Inst_Helper(TetrahedronMesh<3>);Inst_Helper(HexMesh<2>);Inst_Helper(HexMesh<3>);
	#undef Inst_Helper
};
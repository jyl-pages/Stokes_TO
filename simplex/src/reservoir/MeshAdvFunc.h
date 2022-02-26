//////////////////////////////////////////////////////////////////////////
// Mesh functions
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __MeshAdvFunc_h__
#define __MeshAdvFunc_h__
#include "Common.h"

template<int> class Box;
template<int,int> class SimplicialMesh;
template<int> class TriangleMesh;
template<int> class TetrahedronMesh;
template<int> class SegmentMesh;

namespace MeshFunc{
	////Triangulation
	void Initialize_Circle_Mesh(const Vector2& c,const real r,TriangleMesh<2>* mesh,const int sub_n);
	void Initialize_Circle_Mesh(const Vector3& c,const real r,const Vector3& normal,TriangleMesh<3>* mesh,const int sub_n);
	void Initialize_Ring_Mesh(const Vector2& c,const real R,const real r,TriangleMesh<2>* mesh,const int sub_n);
	void Initialize_Polygon_Mesh(const Array<Vector2>& vtx,TriangleMesh<2>* mesh,const real dx);
	void Initialize_Polygon_Mesh(const Array<Array<Vector2> >& vtx,TriangleMesh<2>* mesh,const real dx);
	void Initialize_Rectangle_Mesh(const Vector2& domain_min,const Vector2& domain_size,TriangleMesh<2>* mesh,const real dx);
	
	////////Ellipsoid mesh: TOIMPL
	//void Initialize_Ellipsoid_Mesh(const real a,const real b,const real c,TriangleMesh<3>* mesh){}

	////Tetrahedronization
	void Initialize_Sphere_Mesh(const Vector3& c,const real r,TetrahedronMesh<3>& volume_mesh,const real dx);
	void Initialize_Cuboid_Mesh(const Vector3& domain_min, const Vector3& domain_size, TetrahedronMesh<3>& volume_mesh, const real dx);
	void Initialize_Cylinder_Mesh(const Vector3& c,const real r,const real h,TetrahedronMesh<3>& volume_mesh,const real dx);
	void Initialize_Polygon3d_Mesh(const Array<Array<Vector3>>& vtx,TetrahedronMesh<3>& volume_mesh, const real dx);

	void Initialize_Polygon_Mesh_From_File(const char* file_name, TriangleMesh<2>* mesh);
	void Initialize_Polygon_Mesh_From_File(const char* file_name, TetrahedronMesh<3>* mesh);

	template<typename T> bool Read_Vtx_From_File(const char* file_name, Array<Array<T>>& vtx, int d, real& dx, int mode = 0); // mode 0 : read contour
};
#endif
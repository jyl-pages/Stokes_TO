//////////////////////////////////////////////////////////////////////////
// Mesh
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Mesh_h__
#define __Mesh_h__
#include <fstream>
#include "Common.h"
#include "TypeFunc.h"

////Simplicial mesh
template<int d,int e_d> class SimplicialMesh
{Typedef_VectorDii(d);Typedef_VectorEi(e_d);
public:
	ArrayPtr<VectorD> vertices;
	Array<VectorEi> elements;
	ArrayPtr<VectorD> normals;		////default=nullptr
	ArrayPtr<Vector2> textures;		////default=nullptr

	////constructors
	SimplicialMesh(const ArrayPtr<VectorD> _vertices=nullptr);
	SimplicialMesh(const SimplicialMesh<d,e_d>& copy){*this=copy;}
	SimplicialMesh<d,e_d>& operator=(const SimplicialMesh<d,e_d>& copy);

	////attribute access
	static constexpr int Dim() {return d;}
	static constexpr int Element_Dim() {return e_d;}
	virtual Array<VectorD>& Vertices(){return *vertices.get();}
	virtual const Array<VectorD>& Vertices() const {return *vertices.get();}
	virtual Array<VectorEi>& Elements(){return elements;}
	virtual const Array<VectorEi>& Elements() const {return elements;}

	virtual Array<VectorD>* Normals(){return normals.get();}
	virtual const Array<VectorD>* Normals() const {return normals.get();}
	virtual Array<Vector2>* Textures(){return textures.get();}
	virtual const Array<Vector2>* Textures() const {return textures.get();}

	virtual void Clear(){if(vertices)vertices->clear();elements.clear();if(normals)normals->clear();if(textures)textures->clear();}

	////IO
	virtual void Write_Binary(std::ostream& output) const;
	virtual void Read_Binary(std::istream& input);
	virtual void Write_To_File_3d(const std::string& file_name) const;
	virtual void Write_Text(std::ostream& output) const;
	virtual void Read_Text(std::istream& input);
};

template<int d> class TetrahedronMesh : public SimplicialMesh<d,4>
{Typedef_VectorDii(d);using VectorEi=Vector4i;using Base=SimplicialMesh<d,4>;
public:
	using Base::vertices;using Base::elements;
	TetrahedronMesh(const ArrayPtr<VectorD> _vertices=nullptr);
};

template<int d> class TriangleMesh : public SimplicialMesh<d,3>
{Typedef_VectorDii(d);using VectorEi=Vector3i;using Base=SimplicialMesh<d,3>;
public:
	using Base::vertices;using Base::elements;using Base::normals;
	TriangleMesh(const ArrayPtr<VectorD> _vertices=nullptr);
};

template<int d> class SegmentMesh : public SimplicialMesh<d,2>
{Typedef_VectorDii(d);using VectorEi=Vector2i;using Base=SimplicialMesh<d,2>;
public:
	using Base::vertices;using Base::elements;
	SegmentMesh(const ArrayPtr<VectorD> _vertices=nullptr);
};

template<int d> class QuadMesh : public SimplicialMesh<d,4>
{Typedef_VectorDii(d);Typedef_VectorEi(4);using Base=SimplicialMesh<d,4>;
public:
	using Base::vertices;using Base::elements;
	QuadMesh(const ArrayPtr<VectorD> _vertices=nullptr);
};

////HexMesh<2> is a QuadMesh
template<int d> class HexMesh : public SimplicialMesh<d,Pow(2,d)>
{Typedef_VectorDii(d);Typedef_VectorEi(Pow(2,d));using Base=SimplicialMesh<d,Pow(2,d)>;
public:
	using Base::vertices;using Base::elements;
	HexMesh(const ArrayPtr<VectorD> _vertices=nullptr);
};

template<int d> class SurfaceQuadMesh : public SimplicialMesh<d,Pow(2,d-1)>		////segment mesh in 2D, quad mesh in 3D
{Typedef_VectorDii(d);Typedef_VectorEi(Pow(2,d-1));using Base=SimplicialMesh<d,Pow(2,d-1)>;
public:
	using Base::vertices;using Base::elements;
	SurfaceQuadMesh(const ArrayPtr<VectorD> _vertices=nullptr);
};

////CODESAMPLE: this is an example of template alias instantiation, should use 2, 3 instead of d
template<int d> using VolumetricMesh=typename If<d==2,TriangleMesh<2>,TetrahedronMesh<3> >::Type;
template<int d> using SurfaceMesh=typename If<d==2,SegmentMesh<2>,TriangleMesh<3> >::Type;

////Aux functions
////Dimension and type conversion
template<class MESH_T1,class MESH_T2> void Dim_Conversion(const MESH_T1& mesh2,/*rst*/MESH_T2& mesh3);

#endif


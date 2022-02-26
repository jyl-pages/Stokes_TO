//////////////////////////////////////////////////////////////////////////
// Linear Hex FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyLinearFemHex_h__
#define __SoftBodyLinearFemHex_h__
#include "Hashtable.h"
#include "LinearFemFunc.h"
#include "BoundaryCondition.h"
#include "Field.h"
#include "TypeFunc.h"
#include "Mesh.h"

template<int d> class SoftBodyLinearFemHex
{Typedef_VectorDii(d);Typedef_MatrixD(d);Typedef_VectorEi(Pow(2,d));
public:
	static constexpr int ei=Pow(2,d);
	HexMesh<d>* mesh=nullptr;
	SparseMatrixT K;
	VectorX u;
	VectorX f;

	Array<ElasticParam> materials;
	Array<int> material_id;
	BoundaryConditionMesh<d> bc;

	void Initialize(HexMesh<d>& _mesh);
	void Allocate_K();
	void Update_K_And_f();
	void Solve();

	void Add_Material(real youngs,real poisson);

	void Set_Fixed(const int node);
	void Set_Displacement(const int node,const VectorD& dis);
	void Set_Force(const int node,const VectorD& force);
	void Add_Force(const int node,const VectorD& force);
};

#endif
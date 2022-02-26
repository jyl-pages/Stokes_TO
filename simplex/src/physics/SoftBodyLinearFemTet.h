//////////////////////////////////////////////////////////////////////////
// Linear Tetrahedral FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyLinearFemTet_h__
#define __SoftBodyLinearFemTet_h__
#include "Hashtable.h"
#include "LinearFemFunc.h"
#include "BoundaryCondition.h"
#include "Mesh.h"
#include "Field.h"

template<int d> class SoftBodyLinearFemTet
{Typedef_VectorDii(d);Typedef_MatrixD(d);Typedef_VectorEi(d+1);
public:
	VolumetricMesh<d>* mesh=nullptr;
	SparseMatrixT K;
	VectorX u;
	VectorX f;

	Array<ElasticParam> materials;
	Array<int> material_id;
	BoundaryConditionMesh<d> bc;

	bool use_dynamics=false;
	Array<VectorD> X0;		////material space X
	DiagonalMatrixT M;		////mass matrix
	VectorX v;				////velocity
	VectorX a;				////acceleration
	VectorX bf;
	bool use_damping=false;
	real damping=1e-3;

	void Initialize(VolumetricMesh<d>& _mesh);
	void Allocate_K();
	void Update_K_And_f();
	void Solve();
	void Add_Material(real youngs,real poisson);
	void Set_Fixed(const int node);
	void Set_Displacement(const int node,const VectorD& dis);
	void Set_Force(const int node,const VectorD& force);
	void Add_Force(const int node,const VectorD& force);

	////dynamic
	void Advance(const real dt,const real time){/*TOIMPL*/}

	////Helper functions
	int Node_Num() const {return (int)material_id.size();}
	void Compute_Strain(Array<MatrixD>& strains) const;
    void Compute_Stress(Array<MatrixD>& stresses) const;
	void Compute_Von_Mises_Stress(Array<real>& von_mises,const Array<MatrixD>* stresses=nullptr) const;
protected:
	void Compute_Tet_Tensor_Helper(const VectorX& u,Array<MatrixD>& tensors,Array<MatrixX>* E=nullptr) const;
};

#endif
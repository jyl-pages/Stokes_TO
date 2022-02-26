//////////////////////////////////////////////////////////////////////////
// Nonlinear Tetrahedral FEM
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SoftBodyNonlinearFemTet_h__
#define __SoftBodyNonlinearFemTet_h__
#include "Hashtable.h"
#include "LinearFemFunc.h"
#include "BoundaryCondition.h"
#include "Mesh.h"
#include "Field.h"
#include "Particles.h"

template<int d> class SoftBodyNonlinearFemTet
{Typedef_VectorDii(d);Typedef_MatrixD(d);Typedef_VectorEi(d+1);
public:
	Particles<d> particles;
	std::shared_ptr<VolumetricMesh<d> > mesh=nullptr;

	real density=(real)1;
	DiagonalMatrix<real> M;
	SparseMatrixT K;
	VectorX u;
	VectorX f;

	Array<ElasticParam> materials;
	Array<int> material_id;
	BoundaryConditionMesh<d> bc;
	Array<MatrixD> Dm_inv;
	Array<VectorD> X0;
	Array<real> vol;

	bool use_explicit=true;
	bool use_body_force=true;
	VectorD g=VectorD::Unit(1)*-1;
	bool use_damping=false;
	real damping=(real)10;	
	std::function<void(const real dt)> Collision;
	std::function<void(const real dt,const real time)> Kinematic_Boundary_Condition;

	////implicit time integration
	VectorX v;				////velocity
	VectorX a;				////acceleration
	VectorX bf;

	void Initialize(VolumetricMesh<d>& _mesh);
	virtual void Advance(const real dt,const real time);
	void Advance_Explicit(const real dt,const real time);
	void Advance_Implicit(const real dt,const real time){/*TOIMPL*/}

	void Add_Material(real youngs,real poisson);
	void Set_Fixed(const int node);
	void Set_Displacement(const int node,const VectorD& dis);
	void Set_Force(const int node,const VectorD& force);
	void Add_Force(const int node,const VectorD& force);
	void Clear_Force();
	void Set_Rest_Shape(const Array<VectorD>& _X0);
	
	////Helper functions
	void Strain_To_Stress(const MatrixD& strain,const MatrixX& E,MatrixD& stress);
	void Area_Weighted_Normals(const Array<VectorD>& X,const int e,MatrixD& normals);
	inline real Mass(const int i) const {return M.diagonal()[i*d];}
	inline Array<VectorD>& X(){return particles.XRef();}
	inline Array<VectorD>& V(){return particles.VRef();}
	inline Array<VectorD>& F(){return particles.FRef();}
	inline auto& E(){return mesh->Elements();}
	inline const Array<VectorD>& X() const {return particles.XRef();}
	inline const Array<VectorD>& V() const {return particles.VRef();}
	inline const Array<VectorD>& F() const {return particles.FRef();}
	inline const auto& E() const {return mesh->Elements();}

	inline int Vtx_Num() const {return particles.Size();}
	inline int Ele_Num() const {return (int)mesh->Elements().size();}

	////omp acc
	bool use_omp=false;
	void Advance_Explicit_Omp(const real dt,const real time);
};

#endif
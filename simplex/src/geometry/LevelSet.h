//////////////////////////////////////////////////////////////////////////
// Level set
// Copyright (c) (2018-), Bo Zhu, Xingyu Ni
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __LevelSet_h__
#define __LevelSet_h__

#include "Grid.h"
#include "Field.h"
#include "FaceField.h"
#include "Interpolation.h"
#include "GeometryPrimitives.h"

template<int d> class LevelSet
{Typedef_VectorDii(d);
public:
	Grid<d> grid;
	Field<real,d> phi;
protected: 
	std::unique_ptr<Interpolation<d> > intp=nullptr;
public:
	LevelSet(){}
	LevelSet(const Grid<d>& _grid);
	void Initialize(const Grid<d>& _grid);
	void Set_By_Geom(ImplicitGeometry<d>& geom);

	real Phi(const VectorD& pos) const;
	VectorD Normal(const VectorD& pos) const;	////TOFIX: fix finite difference on the boundary
	VectorD Gradient(const VectorD& pos) const;	////TOFIX: fix finite difference on the boundary

	VectorD Closest_Point(const VectorD& pos,const real epsilon=(real)0) const;
	VectorD Closest_Point_With_Iterations(const VectorD& pos,const int max_iter=5) const;
	
	void Update_Normals(FaceField<real,d>& normals) const;
	void Update_Normals(Field<VectorD,d>& normals) const;
	real Curvature(const VectorD& pos) const;
	
	real Cell_Fraction(const VectorDi& cell) const;		////approximate cell volume using phi value
	real Total_Volume() const;

	////Helper functions
	static real Sign(const real phi){return phi<=(real)0?(real)-1:(real)1;}
	static bool Interface(const real phi_1,const real phi_2){return Sign(phi_1)!=Sign(phi_2);}
	static real Theta(const real phi_1,const real phi_2){return phi_1/(phi_1-phi_2);}

	//////////////////////////////////////////////////////////////////////////
	////Fast marching
	void Fast_Marching(const real band_width=(real)-1);
protected:
	////Fast marching auxiliary data structures
	Field<real,d> tent;
	Array<char> done;
	Array<int> intf_cells;
	using PRI=std::pair<real,int>;
	std::priority_queue<PRI,Array<PRI>,std::greater<PRI> > heap;
	////Fast marching helper functions
	void Initialize_Fast_Marching();
	void Precondition(const real band_width);
	real Solve_Eikonal(const VectorDi &cell);
	bool Solve_Quadratic(const real p1,const real p2,const real dx,real& rst);
	bool Solve_Quadratic(const real p1,const real p2,const real p3,const real dx,real& rst);
};

#endif
//////////////////////////////////////////////////////////////////////////
// Marching cube
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __MarchingCubes__
#define __MarchingCubes__
#include "Grid.h"
#include "Mesh.h"
template<int> class LevelSet;

template<int d> class MarchingCubes{};

template<> class MarchingCubes<2>
{const int d=2;Typedef_VectorDii(2);
public:
    LevelSet<2>& levelset;
	Grid<2>& grid;
    std::shared_ptr<SurfaceMesh<2> > mesh;
    real contour_value;
public:
    MarchingCubes(LevelSet<2>& _levelset,std::shared_ptr<SurfaceMesh<2> > _mesh=nullptr,const real _contour_value=(real)0);
	void Marching();
protected:
	////Helper data structures and functions
    unsigned int Square_Type(real v0,real v1,real v2,real v3,real iso_value=(real)0);
	bool Edge_Has_Vertex(real phi1,real phi2,real iso_value=(real)0);
};

template<> class MarchingCubes<3>
{const int d=3;Typedef_VectorDii(3);
public:
    LevelSet<3>& levelset;
	Grid<3>& grid;
    std::shared_ptr<SurfaceMesh<3> > mesh;
    real contour_value;
public:
    MarchingCubes(LevelSet<3>& _levelset,std::shared_ptr<SurfaceMesh<3> > _mesh=nullptr,const real _contour_value=(real)0);
	void Marching();
protected:
	////Helper data structures and functions
	static const int edge_table[256];
	static const int triangle_table[256][16];
	Grid<3> edge_grid[3];
    Field<int,3> edge_grid_array[3];
	bool Has_Edge(const int edge_type,const int edge_index);
	bool Has_Particle_On_Edge(const VectorDi& cell_index,const int edge_index);
	void Get_Edge_Vertex_Index(const VectorDi& cell_index,int edge_index,/*result*/VectorDi& v1,/*result*/VectorDi& v2);
	int Get_Particle_Index_On_Edge(const VectorDi& cell_index,const int edge_index);
	void Set_Particle_Index_On_Edge(const VectorDi& cell_index,const int edge_index,const int value_input);
	unsigned int Get_Cell_Type(real v0,real v1,real v2,real v3,real v4,real v5,real v6,real v7,real iso_value);
};
#endif

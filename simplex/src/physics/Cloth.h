//////////////////////////////////////////////////////////////////////////
// Soft body mass spring
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Cloth_h__
#define __Cloth_h__

#include "Hashtable.h"
#include "Grid.h"
#include "SoftBodyMassSpring.h"

class ClothMassSpring: public SoftBodyMassSpring<3>
{public:
	Typedef_VectorDii(3);Typedef_MatrixD(3);typedef SoftBodyMassSpring<3> Base;using Base::springs;

	Grid<2> grid;

	void Initialize_Springs()
	{
		Hashset<Vector2i> edge_hashset;
		iterate_cell_d(iter,grid,2){Vector2i cell=iter.Coord();
			ArrayF<Vector2i,6> cell_edges;Get_Cell_Edges(cell,grid,cell_edges);
			for(int i=0;i<6;i++)edge_hashset.insert(Unique_Ordered(cell_edges[i]));}
		springs.clear();
		for(auto& e:edge_hashset)springs.push_back(e);
	}

	void Get_Cell_Edges(const Vector2i& cell,const Grid<2>& grid,ArrayF<Vector2i,6>& edges)
	{
		int p0=grid.Node_Index(cell);
		int p1=grid.Node_Index(cell+Vector2i(0,1));
		int p2=grid.Node_Index(cell+Vector2i(1,0));
		int p3=grid.Node_Index(cell+Vector2i(1,1));

		edges[0]=Vector2i(p0,p1);
		edges[1]=Vector2i(p2,p3);
		edges[2]=Vector2i(p0,p2);
		edges[3]=Vector2i(p1,p3);
		edges[4]=Vector2i(p0,p3);
		edges[5]=Vector2i(p1,p2);
	}
};
#endif


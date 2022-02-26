//////////////////////////////////////////////////////////////////////////
// Color Multigrid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __ColorGrid_h__
#define __ColorGrid_h__

#include "Grid.h"

template<int d> class ColorGrid
{	Typedef_VectorDii(d);
public:
	static int Number_Of_Colors(){return (int)pow(2,d);}

	////return an index from 0-(2^d-1)
	static int Color(const Vector<int,d>& node)
	{
		static Vector<int,d> counts=Vector<int,d>::Ones()*2;
		Vector<int,d> node0;for(int i=0;i<d;i++)node0[i]=node[i]%2;
		return Grid<d>::Index(node0,counts);
	}
	
	////this function is used for generating regular domain colors
	static void Color(const Vector<int,d>& counts,Array<int>& colored_node_ptr,Array<int>& colored_node_indices)
	{	
		Grid<d> grid(counts,1);
		int color_n=Number_Of_Colors();

		Array<Array<int> > buckets(color_n);
		for(int i=0;i<color_n;i++)buckets[i].reserve(grid.Number_Of_Cells()/color_n+1);

		iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();int c=Color(cell);
			buckets[c].push_back(grid.Cell_Index(cell));}

		colored_node_ptr.resize(color_n+1);
		colored_node_indices.resize(grid.Number_Of_Cells());

		int i=0;for(int c=0;c<color_n;c++){
			colored_node_ptr[c]=i;
			std::copy(buckets[c].begin(),buckets[c].end(),&colored_node_indices[i]);
			i+=(int)buckets[c].size();}
		colored_node_ptr[color_n]=i;
	}

	////irregular domain colors: the way to map from mat_id cells to matrix index needs to be consistent with the multigrid solver
	static void Color(const Vector<int,d>& counts,const Array<short>& mat_id,Array<int>& colored_node_ptr,Array<int>& colored_node_indices)
	{	
		Grid<d> grid(counts,1);
		int color_n=Number_Of_Colors();

		Array<Array<int> > buckets(color_n);
		for(int i=0;i<color_n;i++)buckets[i].reserve(grid.Number_Of_Cells()/color_n+1);

		int cell_num=counts.prod();int n=0;
		for(int i=0;i<cell_num;i++){if(mat_id[i]==-1)continue;
			const VectorDi cell=grid.Cell_Coord(i);int c=Color(cell);
			buckets[c].push_back(n);n++;}

		colored_node_ptr.resize(color_n+1);
		colored_node_indices.resize(n);

		int i=0;for(int c=0;c<color_n;c++){
			colored_node_ptr[c]=i;
			std::copy(buckets[c].begin(),buckets[c].end(),&colored_node_indices[i]);
			i+=(int)buckets[c].size();}
		colored_node_ptr[color_n]=i;
	}

	////this function is not used
	static void Color(const Vector<int,d>& counts,const Array<int>& node_indices,Array<int>& colored_node_ptr,Array<int>& colored_node_indices)
	{
		int color_n=Number_Of_Colors();

		Array<Array<int> > buckets(color_n);
		for(int i=0;i<color_n;i++)buckets[i].reserve(((int)node_indices.size())/color_n+1);
		
		for(int i=0;i<(int)node_indices.size();i++){
			VectorDi node=Grid<d>::Coord(node_indices[i],counts);int c=Color(node);
			buckets[c].push_back(node_indices[i]);}
		
		colored_node_ptr.resize(color_n+1);
		colored_node_indices.resize((int)node_indices.size());

		int i=0;for(int c=0;c<color_n;c++){
			colored_node_ptr[c]=i;
			std::copy(buckets[c].begin(),buckets[c].end(),&colored_node_indices[i]);
			i+=(int)buckets[c].size();}
		colored_node_ptr[color_n]=i;
	}
};

#endif
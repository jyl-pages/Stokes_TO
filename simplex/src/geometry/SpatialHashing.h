//////////////////////////////////////////////////////////////////////////
// Spatial hashing grid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __SpatialHashing_h__
#define __SpatialHashing_h__
#include <limits>
#include <iostream>
#include "Common.h"
#include "Hashtable.h"
#include "Grid.h"
#include "GeometryPrimitives.h"
#include "MeshFunc.h"

template<int d,int bs=512> class SpatialHashing
{Typedef_VectorDii(d);using Bucket=ArrayF<int,bs>;
public:
	Grid<d> grid;
	Hashtable<VectorDi,Bucket> voxels;

	SpatialHashing(){grid=Grid<d>(VectorDi::Ones()*16,(real)1);}
	SpatialHashing(const Grid<d>& _grid){Initialize(_grid);}
	void Initialize(const real dx){grid=Grid<d>(VectorDi::Ones()*16,dx);}
	void Initialize(const Grid<d>& _grid){grid=_grid;}

	static constexpr int Bs(){return bs;}

	void Update_Voxels(const Array<VectorD>& points)
	{Clear_Voxels();for(auto i=0;i<points.size();i++)Add_Point(i,points[i]);}

	template<typename F> void Update_Voxels(const Array<VectorD>& points,F& valid)
	{Clear_Voxels();for(auto i=0;i<points.size();i++){if(!valid(i))continue;Add_Point(i,points[i]);}}

	void Clear_Voxels(){voxels.clear();}

	void Refresh_Voxels(){for(auto& iter:voxels)iter->second[0]=0;}	////Clear each bucket

	bool Add_Point(const int idx,const VectorD& pos)
	{
		VectorDi cell=grid.Cell_Coord(pos);
		return Add_Element_To_Voxel(idx,cell);
	}

	bool Add_Points(const int idx,const Array<VectorD>& pos)
	{
		Box<d> bounding_box=MeshFunc::Bounding_Box<d>(pos);
		VectorDi start_cell=grid.Cell_Coord(bounding_box.min_corner);
		VectorDi end_cell=grid.Cell_Coord(bounding_box.max_corner);
		iterate_range(iter,start_cell,end_cell){
			bool rst=Add_Element_To_Voxel(idx,iter.Coord());if(!rst)return false;}return true;
	}

	template<typename F> bool Add_Points(const int idx,const Box<d>& box,F& valid)
	{
		VectorDi start_cell=grid.Cell_Coord(box.min_corner);
		VectorDi end_cell=grid.Cell_Coord(box.max_corner);
		iterate_range(iter,start_cell,end_cell){
			if(!valid(grid.Center(iter.Coord())))continue;
			bool rst=Add_Element_To_Voxel(idx,iter.Coord());if(!rst)return false;}return true;
	}

	////return an array with fixed size (nbs), the number of elements in nbs is stored in its first element 
	////this function assumes all nb points lie within the 9/27 neighboring cells
	bool Find_Nbs(const VectorD& pos,const Array<VectorD>& points,const real R,ArrayF<int,bs*4>& nbs) const
	{
		real r_sqr=R*R;nbs[0]=0;
		VectorDi cell=grid.Cell_Coord(pos);
		for (int i = 0; i < grid.Number_Of_Nb_R(); i++) {
			const VectorDi& nb_cell = grid.Nb_R(cell, i);
			auto iter = voxels.find(nb_cell); if (iter == voxels.end())continue;
			const Bucket& bucket = iter->second;
			for (int j = 0; j < bucket[(size_t)0]; j++) {
				int point_idx = bucket[j + 1];
				const VectorD& point_pos = points[point_idx];
				if ((pos - point_pos).squaredNorm() < r_sqr) {
					if (nbs[0] >= bs * 4 - 1) { std::cerr << "Error: [SpatialHashingGrid] nbs full" << std::endl; return false; }
					nbs[nbs[0] + 1] = point_idx; nbs[0]++;}}}
		return true;
	}

	int Find_Nearest_Nb(const VectorD& pos,const Array<VectorD>& points)
	{
		int min_dis_idx=-1;real min_dis_sq=std::numeric_limits<real>::max();
		VectorDi cell=grid.Cell_Coord(pos);
		for(int i=0;i<grid.Number_Of_Nb_R();i++){
			const VectorDi& nb_cell=grid.Nb_R(cell,i);
			auto iter=voxels.find(nb_cell);if(iter==voxels.end())continue;
			Bucket& bucket=iter->second;
			for(int j=0;j<bucket[0];j++){int point_idx=bucket[j+1];
				const VectorD& point_pos=points[point_idx];
				real dis_sq=(pos-point_pos).squaredNorm();
				if(dis_sq<min_dis_sq){min_dis_sq=dis_sq;min_dis_idx=point_idx;}}}	
		return min_dis_idx;
	}

	void Collect_Elements_From_Cell(const VectorDi& cell,Hashset<int>& elements) const
	{
		auto iter=voxels.find(cell);if(iter==voxels.end())return;
		const Bucket& bucket=iter->second;
		for(int j=0;j<bucket[0];j++){elements.insert(bucket[j+1]);}
	}

	void Collect_Elements_From_Cells(const Array<VectorDi>& cells,Hashset<int>& elements) const
	{
		for(auto& cell:cells){Collect_Elements_From_Cell(cell,elements);}		
	}

	bool Add_Element_To_Voxel(const int idx,const VectorDi& cell)
	{
		auto iter=voxels.find(cell);
		if(iter==voxels.end())iter=voxels.insert(std::make_pair(cell,Bucket{0})).first;
		Bucket& bucket=iter->second;
		if (bucket[0] >= bs - 1) {
			std::cerr << "Error: [SpatialHashingGrid] bucket full: " << cell.transpose() << ", " << bucket[0] << std::endl; return false;
		}
		bucket[bucket[0]+1]=idx;bucket[0]++;return true;	
	}

	//////////////////////////////////////////////////////////////////////////
	////n-ring neighbors
	Array<VectorDi> n_ring_nbs;

	void Initialize_NR_Nbs(const real R)
	{
		int n=std::ceil(R/grid.dx);
		grid.Nb_NR(n,n_ring_nbs);		
	}

	bool Find_Nbs_NR(const VectorD& pos,const Array<VectorD>& points,const real R,ArrayF<int,bs*4>& nbs) const
	{
		real r_sqr=R*R;nbs[0]=0;
		VectorDi cell=grid.Cell_Coord(pos);
		int n_nbs=(int)(n_ring_nbs.size());

		for(int i=0;i<n_nbs;i++){
			const VectorDi nb_cell=cell+n_ring_nbs[i];
			auto iter=voxels.find(nb_cell);if(iter==voxels.end())continue;
			const Bucket& bucket=iter->second;
			for(int j=0;j<bucket[0];j++){int point_idx=bucket[j+1];
				const VectorD& point_pos=points[point_idx];
				if((pos-point_pos).squaredNorm()<r_sqr){
					if(nbs[0]>=bs*4-1){std::cerr<<"Error: [SpatialHashingGrid] n-ring nbs full"<<std::endl;return false;}
					nbs[nbs[0]+1]=point_idx;nbs[0]++;}}}
		return true;
	}
};
#endif

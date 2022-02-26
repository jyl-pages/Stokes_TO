#ifndef __LevelSetFunc_h__
#define __LevelSetFunc_h__
#include "LevelSet.h"
#include "PointSet.h"
#include "NeighborSearcher.h"
#include "SpatialHashing.h"

namespace LevelSetFunc
{
	// Linux compilation error: non-local lambda expression cannot have a capture-default
	//template<int d> void Points_To_Level_Set(Array<Vector<real,d> >& x,/*result*/LevelSet<d>& levelset,
	//	std::function<bool(const int)> valid = [=](const int)->bool{return true;})
	//{
	//	////TODO: implement Bridson's paper		
	//}

	template<int d> void Point_Set_To_Level_Set(PointSet<d>& point_set,/*result*/LevelSet<d>& levelset,real pointset_nb_width,real leveset_nb_width,const bool reinit=true)
	{	Typedef_VectorDii(d);
		auto& grid=levelset.grid;
		iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
			VectorD pos=grid.Center(cell);
			if(abs(levelset.phi(cell))>pointset_nb_width)continue;
			VectorD proj_pos=point_set.Project_To_Surface(pos);
			VectorD normal=point_set.Normal(proj_pos);
			VectorD vec=proj_pos-pos;
			real sign=normal.dot(vec)>=(real)0?(real)1:(real)0;
			levelset.phi(cell)=vec.norm()*sign;}
		if(reinit)levelset.Fast_Marching(/*leveset_nb_width*/);
	}


	//Linux compilation error: error: ¡®point_set¡¯ was not declared in this scope; did you mean ¡®PointSet¡¯?
	//template<int d> void Points_To_Level_Set(Array<Vector<real,d> >& x,/*result*/LevelSet<d>& levelset,real pointset_nb_width,real leveset_nb_width,const bool reinit=true)
	//{	Typedef_VectorDii(d);
	//	auto& grid=levelset.grid;
	//	iterate_cell(iter,grid){const VectorDi& cell=iter.Coord();
	//		VectorD pos=grid.Center(cell);
	//		if(abs(levelset.phi(cell))>pointset_nb_width)continue;
	//		VectorD proj_pos=point_set.Project_To_Surface(pos);
	//		VectorD normal=point_set.Normal(proj_pos);
	//		VectorD vec=proj_pos-pos;
	//		real sign=normal.dot(vec)>=(real)0?(real)1:(real)0;
	//		levelset.phi(cell)=vec.norm()*sign;}
	//	if(reinit)levelset.Fast_Marching(leveset_nb_width);
	//}
};

#endif
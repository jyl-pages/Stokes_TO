//////////////////////////////////////////////////////////////////////////
// Solve the viscosity equation
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __Vorticity_h__
#define __Vorticity_h__
#include "MacGrid.h"
#include "FaceField.h"
#include "Field.h"
#include "PoissonOld.h"
#include "FluidFunc.h"
#include "BoundaryCondition.h"

namespace Vorticity
{
	void Vorticity_To_Velocity(MacGrid<2>& mac_grid,const Field<real,2>& vorticity_on_node,FaceField<real,2>& velocity,
			const Hashtable<Vector2i,real>* psi_N_values=nullptr,const Hashtable<int,real>* psi_D_values=nullptr)
	{static const int d=2;Typedef_VectorDii(d);
		PoissonOld<d> poisson;
		Grid<d> node_grid(mac_grid.grid.node_counts,mac_grid.grid.dx,mac_grid.grid.domain_min-VectorD::Ones()*mac_grid.grid.dx*(real).5);
		poisson.Initialize(node_grid);
		poisson.use_multigrid_solver=true;
		iterate_cell(iter,node_grid){const VectorDi& cell=iter.Coord();
			poisson.rhs(cell)=-vorticity_on_node(cell);}

		if(psi_N_values!=nullptr)poisson.bc.psi_N_values=*psi_N_values;
		//if(psi_D_values!=nullptr)poisson.bc.psi_D_values=*psi_D_values;
		for (const auto& p : *psi_D_values) {
			poisson.bc.Set_Psi_D(mac_grid.grid.Cell_Coord(p.first), (ushort)CellType::Fluid, p.second);
		}

		poisson.Build_And_Solve();
		FluidFunc::Curl_On_Face(mac_grid,poisson.p,velocity);
	}
};

template<int d> class Vorticity_To_Velocity{};

template<> class Vorticity_To_Velocity<2>
{static const int d=2;Typedef_VectorDii(d);
public:
	bool initialized=false;
	MacGrid<d> mac_grid;
	Grid<d> node_grid;
	PoissonOld<d> poisson;

	void Initialize(MacGrid<2>& _mac_grid,const Field<real,2>& vorticity_on_node,FaceField<real,2>& velocity,
			const Hashtable<Vector2i,real>* psi_N_values=nullptr,const Hashtable<int,real>* psi_D_values=nullptr)
	{
		mac_grid=_mac_grid;
		node_grid=mac_grid.grid.Cell_Grid_To_Node_Grid();

		poisson.Initialize(node_grid);
		poisson.use_multigrid_solver=true;
		iterate_cell(iter,node_grid){const VectorDi& cell=iter.Coord();
			poisson.rhs(cell)=-vorticity_on_node(cell);}

		if(psi_N_values!=nullptr)poisson.bc.psi_N_values=*psi_N_values;
		//if(psi_D_values!=nullptr)poisson.bc.psi_D_values=*psi_D_values;	
		for (const auto& p : *psi_D_values) {
			poisson.bc.Set_Psi_D(mac_grid.grid.Cell_Coord(p.first), (ushort)CellType::Fluid, p.second);
		}
		poisson.Build();
		initialized=true;
	}

	void Solve(const Field<real,2>& vorticity_on_node,FaceField<real,2>& velocity)
	{
		iterate_cell(iter,node_grid){const VectorDi& cell=iter.Coord();
			poisson.rhs(cell)=-vorticity_on_node(cell);}
		poisson.Update_Rhs_With_BC();
		poisson.Solve();
		FluidFunc::Curl_On_Face(mac_grid,poisson.p,velocity);
	}
};

#endif

//////////////////////////////////////////////////////////////////////////
// Incompressible fluid solver on an Eulerian grid, supporting free interface and surface tension
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidEulerFreeSurface_h__
#define __FluidEulerFreeSurface_h__
#include "FluidEuler.h"
#include "LevelSet.h"
#include "ProjectionIrregular.h"
#include "Timer.h"

template<int d> class FluidEulerFreeSurface : public FluidEuler<d>
{Typedef_VectorDii(d);using Base=FluidEuler<d>;
public:
	using Base::mac_grid;
	using Base::velocity;
	//using Base::alpha;		////not using alpha now
	using Base::type;			////supporting type: Fluid, NB, Air
	using Base::bc;
	using Base::g;
	using Base::use_body_force;	////default false
	using Base::Enforce_Boundary_Conditions;

	LevelSet<d> levelset;
	int narrow_band_cell_num = 5;
	real narrow_band_width;
	int dirac_band_cell_num = 3;
	real dirac_band_width;
	ProjectionIrregular<d> projection;		////use_jump_condition (for surface tension) and sigma were moved to the projection class
	
	bool verbose=false;
public:
	FluidEulerFreeSurface(const SolverType &_solver_mode = SolverType::MULTIGRID_AUTO): projection(mac_grid, velocity, levelset, bc, _solver_mode) {}

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		Base::Initialize(cell_counts,dx,domain_min);
		Initialize_Interface();

		////setup projection callback functions
		//projection.Is_Fluid_Cell=std::bind(&FluidEulerFreeSurface<d>::Is_Fluid_Cell,this,std::placeholders::_1);
		//projection.Is_Fluid_Cell_Index=std::bind(&FluidEulerFreeSurface<d>::Is_Fluid_Cell_Index,this,std::placeholders::_1);
		projection.Is_Fluid_Interior_Cell_Index=std::bind(&FluidEulerFreeSurface<d>::Is_Fluid_Interior_Cell_Index,this,std::placeholders::_1);
		projection.Is_Interface_Face_Index=std::bind(&FluidEulerFreeSurface<d>::Is_Interface_Face_Index, this, std::placeholders::_1);
		projection.Dirac = std::bind(&FluidEulerFreeSurface<d>::Dirac, this, std::placeholders::_1);
	}

	virtual void Initialize_Interface()
	{
		levelset.Initialize(mac_grid.grid);
		narrow_band_width=mac_grid.grid.dx*(real)narrow_band_cell_num;	
		dirac_band_width = mac_grid.grid.dx * (real)dirac_band_cell_num;
	}

	virtual void Advance(const real dt)
	{
		Timer<real> timer;					if(verbose)timer.Reset();
		Advection(dt);						if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: advection");
		Apply_Body_Forces(dt);				if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: body force");
		Enforce_Incompressibility(dt);		if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: projection");
		Extrapolation(velocity);			if(verbose)timer.Elapse_And_Output_And_Reset("\t\ttimestep: extrapolation");
	}

	virtual void Advection(const real dt)
	{
		Timer<real> timer;if(verbose)timer.Reset();

		////advect velocity
		MacGrid<d> ghost_grid=mac_grid;FaceField<real,d> ghost_velocity=velocity;
		std::function<bool(const int,const VectorDi&)> valid_face=[&](const int axis,const VectorDi& face)->bool{return Is_Fluid_Nb_Face(axis,face);};
		Advection::Semi_Lagrangian<d>(dt,ghost_grid,ghost_velocity,mac_grid,velocity,valid_face);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: vel adv");

		////advect interface
		Field<real,d> ghost_phi=levelset.phi;
		std::function<bool(const VectorDi&)> valid_cell=[&](const VectorDi& cell)->bool{return Is_Fluid_Nb_Cell(cell);};
		Advection::Semi_Lagrangian(dt,ghost_velocity,ghost_grid,ghost_phi,mac_grid,levelset.phi,false,valid_cell);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: ls adv");

		levelset.Fast_Marching(narrow_band_width);
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: ls fmm");

		Update_Cell_Types();
		if(verbose)timer.Elapse_And_Output_And_Reset("\tadv: celltype");
	}

	virtual void Update_Cell_Types()
	{
		int cell_num=mac_grid.grid.Number_Of_Cells();
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			if(levelset.phi.array[i]<(real)0)type.array[i]=(ushort)CellType::Fluid;
			else if(levelset.phi.array[i]<narrow_band_width)type.array[i]=(ushort)CellType::NB;
			else type.array[i]=(ushort)CellType::Air;}

		////update boundary cell types
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			if(type.array[i]!=(ushort)CellType::Fluid)continue;
			const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			for(int j=0;j<mac_grid.grid.Number_Of_Nb_C();j++){
				const VectorDi& nb=mac_grid.grid.Nb_C(cell,j);
				if(!mac_grid.grid.Valid_Cell(nb)||(type(nb)==(ushort)CellType::NB||type(nb)==(ushort)CellType::Air)){
					type.array[i]=(ushort)CellType::BD;}}}
	}

	virtual void Apply_Body_Forces(const real dt)
	{
		if(!use_body_force)return;
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				VectorD pos=mac_grid.Face_Center(axis,face);real phi=levelset.Phi(pos);
				if(phi<narrow_band_width){velocity.face_fields[axis](face)+=g[axis]*dt;}}}
	}

	virtual void Enforce_Incompressibility(const real dt)
	{
		Enforce_Boundary_Conditions();
		projection.current_dt=dt;
		projection.Project();
	}

	////extrapolation with levelset
	virtual void Extrapolation(FaceField<real,d>& field_q)
	{
		Interpolation<d> intp(mac_grid);FaceField<real,d> ghost_field_q=field_q;

		////Fill ghost fields with fluid velocities
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				if(Is_Fluid_Face(axis,face))continue;

				int n=0;real fluid_v=(real)0;
				for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
					VectorDi nb_face=Grid<d>::Nb_C(face,i);
					if(!mac_grid.face_grids[axis].Valid_Node(nb_face))continue;
					if(Is_Fluid_Face(axis,nb_face)){fluid_v+=ghost_field_q(axis,nb_face);n++;}}
				if(n>0)ghost_field_q(axis,face)=fluid_v/(real)n;}}

		
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				if(bc.Is_Psi_N(axis,face))continue;

				VectorD pos=mac_grid.Face_Center(axis,face);real phi=levelset.Phi(pos);
				if(phi>(real)0){
					if(phi<narrow_band_width){
						VectorD interface_pos=levelset.Closest_Point(pos);
						field_q(axis,face)=intp.Interpolate_Faces(ghost_field_q,interface_pos,axis);}
					else field_q(axis,face)=(real)0;}}}
	}

protected:
	//////////////////////////////////////////////////////////////////////////
	////helper functions

	inline bool Is_Fluid_Interior_Cell(const VectorDi& cell) const 
	{return mac_grid.grid.Valid_Cell(cell)&&(type(cell)==(ushort)CellType::Fluid);}

	inline bool Is_Fluid_Interior_Cell_Index(const int cell_idx) const	////assuming index is always valid
	{return (type.array[cell_idx]==(ushort)CellType::Fluid);}

	inline bool Is_Fluid_Cell(const VectorDi& cell) const 
	{return mac_grid.grid.Valid_Cell(cell)&&(type(cell)==(ushort)CellType::Fluid||type(cell)==(ushort)CellType::BD);}

	inline bool Is_Fluid_Cell_Index(const int cell_idx) const	////assuming index is always valid
	{return (type.array[cell_idx]==(ushort)CellType::Fluid||type.array[cell_idx]==(ushort)CellType::BD);}

	inline bool Is_Fluid_Nb_Cell(const VectorDi& cell) const
	{return mac_grid.grid.Valid_Cell(cell)&&((type(cell)==(ushort)CellType::Fluid||type(cell)==(ushort)CellType::BD)||type(cell)==(ushort)CellType::NB);}

	inline bool Is_Fluid_Nb_Cell_Index(const int cell_idx) const 
	{return (type.array[cell_idx]==(ushort)CellType::Fluid||type.array[cell_idx]==(ushort)CellType::BD)||type.array[cell_idx]==(ushort)CellType::NB;}

	////a face is fluid if it is incident to a fluid cell
	inline bool Is_Fluid_Face(const int axis,const VectorDi& face) const
	{{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,0);if(Is_Fluid_Cell(cell))return true;}
	{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,1);if(Is_Fluid_Cell(cell))return true;}return false;}

	inline bool Is_Fluid_Nb_Face(const int axis,const VectorDi& face) const
	{{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,0);if(Is_Fluid_Nb_Cell(cell))return true;}
	{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,1);if(Is_Fluid_Nb_Cell(cell))return true;}return false;}

	inline bool Is_Interface_Face(const int axis, const VectorDi& face) const
	{
		VectorD pos = mac_grid.Face_Center(axis, face);
		if (bc.Is_Psi_N(axis, face)) return false;
		real phi = levelset.Phi(pos);
		return (phi > -narrow_band_width && phi < narrow_band_width);
	}

	inline bool Is_Interface_Face_Index(const std::pair<int, int> face_idx) const
	{
		int axis = face_idx.first;
		VectorDi face = mac_grid.face_grids[axis].Node_Coord(face_idx.second);
		VectorD pos = mac_grid.Face_Center(axis, face);
		if (bc.Is_Psi_N(axis, face)) return false;
		real phi = levelset.Phi(pos);
		return (phi > -narrow_band_width && phi < narrow_band_width);
	}

	inline real Dirac(const real phi) const
	{
		if (phi < -dirac_band_width) return 0;
		else if (phi > dirac_band_width) return 0;
		else return 0.5 * (1.0 + cos(pi * phi / dirac_band_width)) / dirac_band_width;
	}
};
#endif
//////////////////////////////////////////////////////////////////////////
// PIC-FLIP fluid
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __FluidPicFlip_h__
#define __FluidPicFlip_h__
#include <set>
#include "FluidEuler.h"
#include "Particles.h"
#include "RandomNumber.h"

template<int d> class FluidPicFlip : public FluidEuler<d>
{Typedef_VectorDii(d);using Base=FluidEuler<d>;
public:
	using Base::mac_grid;
	using Base::velocity;
	using Base::g;
	using Base::type;
	using Base::Initialize_Fields;
	using Base::Enforce_Incompressibility;
	using Base::Enforce_Boundary_Conditions;
	
	TrackerPoints<d> points;			////stores only position and velocity
	FaceField<real,d> face_mass;		////mass on faces for interpolation weights
	
	RandomNumber random_number;
	int cell_num=0;
	int npc=(int)pow(2,d);				////number of points per cell
	int npc_tol=1;						////number of points per cell tolerance

	virtual void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts,dx,domain_min);
		Initialize_Fields();
		face_mass.Resize(mac_grid.grid.cell_counts,(real)0);
		cell_num=mac_grid.grid.cell_counts.prod();
	}

	virtual void Advance(const real dt)
	{
		Timer<real> timer;
		timer.Reset();
		Advection(dt);
		timer.Elapse_And_Output_And_Reset("Advection");
		Update_Cell_Types();
		Transfer_Velocities_From_Particle_To_Grid();
		Velocity_Extrapolation();
		timer.Elapse_And_Output_And_Reset("P2G");
		Apply_Body_Forces(dt);
		Enforce_Incompressibility();
		timer.Elapse_And_Output_And_Reset("Projection");
		Reseed_Particles();
		timer.Elapse_And_Output_And_Reset("Reseed");

		Transfer_Velocities_From_Grid_To_Particle();
		timer.Elapse_And_Output_And_Reset("G2P");
	}

	virtual void Advection(const real dt)
	{
		int p_size=points.Size();
		#pragma omp parallel for
		for(int i=0;i<p_size;i++){
			if(!Valid_Particle(i))continue;
			points.X(i)+=points.V(i)*dt;}
	}

	virtual void Apply_Body_Forces(const real dt)
	{
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				if(face_mass(axis,face)==(real)0)continue;
				VectorD pos=mac_grid.Face_Center(axis,face);
				velocity.face_fields[axis](face)+=g[axis]*dt;}}
	}

	void Transfer_Velocities_From_Grid_To_Particle()
	{
		Interpolation<d> intp(mac_grid.grid);
		int p_size=points.Size();
		#pragma omp parallel for
		for(int i=0;i<p_size;i++){
			if(!Valid_Particle(i))continue;
			VectorD v=intp.Interpolate_Face_Vectors_Cubic(velocity,points.X(i));
			points.V(i)=v;}
	}

	void Transfer_Velocities_From_Particle_To_Grid()
	{
		velocity.Fill((real)0);face_mass.Fill((real)0);
		Interpolation<d> intp(mac_grid.grid);

		int p_size=points.Size();
		//#pragma omp parallel for
		for(int i=0;i<p_size;i++){
			if(!Valid_Particle(i))continue;
			intp.Interpolate_Point_To_Faces_Cubic(points.X(i),points.V(i),velocity,face_mass);}

		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
			if(face_mass(axis,face)!=(real)0)velocity(axis,face)/=face_mass(axis,face);}}
	}

	virtual void Update_Cell_Types()
	{
		Field<short,d> pc(mac_grid.grid.cell_counts,0);
		type.Fill((ushort)CellType::Air);
		int p_size=points.Size();

		#pragma omp parallel for
		for(int i=0;i<p_size;i++){
			if(!Valid_Particle(i))continue;
			const VectorD& pos=points.X(i);
			VectorDi cell=mac_grid.grid.Cell_Coord(pos);
			if(!mac_grid.grid.Valid_Cell(cell))continue;	////particle outside grid
			type(cell)=(ushort)CellType::Fluid;}
		
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			const VectorDi cell=mac_grid.grid.Cell_Coord(i);
			if(Is_Air_Cell(cell)){
				bool is_bubble_cell=true;
				for(int j=0;j<mac_grid.grid.Number_Of_Nb_C();j++){
					const VectorDi nb=mac_grid.grid.Nb_C(cell,j);
					if(mac_grid.grid.Valid_Cell(nb)&&Is_Air_Cell(nb)){is_bubble_cell=false;break;}}
				if(is_bubble_cell){
					type(cell)=(ushort)CellType::Fluid;}}}
	}

	virtual void Reseed_Particles()
	{
		Field<short,d> pc(mac_grid.grid.cell_counts,0);

		////three cases for removing:
		////1. outside grid; 2. inside solid; 3. #points exceeds max;
		int p_size=points.Size();

		#pragma omp parallel for
		for(int i=0;i<p_size;i++){
			if(!Valid_Particle(i))continue;
			const VectorD& pos=points.X(i);
			VectorDi cell=mac_grid.grid.Cell_Coord(pos);
			bool outside=!mac_grid.grid.Valid_Cell(cell);
			
			bool inside_solid=!outside&&Is_Solid_Cell(cell);
			bool exceed_max_npc=!outside&&(pc(cell)>npc+npc_tol);

			if(outside||inside_solid||exceed_max_npc){
				#pragma omp critical
				{Remove_Particle(i);}}
			else pc(cell)++;}
		
		////two cases for reseeding: 
		////1. fluid cell with insufficient points; 2. interior air cell
		#pragma omp parallel for
		for(int i=0;i<cell_num;i++){
			const VectorDi cell=mac_grid.grid.Cell_Coord(i);
			int reseed_n=0;
			
			if(Is_Fluid_Cell(cell)&&pc(cell)<npc-npc_tol){
				bool is_interior_cell=true;
				for(int j=0;j<mac_grid.grid.Number_Of_Nb_C();j++){
					const VectorDi nb=mac_grid.grid.Nb_C(cell,j);
					if(mac_grid.grid.Valid_Cell(nb)&&!Is_Fluid_Cell(nb)){is_interior_cell=false;break;}}
				if(!is_interior_cell)continue;

				reseed_n=npc-pc(cell);
				Array<int> idx;
				#pragma omp critical
				{Calc_Seed_Particle_Index(npc-pc(cell),idx);}
				Seed_Particles(cell,idx);}}
	}

	inline void Calc_Seed_Particle_Index(const int n,Array<int>& idx)
	{idx.resize(n);for(int i=0;i<n;i++){idx[i]=Add_Particle();}}

	////parallel seeding
	void Seed_Particles(const VectorDi& cell,const int n)
	{
		Array<int> idx;
		Calc_Seed_Particle_Index(n,idx);
		Seed_Particles(cell,idx);
	}

	void Seed_Particles(const VectorDi& cell,const Array<int>& idx)
	{
		VectorD node_pos=mac_grid.grid.Node(cell);
		int n=(int)idx.size();
		for(int i=0;i<n;i++){int p=idx[i];
			VectorD pos=node_pos+random_number.VectorValue<d>()*mac_grid.grid.dx;
			points.X(p)=pos;
			points.V(p)=VectorD::Zero();}	////velocity will be interpolated from grid
	}

	void Velocity_Extrapolation()
	{
		FaceField<real,d> ghost_velocity=velocity;

		////Fill ghost fields with fluid velocities
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			#pragma omp parallel for
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				if(Is_Fluid_Face(axis,face)){continue;}
				int n=0;real fluid_v=(real)0;
				for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){
					VectorDi nb_face=Grid<d>::Nb_C(face,i);
					if(!mac_grid.face_grids[axis].Valid_Node(nb_face))continue;
					if(Is_Fluid_Face(axis,nb_face)){fluid_v+=velocity(axis,nb_face);n++;}}
				if(n>0){ghost_velocity(axis,face)=fluid_v/(real)n;}}}
		velocity=ghost_velocity;
	}

protected:
	//////////////////////////////////////////////////////////////////////////
	////particle system operations
	Heap<int,std::greater<int> > invalid_point_heap;			////for particle reseeding
	inline bool Valid_Particle(const int i) const {return points.I(i)!=-1;}
	inline int Append_Particle(){points.Add_Element();return points.Size()-1;}
	int Add_Particle(){if(!invalid_point_heap.empty()){int p=invalid_point_heap.top();invalid_point_heap.pop();points.I(p)=0;return p;}else return Append_Particle();}
	void Remove_Particle(const int i){invalid_point_heap.push(i);points.I(i)=-1;points.V(i)=VectorD::Zero();points.X(i)=VectorD::Zero();}

	inline bool Is_Fluid_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Fluid;}
	inline bool Is_Solid_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Solid;}
	inline bool Is_Air_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Air;}
	inline bool Is_Source_Cell(const VectorDi& cell) const {return type(cell)==(ushort)CellType::Source;}

	bool Is_Fluid_Face(const int axis,const VectorDi& face) const
	{
		{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,0);if(mac_grid.grid.Valid_Cell(cell)&&Is_Fluid_Cell(cell))return true;}
		{VectorDi cell=MacGrid<d>::Face_Incident_Cell(axis,face,1);if(mac_grid.grid.Valid_Cell(cell)&&Is_Fluid_Cell(cell))return true;}return false;
	}
};
#endif
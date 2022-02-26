#ifndef __BoundaryCondition_h__
#define __BoundaryCondition_h__
#include "Common.h"
#include "Hashtable.h"	
#include "MacGrid.h"
#include "FaceField.h"
#include "GeometryPrimitives.h"

enum class CellType: ushort {Fluid=0,Air,Solid,IB,Source,NB,BD};	////NB-narrowband,BD-boundary

template<int d> class BoundaryConditionGrid
{Typedef_VectorDii(d);
public:
	Hashtable<VectorDi,VectorD> psi_D_values;	////fixed points
	Hashtable<VectorDi,VectorD> forces;			////force

	BoundaryConditionGrid(){}
	void Set_Psi_D(const VectorDi& node,const VectorD& displacement)
	{psi_D_values[node]=displacement;}
	void Set_Force(const VectorDi& node,const VectorD& force)
	{forces[node]=force;}
	void Clear(){psi_D_values.clear();forces.clear();}
};

template<int d> class BoundaryConditionMacGrid
{Typedef_VectorDii(d);
public:
	MacGrid<d>& mac_grid;
	Hashtable<Vector2i,real> psi_N_values;					////[axis,face_index]->value
	Hashtable<int,std::pair<ushort,real> > psi_D_values;	////cell_index->[type,value]

	BoundaryConditionMacGrid(MacGrid<d>& _mac_grid):mac_grid(_mac_grid){}
	void Set_Psi_N(const int axis,const VectorDi& face,const real value=(real)0)
	{psi_N_values[Vector2i(axis,mac_grid.Face_Index(axis,face))]=value;}
	void Set_Psi_N(const int axis,const int face_index,const real value=(real)0)
	{psi_N_values[Vector2i(axis,face_index)]=value;}
	void Del_Psi_N(const int axis,const VectorDi& face)
	{psi_N_values.erase(Vector2i(axis,mac_grid.Face_Index(axis,face)));}
	void Del_Psi_N(const int axis,const int face_index)
	{psi_N_values.erase(Vector2i(axis,face_index));}
	void Set_Psi_D(const VectorDi& cell,const ushort type,const real value=(real)0)
	{psi_D_values[mac_grid.grid.Cell_Index(cell)]=std::make_pair(type,value);}
	void Del_Psi_D(const VectorDi& cell)
	{psi_D_values.erase(mac_grid.grid.Cell_Index(cell));}
	void Set_Psi_D(ImplicitGeometry<d>* geom, const ushort type = 0, const real value = 0);
	bool Is_Psi_N(const int axis,const VectorDi& face) const 
	{return psi_N_values.find(Vector2i(axis,mac_grid.Face_Index(axis,face)))!=psi_N_values.end();}
	bool Is_Psi_D(const VectorDi& cell) const
	{return psi_D_values.find(mac_grid.grid.Cell_Index(cell))!=psi_D_values.end();}
	real Psi_N_Value(const int axis, const VectorDi& face)const
	{return psi_N_values.at(Vector2i(axis, mac_grid.Face_Index(axis, face)));}
	ushort Psi_D_Type(const VectorDi& cell) 
	{return psi_D_values[mac_grid.grid.Cell_Index(cell)].first;}
	real Psi_D_Value(const VectorDi& cell)const
	{return psi_D_values.at(mac_grid.grid.Cell_Index(cell)).second;}
	void Enforce_Boundary_Conditions(FaceField<real, d>& v);
	void Enforce_Boundary_Conditions(Field<real, d>& v);
	void Write_Psi_D_To_File_3d(std::string file_name);
	void Write_Psi_N_To_File_3d(std::string file_name);

	////acceleration arrays: maintaining an array-copy of the hashtable for random access
	Array<Vector2i> psi_N_key_array;
	Array<real> psi_N_val_array;
	Array<int> psi_D_key_array;
	Array<std::pair<ushort,real> > psi_D_val_array;
	bool acc_arrays_init=false;

	void Build_Acceleration_Arrays()
	{
		int psi_N_n=(int)psi_N_values.size();
		psi_N_key_array.resize(psi_N_n);
		psi_N_val_array.resize(psi_N_n);
		int i=0;for(auto p:psi_N_values){psi_N_key_array[i]=p.first;psi_N_val_array[i]=p.second;i++;}

		int psi_D_n=(int)psi_D_values.size();
		psi_D_key_array.resize(psi_D_n);
		psi_D_val_array.resize(psi_D_n);
		i=0;for(auto p:psi_D_values){psi_D_key_array[i]=p.first;psi_D_val_array[i]=p.second;i++;}

		acc_arrays_init=true;
	}
};

template<int d> class BoundaryConditionMacGridViscosity
{Typedef_VectorDii(d);
public:
	MacGrid<d>& mac_grid;
	Hashtable<Vector2i,real> psi_D_values;	////[axis,face_index]->value, store psi_D velocity value on faces

	BoundaryConditionMacGridViscosity(MacGrid<d>& _mac_grid):mac_grid(_mac_grid){}
	void Set_Psi_D(const int axis,const VectorDi& face,const real value=(real)0)
	{psi_D_values[Vector2i(axis,mac_grid.Face_Index(axis,face))]=value;}
	bool Is_Psi_D(const int axis,const VectorDi& face) const 
	{return psi_D_values.find(Vector2i(axis,mac_grid.face_grids[axis].Node_Index(face)))!=psi_D_values.end();}

	void Enforce_Boundary_Conditions(FaceField<real,d>& v);
};

template<int d> class BoundaryConditionMesh
{Typedef_VectorDii(d);
public:
	Hashtable<int,VectorD> psi_D_values;	////fixed points
	Hashtable<int,VectorD> forces;			////force

	BoundaryConditionMesh(){}
	void Set_Psi_D(const int& node,const VectorD& displacement)
	{psi_D_values[node]=displacement;}
	void Set_Force(const int& node,const VectorD& force)
	{forces[node]=force;}	
	bool Is_Psi_D(const int& node) const
	{return psi_D_values.find(node)!=psi_D_values.end();}
};

#endif

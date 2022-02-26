//#####################################################################
// Face field
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __BitFaceField_h__
#define __BitFaceField_h__
#include "BitField.h"
#include "MacGrid.h"

template<int d> class BitFaceField
{Typedef_VectorDi(d);
public:
	MacGrid<d> mac_grid;
	ArrayF<BitField<d>,d> face_fields;

	BitFaceField(const VectorDi& cell_counts=VectorDi::Zero())
	{Resize(cell_counts);}
	
	BitFaceField(const VectorDi& cell_counts,const bool& value)
	{Resize(cell_counts,value);}
	
	BitFaceField(const BitFaceField<d>& copy)
	{*this=copy;}
	
	BitFaceField<d>& operator=(const BitFaceField<d>& copy)
	{mac_grid=copy.mac_grid;face_fields=copy.face_fields;return *this;}

	void Resize(const VectorDi& cell_counts)
	{
		if(mac_grid.grid.cell_counts==cell_counts)return;
		mac_grid.Initialize(cell_counts);
		for(int i=0;i<d;i++)face_fields[i].Resize(mac_grid.face_counts[i]);	
	}

	void Resize(const VectorDi& cell_counts,const bool& value)
	{
		if(mac_grid.grid.cell_counts==cell_counts)return;
		mac_grid.Initialize(cell_counts);
		for(int i=0;i<d;i++)face_fields[i].Resize(mac_grid.face_counts[i],value);
	}

	void Fill(const bool& value)
	{for(int i=0;i<d;i++)face_fields[i].Fill(value);}
	
	void Fill(const bool& value,const int axis)
	{face_fields[axis].Fill(value);}
	
	void Fill(const Vector<bool,d>& value)
	{for(int i=0;i<d;i++)Fill(value[i],i);}

    bool operator() (const int axis,const VectorDi& coord) const
	{return face_fields[axis](coord);}

	void Set(const int axis,const VectorDi& coord,const bool val=true)
	{face_fields[axis].Set(coord,val);}
};

template<int d,int vd> class BitVectorFaceField
{Typedef_VectorDi(d);
public:
	MacGrid<d> mac_grid;
	ArrayF<BitVectorField<d,vd>,d> face_fields;

	BitVectorFaceField(const VectorDi& cell_counts=VectorDi::Zero())
	{Resize(cell_counts);}
	
	BitVectorFaceField(const VectorDi& cell_counts,const bool& value)
	{Resize(cell_counts,value);}
	
	BitVectorFaceField(const BitVectorFaceField<d,vd>& copy)
	{*this=copy;}
	
	BitVectorFaceField<d,vd>& operator=(const BitVectorFaceField<d,vd>& copy)
	{mac_grid=copy.mac_grid;face_fields=copy.face_fields;return *this;}

	void Resize(const VectorDi& cell_counts)
	{
		if(mac_grid.grid.cell_counts==cell_counts)return;
		mac_grid.Initialize(cell_counts);
		for(int i=0;i<d;i++)face_fields[i].Resize(mac_grid.face_counts[i]);	
	}

	void Resize(const VectorDi& cell_counts,const bool& value)
	{
		if(mac_grid.grid.cell_counts==cell_counts)return;
		mac_grid.Initialize(cell_counts);
		for(int i=0;i<d;i++)face_fields[i].Resize(mac_grid.face_counts[i],value);
	}

	void Fill(const bool& value)
	{for(int i=0;i<d;i++)face_fields[i].Fill(value);}
	
	void Fill(const bool& value,const int axis)
	{face_fields[axis].Fill(value);}
	
	void Fill(const Vector<bool,d>& value)
	{for(int i=0;i<d;i++)Fill(value[i],i);}

    bool operator() (const int axis,const VectorDi& coord,const int v) const
	{return face_fields[axis](coord,v);}

	void Set(const int axis,const VectorDi& coord,const int v,const bool val=true)
	{face_fields[axis].Set(coord,v,val);}
};
#endif
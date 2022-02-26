//#####################################################################
// BitField
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//#####################################################################
#ifndef __BitField_h__
#define __BitField_h__
#include "Grid.h"

template<int d> class BitField
{Typedef_VectorDi(d);
public:
	VectorDi counts;
	Array<bool> array;

	////Constructors
	BitField(const VectorDi& _counts=VectorDi::Zero())
	{Resize(_counts);}
	
	BitField(const VectorDi& _counts,const bool value)
	{Resize(_counts,value);}
    
	BitField(const BitField<d>& copy)
	{*this=copy;}
    
	BitField<d>& operator=(const BitField<d>& copy)
	{counts=copy.counts;array=copy.array;return *this;}

	////Resize
	void Resize(const VectorDi& _counts)
	{counts=_counts;int n=counts.prod();if(n>0)array.resize(n);}
	
	void Resize(const VectorDi& _counts,const bool value)
	{counts=_counts;int n=counts.prod();if(n>0)array.resize(n,value);}
	
	void Resize(const int size)
	{counts=VectorDi::Ones();counts[0]=size;Resize(counts);}
	
	void Resize(const int size,const bool value)
	{counts=VectorDi::Ones();counts[0]=size;Resize(counts,value);}
    
	////element set
	inline void Set(const VectorDi& coord,bool val=true){array[Index(coord)]=val;}
	inline void Set(const int i,bool val=true){array[i]=val;}
    inline void Set(const int i,const int j,bool val=true){array[Index(i,j)]=val;}
	inline void Set(const int i,const int j,const int k,bool val=true){array[Index(i,j,k)]=val;}

	////element get
	inline bool operator()(const VectorDi& coord) const {return array[Index(coord)];} 
	inline bool operator()(const int i) const {return array[i];}
    inline bool operator()(const int i,const int j) const {return array[i*counts[1]+j];}
	inline bool operator()(const int i,const int j,const int k) const {return array[i*counts[1]*counts[2]+j*counts[2]+k];}\

	////global operations
	void Fill(const bool value)
	{std::fill(array.begin(),array.end(),value);}
	
	bool Has(const bool value) const 
	{auto find=std::find(array.begin(),array.end(),value);return find!=array.end();}

	////Helper functions
protected:
	inline int Index(const Vector1i& c) const {return c[0];}
	inline int Index(const Vector2i& c) const {return c[0]*counts[1]+c[1];}
	inline int Index(const Vector3i& c) const {return c[0]*counts[1]*counts[2]+c[1]*counts[2]+c[2];}
	inline int Index(const int i) const {return i;}
	inline int Index(const int i,const int j) const {return i*counts[1]+j;}
	inline int Index(const int i,const int j,const int k) const {return i*counts[1]*counts[2]+j*counts[2]+k;}
};

template<int d,int vd> class BitVectorField
{Typedef_VectorDi(d);
public:
	VectorDi counts;
	Array<bool> array;

	////Constructors
	BitVectorField(const VectorDi& _counts=VectorDi::Zero())
	{Resize(_counts);}
	
	BitVectorField(const VectorDi& _counts,const bool value)
	{Resize(_counts,value);}
    
	BitVectorField(const BitVectorField& copy)
	{*this=copy;}
    
	BitVectorField<d,vd>& operator=(const BitVectorField<d,vd>& copy)
	{counts=copy.counts;array=copy.array;return *this;}

	////Resize
	void Resize(const VectorDi& _counts)
	{counts=_counts;int n=counts.prod()*vd;if(n>0)array.resize(n);}
	
	void Resize(const VectorDi& _counts,const bool value)
	{counts=_counts;int n=counts.prod()*vd;if(n>0)array.resize(n,value);}
	
	void Resize(const int size)
	{counts=VectorDi::Ones();counts[0]=size;Resize(counts);}
	
	void Resize(const int size,const bool value)
	{counts=VectorDi::Ones();counts[0]=size;Resize(counts,value);}
    
	////element set
	inline void Set(const VectorDi& coord,const int v,bool val=true){array[Index(coord,v)]=val;}
	inline void Set(const int i,const int v,bool val=true){array[i*vd+v]=val;}
    inline void Set(const int i,const int j,const int v,bool val=true){array[Index(i,j,v)]=val;}
	inline void Set(const int i,const int j,const int k,const int v,bool val=true){array[Index(i,j,k,v)]=val;}

	////element get
	inline bool operator()(const VectorDi& coord,const int v) const {return array[Index(coord,v)];} 
	inline bool operator()(const int i,const int v) const {return array[i*vd+v];}
    inline bool operator()(const int i,const int j,const int v) const {return array[(i*counts[1]+j)*vd+v];}
	inline bool operator()(const int i,const int j,const int k,const int v) const {return array[(i*counts[1]*counts[2]+j*counts[2]+k)*vd+v];}

	////global operations
	void Fill(const bool value)
	{std::fill(array.begin(),array.end(),value);}
	
	bool Has(const bool value) const 
	{auto find=std::find(array.begin(),array.end(),value);return find!=array.end();}

	////Helper functions
protected:
	inline int Index(const Vector1i& c,const int v) const {return c[0]*vd+v;}
	inline int Index(const Vector2i& c,const int v) const {return (c[0]*counts[1]+c[1])*vd+v;}
	inline int Index(const Vector3i& c,const int v) const {return (c[0]*counts[1]*counts[2]+c[1]*counts[2]+c[2])*vd+v;}
	inline int Index(const int i,const int v) const {return i*vd+v;}
	inline int Index(const int i,const int j,const int v) const {return (i*counts[1]+j)*vd+v;}
	inline int Index(const int i,const int j,const int k,const int v) const {return (i*counts[1]*counts[2]+j*counts[2]+k)*vd+v;}
};

#endif

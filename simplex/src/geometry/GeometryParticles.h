//////////////////////////////////////////////////////////////////////////
// Geometry particles
// Copyright (c) (2018-),Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
//// Geometry particles refer to (dynamic) particles enhanced with local frames
//////////////////////////////////////////////////////////////////////////

#ifndef __GeometryParticles_h__
#define __GeometryParticles_h__

#include "Particles.h"

template<int d> class GeometryParticles : public Points<d>
{Typedef_VectorDii(d);Typedef_MatrixD(d);using Base=Points<d>;
public:
	//NOTE: seems these "using" are necessary for linux compilation - Mengdi
	using Base::Size; using Base::Resize; using Base::Add_Element; using Base::Add_Elements; using Base::Join; using Base::Copy_Element_From; using Base::Delete_Elements; using Base::Print_Attributes; using Base::Save_Snapshot; using Base::Load_Snapshot;
	using Base::X;	////use position from the base class
	using MatrixT=Matrix<real,d-1>;						////tangential matrix type
	using VectorT=Vector<real,d-1>;						////tangential vector type

	GeometryParticles(){New_Attributes();}

	Declare_Attribute(VectorD,V,v,Rebind_V);			////velocity
	Declare_Attribute(VectorD,F,f,Rebind_F);			////force
	Declare_Attribute(real,M,m,Rebind_M);				////mass
	Declare_Attribute(int,I,idx,Rebind_I);				////idx, -1 invalid particle

	Declare_Attribute(MatrixD,E,e,Rebind_E);			////local frame
	Declare_Attribute(MatrixT,G,g,Rebind_G);			////metric tensor
	Declare_Attribute(VectorT,dH,dh,Rebind_dH);			////dhdx

	Declare_Attribute_Inherent_Func(v,f,m,idx,e,g,dh);

	VectorD Normal(const int i) const {return E(i).col(d-1);}
};

#endif
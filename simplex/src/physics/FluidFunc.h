#ifndef __FluidFunc_h__
#define __FluidFunc_h__

#include "MacGrid.h"
#include "FaceField.h"
#include "Interpolation.h"

namespace FluidFunc
{
	//////////////////////////////////////////////////////////////////////////
	////vorticity
	void Curl_On_Cell(const MacGrid<2>& mac_grid,const FaceField<real,2>& v,Field<real,2>& vorticity);
	void Curl_On_Cell(const MacGrid<3>& mac_grid,const FaceField<real,3>& v,Field<Vector3,3>& vorticity);
	void Curl_On_Node(const MacGrid<2>& mac_grid,const FaceField<real,2>& v,Field<real,2>& vorticity);
	////TODO: 3D vorticity on edge
	//void Curl_On_Edge(const MacGrid<3>& mac_grid,const FaceField<real,3>& v,FaceField<real,3>& vorticity){}
	void Curl_On_Face(const MacGrid<2>& mac_grid,Field<real,2>& psi,FaceField<real,2>& vorticity);
	/////TODO: 3D vorticity on face
	template<int d> void Enforce_Vorticity_Confinement(const MacGrid<d>& mac_grid, const real dt, FaceField<real, d>& velocity, const real vor_conf_coef = 15);

	//////////////////////////////////////////////////////////////////////////
	////projection
	template<int d> void Correct_Velocity_With_Pressure(const MacGrid<d>& mac_grid, FaceField<real, d>& velocity, const Field<real, d>& pressure, const FaceField<real, d>& alpha, const Field<ushort, d>& type);

	//////////////////////////////////////////////////////////////////////////
	////diffusion
	////Diffusion on a field: this function implements Jos Stam's diffusion function in "REAL-TIME FLUID DYNAMICS FOR GAMES". 
	////The function solves the diffusion equation approximately with a fixed number of G-S iterations.
	////a=alpha*dt/(dx^2)
	////diffusion on a cell field
	template<class T,int d> void Diffusion(Field<T,d>& field,const Grid<d>& grid,const real a,const int iter_num/*=10*/);
	////diffusion on a face field
	template<int d> void Diffusion(FaceField<real,d>& field,const MacGrid<d>& mac_grid,const real a,std::function<real(const Vector<real,d>&)> Phi/*=nullptr*/,const int iter_num/*=10*/);

	//////////////////////////////////////////////////////////////////////////
	////Jacobian
	template<int d> void Jacobian(const FaceField<real,d>& field,const Interpolation<d>& intp,const Vector<real,d>& pos,
		/*rst*/Matrix<real,d>& jac,std::function<real(const Vector<real,d>&)> Phi=nullptr);

	////Deformation gradient
	template<int d> void Deformation_Gradient(const FaceField<real,d>& field,const Interpolation<d>& intp,const Vector<real,d>& pos,/*rst*/Matrix<real,d>& F);

};

#endif
//////////////////////////////////////////////////////////////////////////
// Point set
// Copyright (c) (2018-), Bo Zhu, Hui Wang, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __PointSet_h__
#define __PointSet_h__
#include "Common.h"
#include "GeometryParticles.h"
#include "TypeFunc.h"
#include "SpatialHashing.h"
#include "Kernels.h"
#include "LeastSquares.h"
#include "NeighborSearcher.h"

template<int d> class PointSet						////the second template parameter is the size of nb array for each particle
{Typedef_VectorDii(d);Typedef_MatrixD(d);
public:
	using VectorT=Vector<real,d-1>;								////tangential vector type
	using VectorTi=Vector<int,d-1>;								////tangential vector int
    using MatrixT=Matrix<real,d-1>;								////metric tensor matrix type

	GeometryParticles<d>* points=nullptr;
	bool own_data=false;
	enum class DifType : int {SPH=0,MPS,MLS};					////differential operator type
	DifType dif_type=DifType::SPH;								////default is SPH

	Heap<int,std::greater<int> > invalid_point_heap;			////for particle reseeding
	Array<int> reseeded_points;									////index array to record the reseeded points for the current frame

	Grid<d-1> t_grid;											////local tangential grid, currently not used
	KernelSPH v_kernel;										////volumetric SPH kernel
	KernelSPH t_kernel;									////tangential SPH kernel

	real v_dx=(real)1;											////dx for the global volumetric grid, initialized based on _dx
	real t_dx=(real)1;											////dx of the local tangential grid, initialized based on _dx
	real v_r=(real)1;											////volumetric radius, the supporting radius to search for the volumetric neighbors, initialized based on _dx			
	real t_r=(real)1;											////local tangential radius, the supporting radius to search for the tangential neighbors, initialized based on _dx
	real t_dot=(real).2;										////threshold for the dot product between two normals

private:
	Array<real> nden;											////number density, for all SPH-related functions
	real init_avg_nden = (real)0;								////initial number density
	real avg_nden = (real)0;									////current number density
public:
	std::shared_ptr<NeighborSearcher<d> > nbs_searcher;         ////to perform knn and radius search 
	bool use_kd_tree_nb_search=true;							////use spatial hashing if false
	Array<Array<int> > tang_nbs;                                ////tangential neighbors, data stored inside NeighborSearcher

public:
	//////////////////////////////////////////////////////////////////////////
	////Initialization
	PointSet(GeometryParticles<d>* _points=nullptr)
	{
		if(_points==nullptr&&points==nullptr){
			points=new GeometryParticles<d>();own_data=true;}
		else{if(_points!=nullptr){delete points;points=_points;own_data=false;}}	
	}
	~PointSet(){if(own_data&&points!=nullptr)delete points;}

	////The first parameter dx is the typical distance between two particles (dx for the tangential grid)
	////The second parameter t_grid_s is the number of cells for one-fourth of the tangential grid in one axis
	////e.g., if _t_grid_s=2, we project points to an 8x8 2D grid in the tangential space
	////supporing radius = dx*t_grid_s
	void Initialize(const real _dx,const int _t_grid_s,GeometryParticles<d>* _points=nullptr)
	{
		if (_points == nullptr && points == nullptr) {
			points = new GeometryParticles<d>(); 
			own_data = true;
		}
		else {
			if (_points != nullptr) {
				delete points;
				points = _points;
				own_data = false;
			}
		}
		t_dx=_dx;
		t_r=t_dx*(real)_t_grid_s;
		VectorTi t_cell_counts=VectorTi::Ones()*_t_grid_s*4;
		t_grid.Initialize(t_cell_counts,t_dx,-VectorT::Ones()*(real)_t_grid_s*t_dx);

		v_dx=t_r*(real)1.5;		////the volumetric grid dx is 1.5x of the tangential radius
		v_r=t_r*(real)1.5;		

		if(use_kd_tree_nb_search) nbs_searcher = std::make_shared<NeighborKDTree<d> >();
		else nbs_searcher = std::make_shared<NeighborHashing<d, 1024> >(v_dx);

		v_kernel = KernelSPH(v_r);
		t_kernel = KernelSPH(t_r);

		Update_Nbs();
		init_avg_nden=avg_nden=Calculate_Number_Density(points->XRef(),nden);	////initialize init_nden and nden

		std::cout<<"t_grid: "<<t_grid.cell_counts.transpose()<<", "<<t_grid.dx<<std::endl;
	}


	//////////////////////////////////////////////////////////////////////////
	////Interfaces

	////Global updates
	void Update()
	{
		Update_Nbs();
		Update_Number_Density();
		Reinitialize_Local_Frames();
		//Update_Metric_Tensor();
		//Print_Statistics();
	}

	////This function updates all the volumetric neighboring particles (not stored) and tangential neighboring particles (stored in tang_nbs)
	void Update_Nbs();

	////Tangential Neighbor
	inline const Array<int>& Tangential_Neighbor_Of(int idx) const { return tang_nbs[idx]; }

	////Number Density
	real Initial_Number_Density(void) { return init_avg_nden; }
	real Average_Number_Density(void) { return avg_nden; }
	real Number_Density_Of(int idx) { return nden[idx]; }
	
	//////////////////////////////////////////////////////////////////////////
	////Functions about neighbor search
	bool Is_Tangential_Neighbor(const VectorD& pos, const MatrixD& local_frame, const int p) const;

	////Find points that are within a circle with the radius of t_r on the tangential plane. 
	////Points with their directions that form a large angle with the local tangential plane will not be counted (t_dot is used to check the threshold).
	ArraySlice<int> Record_Tangential_Nbs(const VectorD& pos, const MatrixD& local_frame) const;
	Array<int> Find_Tangential_Nbs(const VectorD& pos, const MatrixD& local_frame) const;

	//////////////////////////////////////////////////////////////////////////
	////Local frame evolution
	////TOFIX: this function has some bugs
	void Update_Local_Frame(const real dt);
	////Initialize the local frame at v using PCA
	void Initialize_Local_Frame(const VectorD& v, const Array<int>& nbs, MatrixD& local_frame);
	////Reinitialization local frames using PCA
	void Reinitialize_Local_Frames();

	////Calculate normal of a point in the point set
	VectorD Normal(const int i) const { return points->E(i).col(d - 1); }
	////Calculate normal for an arbitrary position by PCA
	////The current implementation is to average the space
	VectorD Normal(const VectorD& pos) const;

	////Calculate local frame for an arbitrary position
	MatrixD Local_Frame(const VectorD& pos) const;

	//////////////////////////////////////////////////////////////////////////
	////Metric tensor
	////TOFIX: unnecessary after local frame correction. The metric tensor are identity.
	void Update_Metric_Tensor();

	////metric tensor in 2D
	Matrix<real, 1> Metric_Tensor(const Vector1& dzdx) const;
	////metric tensor in 3D
	Matrix<real, 2> Metric_Tensor(const Vector2& dzdx) const;

	//////////////////////////////////////////////////////////////////////////
	////Projection functions between global and local frames
	static VectorD Local_Coords(const VectorD& u, const MatrixD& e);//coords in local system e of world vector u
	static void Local_Coords(const VectorD& u, const MatrixD& e, VectorT& t, real& z);//(t,z). z is local h-axis
	static VectorT Project_To_TPlane(const VectorD& u, const MatrixD& e);
	static real Project_To_TPlane_H(const VectorD& u, const MatrixD& e);

	////Take tangential components
	////2D->3D projection, back to world space
	inline void Unproject_To_World(const Vector2& t, const Matrix3& e, Vector3& u) const { u = e.col(0) * t[0] + e.col(1) * t[1]; }
	////1D->2D projection, back to world space
	inline void Unproject_To_World(const Vector1& t, const Matrix2& e, Vector2& u) const { u = e.col(0) * t[0]; }

	////2D->3D projection including the local h-axis, back to world space
	inline void Unproject_To_World(const Vector3& t, const Matrix3& e, Vector3& u) const { u = e.col(0) * t[0] + e.col(1) * t[1] + e.col(2) * t[2]; }
	////1D->2D projection including the local h-axis, back to world space
	inline void Unproject_To_World(const Vector2& t, const Matrix2& e, Vector2& u) const { u = e.col(0) * t[0] + e.col(1) * t[1]; }

	//////////////////////////////////////////////////////////////////////////
	////Closest point
	inline int Closest_Point(const VectorD& pos) const { return nbs_searcher->Find_Nearest_Nb(pos); }
	inline VectorD Cloest_Normal(const VectorD& pos) const { int p = Closest_Point(pos); return Normal(p); }

	//////////////////////////////////////////////////////////////////////////
	////Implicit surface
	int Nearest_Geometry(VectorD& pos, MatrixD& frame, int minimum_eqns = 1, bool verbose=true)const;
	VectorD Project_To_Surface(const VectorD& pos) const;
	real Unsigned_Distance(const VectorD& pos) const;
	real Signed_Distance(const VectorD& pos) const;
	real Truncated_Unsigned_Distance(const VectorD& pos, const real& max_dist = -1, int minimum_eqns = 1)const;//max_dist<0 means no truncation
	real Truncated_Signed_Distance(const VectorD& pos, const real& max_dist = -1, int minimum_eqns = 1)const;

	//////////////////////////////////////////////////////////////////////////
	////Differential operator implementations
	//////////////////////////////////////////////////////////////////////////
	////Differential operators on an indexed point
	////Tangential gradient on an indexed point
	template<typename F> VectorT Grad_TPlane(const int p, const F& phi, int typ = -1)
	{
		DifType dt;
		if (typ == -1) dt = dif_type;
		else dt = (DifType)typ;
		switch (dt) {
		default:
		case DifType::SPH:return Grad_TPlane_SPH<F>(p, phi);
		case DifType::MLS:return Grad_TPlane_MLS<F>(p, phi);
		case DifType::MPS:return Grad_TPlane_MPS<F>(p, phi);}
	}

	////Tangential Laplacian on an indexed point
	template<typename F> real Laplacian_TPlane(const int p, const F& phi)
	{
		switch(dif_type){
		default:
		case DifType::SPH:return Laplacian_TPlane_SPH<F>(p,phi);
		case DifType::MLS:return Laplacian_TPlane_MLS<F>(p,phi);
		case DifType::MPS:return Laplacian_TPlane_MPS<F>(p,phi);}
	}

	template<typename F> VectorT Laplacian_VecT_TPlane(const int p,const F& phi)
	{
		switch(dif_type){
		default:
		case DifType::SPH:return Laplacian_VecT_TPlane_SPH<F>(p,phi);
		case DifType::MLS:return Laplacian_VecT_TPlane_MLS<F>(p,phi);
		case DifType::MPS:return Laplacian_VecT_TPlane_MPS<F>(p,phi);}
	}

	////////////////////////////////////////////////////////////////////////////
	//////Differential operators on an arbitrary position
	//////Tangential gradient on an arbitrary position 
	template<typename F> VectorT Grad_TPlane(const VectorD& pos,const MatrixD& lf,const F& phi)
	{
		switch(dif_type){
		default:
		case DifType::SPH:return VectorT::Zero();////TOIMPL
		case DifType::MLS:return Grad_TPlane_MLS<F>(pos,lf,phi);
		case DifType::MPS:return VectorT::Zero();}////TOIMPL
	}

	////Tangential Laplacian on an arbitrary position
	template<typename F> real Laplacian_TPlane(const VectorD& pos,const MatrixD& lf,const F& phi)
	{
		switch(dif_type){
		default:
		case DifType::SPH:return 0;////TOIMPL
		case DifType::MLS:return Laplacian_TPlane_MLS<F>(pos,lf,phi);
		case DifType::MPS:return 0;}////TOIMPL
	}

	template<typename F> VectorT Laplacian_VecT_TPlane(const VectorD& pos,const MatrixD& lf,const F& phi)
	{
		switch(dif_type){
		default:
		case DifType::SPH:return VectorT::Zero();////TOIMPL
		case DifType::MLS:return Laplacian_VecT_TPlane_MLS<F>(pos,lf,phi);
		case DifType::MPS:return VectorT::Zero();}////TOIMPL
	}

	//////////////////////////////////////////////////////////////////////////
	////SPH implementation
	inline real Vol(const int idx) const {return (real)1/nden[idx];}

	template<typename F> VectorT Grad_TPlane_SPH(const int i,const F& phi)
	{
		VectorT grad_phi=VectorT::Zero();
		int nbs_num = (int)tang_nbs[i].size();
		for(int k=0;k<nbs_num;k++){
			int j = tang_nbs[i][k];
			VectorT lr_ij = Project_To_TPlane(points->X(i) - points->X(j), points->E(i));
			// real lr2=lr_ij.squaredNorm();
			// real one_over_lr2=(lr2==(real)0?(real)0:(real)1/lr2);
			// VectorT grad=-(phi(i)*pow(Vol(i),2)+phi(j)*pow(Vol(j),2))*t_kernel.Grad_Spiky(lr_ij);
			// VectorT grad=phi(j)*pow(Vol(j),2)*t_kernel.Grad_Spiky(lr_ij);
			VectorT grad = phi(j) * Vol(j) * t_kernel.Grad<d - 1>(lr_ij, KernelType::SPIKY);
			grad_phi+=grad;}
		return grad_phi;
	}

	template<typename F> real Laplacian_TPlane_SPH(const int i,const F& phi)
	{
		real lap_phi=(real)0;
		int nbs_num = (int)tang_nbs[i].size();
		for(int k=0;k<nbs_num;k++){
			int j = tang_nbs[i][k];
			VectorT lr_ij = Project_To_TPlane(points->X(i) - points->X(j), points->E(i));
			real lr2=lr_ij.squaredNorm();
			real one_over_lr2=(lr2==(real)0?(real)0:(real)1/lr2);
			real phi_ij=phi(i)-phi(j);
			// if(i==j) continue;
			// real lap=Vol(i)*Vol(j)*phi_ij*one_over_lr2*lr_ij.dot(t_kernel.Grad_Spiky(lr_ij));
			////Eq 22 in "2019-Smoothed Particle Hydrodynamics Techniques for the Physics Based Simulation of Fluids and Solids"
			if(i==j) continue;
			real lap = -2 * Vol(j) * phi_ij * t_kernel.Grad<d - 1>(lr_ij, KernelType::CUBIC).norm() / sqrt(lr2);
			////Mesh-free Discrete Laplace–Beltrami Operator
			// if(i==j) continue;
			// real lap=2*Vol(j)*phi_ij*one_over_lr2*lr_ij.dot(t_kernel.Grad_Quintic(lr_ij));
			lap_phi+=lap;}
		return lap_phi;
	}

	template<typename F> VectorT Laplacian_VecT_TPlane_SPH(const int i,const F& u)
	{
		VectorT lap_u=VectorT::Zero();
		int nbs_num=tang_nbs[i].size();
		for(int k=0;k<nbs_num;k++){
			int j=tang_nbs[i][k];
			VectorT lr_ij;Project_To_TPlane(points->X(i)-points->X(j),points->E(i),lr_ij);
			VectorT u_ij=u(i)-u(j);
			// VectorT lap_ij=Vol(i)*Vol(j)*u_ij*one_over_lr2*lr_ij.dot(t_kernel.Grad_Spiky(lr_ij));
			// VectorT lap_ij=u(j)*pow(Vol(j),2)*t_kernel.Lap_Vis(lr_ij);
			if(i==j) continue;
			VectorT lap_ij=-2*Vol(j)*u_ij*t_kernel.Grad<d-1>(lr_ij,KernelType::CUBIC).norm()/lr_ij.norm();//sqrt(lr2);
			lap_u+=lap_ij;}
		return lap_u;
	}

	//////////////////////////////////////////////////////////////////////////
	////MLS implementation
	template<typename F> VectorT Grad_TPlane_MLS(const int p,const F& phi)
	{
		VectorX coeff=Fit_Local_MLS(p,phi);
		if constexpr(d==2) return VectorT(coeff(1));
		else if constexpr(d==3) return VectorT(coeff(1),coeff(2));
		return VectorT::Zero();
	}

	template<typename F> real Laplacian_TPlane_MLS(const int p,const F& phi)
	{
		VectorX coeff=Fit_Local_MLS(p,phi);
		if constexpr(d==2) return 2*coeff(2);
		else if constexpr(d==3) return 2*(coeff(3)+coeff(4));
		return (real)0;
	}

	template<typename F> VectorT Laplacian_VecT_TPlane_MLS(const int p,const F& u)
	{
		std::function<real(const int)> u0=[&](const int idx)->real{return u(idx)(0);};
		std::function<real(const int)> u1=[&](const int idx)->real{return u(idx)(1);};
		VectorX coeff0=Fit_Local_MLS(p,u0);
		VectorX coeff1=Fit_Local_MLS(p,u1);
		if constexpr(d==2) return (real)2.0*VectorT(coeff0(2));
		else if constexpr(d==3) return (real)2.0*VectorT(coeff0(3)+coeff0(4),coeff1(3)+coeff1(4));
		return VectorT::Zero();
	}
	
	//////////////////////////////////////////////////////////////////////////
	////MLS implementation on an arbitrary position
	template<typename F> VectorT Grad_TPlane_MLS(const VectorD& pos,const MatrixD& lf,const F& phi)
	{
		VectorX coeff=Fit_Local_MLS(pos,lf,phi);
		if constexpr(d==2) return VectorT(coeff(1));
		else if constexpr(d==3) return VectorT(coeff(1),coeff(2));
		return VectorT::Zero();
	}
	
	template<typename F> real Laplacian_TPlane_MLS(const VectorD& pos,const MatrixD& lf,const F& phi)
	{
		VectorX coeff=Fit_Local_MLS(pos,lf,phi);
		if constexpr(d==2) return 2*coeff(2);
		else if constexpr(d==3) return 2*(coeff(3)+coeff(4));
		return (real)0;
	}


	template<typename F> VectorT Laplacian_VecT_TPlane_MLS(const VectorD& pos,const MatrixD& lf,const F& u)
	{
		std::function<real(const int)> u0=[&](const int idx)->real{return u(idx)(0);};
		std::function<real(const int)> u1=[&](const int idx)->real{return u(idx)(1);};
		VectorX coeff0=Fit_Local_MLS(pos,lf,u0);
		VectorX coeff1=Fit_Local_MLS(pos,lf,u1);
		if constexpr(d==2) return (real)2.0*VectorT(coeff0(2));
		else if constexpr(d==3) return (real)2.0*VectorT(coeff0(3)+coeff0(4),coeff1(3)+coeff1(4));
		return VectorT::Zero();
	}
	
	//////////////////////////////////////////////////////////////////////////
	////MLS implementation
	template<typename F> VectorD Grad_MLS(const int p,const F& phi)
	{
		////following Codim-MLS fluid by Wang
		VectorX coeff=Fit_Local_MLS(p,phi);
		MatrixT g_inv=points->G(p).inverse();
		real g_det=points->G(p).determinant();
		real par[2];
		VectorD tangent[2];
		if constexpr(d==2){
			par[0]=coeff(1);
			tangent[0](0)=1;tangent[0](1)=coeff(1);}
		else if constexpr(d==3){
			par[0]=coeff(1);
			par[1]=coeff(2);
			tangent[0](0)=1;tangent[0](1)=0;tangent[0](2)=coeff(1);
			tangent[1](0)=0;tangent[1](1)=1;tangent[1](2)=coeff(2);}

		VectorD grad=VectorD::Zero();
		for(int tk=0;tk<d-1;tk++){
			for(int tl=0;tl<d-1;tl++){
				grad+=g_inv(tk,tl)*par[tl]*tangent[tk];}}
		return grad;
	}
	
	template<typename F> real Laplacian_MLS(const int p,const F& phi)
	{
		////following Codim-MLS fluid by Wang
		VectorX coeff=Fit_Local_MLS(p,phi);
		MatrixT g_inv=points->G(p).inverse();
		real g_det=points->G(p).determinant();
		MatrixT par2;
		if constexpr(d==2){
			par2(0)=2*coeff(2);}
		else if constexpr(d==3){
			par2(0,0)=2*coeff(3);
			par2(1,1)=2*coeff(4);
			par2(0,1)=par2(1,0)=coeff(5);}

		real lplc=(real)0.0;
		for(int tk=0;tk<d-1;tk++){
			for(int tl=0;tl<d-1;tl++){
				lplc+=g_inv(tk,tl)*par2(tl,tk);}}
		return lplc;
	}

	////helper function, fit phi at p-th point using MLS
	template<typename F> VectorX Fit_Local_MLS(const int p,const F& phi){
		using namespace LeastSquares;
		MatrixD local_frame=points->E(p);
		// MLS<d-1,2> mls;
		MLSPointSet<d-1> mls;

		size_t n = tang_nbs[p].size();
		Array<real> data(n * (size_t)d);
		VectorD x_p=points->X(p);
		int pi=0;
		for (size_t i = 0; i < n; i++) {
			int nb = tang_nbs[p][i];
			if (nb == p) pi = (int)i;
			VectorD th = Local_Coords(points->X(nb) - x_p, local_frame);
			for (size_t j = 0; j < d - 1; j++)data[i * (size_t)d + j] = th[j];
			data[i * (size_t)d + (size_t)d - (size_t)1] = phi(nb);
		}		
		mls.Fit(&data[0], (int)n, pi);		
		return mls.c;
	}

	////helper function, fit local geometry at p-th point using MLS
	VectorX Fit_Local_Geometry_MLS(const int p){
		using namespace LeastSquares;
		MatrixD local_frame=points->E(p);
		// MLS<d-1,2> mls;
		MLSPointSet<d-1> mls;

		int n = (int)tang_nbs[p].size();
		Array<real> data(n*d);
		VectorD x_p=points->X(p);
		int pi=0;
		for(int i=0;i<n;i++){
			int nb = tang_nbs[p][i];
			if(nb==p) pi=i;
			VectorD th = Local_Coords(points->X(nb) - x_p, local_frame);
			for(int j=0;j<d;j++)data[i*d+j]=th[j];}
		mls.Fit(&data[0],n,pi);
		return mls.c;
	}

	////helper function, fit phi on an arbitrary position using MLS
	template<typename F> VectorX Fit_Local_MLS(const VectorD& pos,const MatrixD& lf,const F& phi){
		using namespace LeastSquares;
		Array<int> tang_nbs = Find_Tangential_Nbs(pos, lf);
		MLSPointSet<d-1> mls;
		size_t n = tang_nbs.size();
		if(n==0) std::cout<<"Warning: no neighbor around"<<pos.transpose()<<std::endl;
		
		////TOFIX? discontinuous weight function
		Array<real> data(n * (size_t)d);
		int pi=-1;
		for (size_t i = 0; i < n; i++) {
			int nb = tang_nbs[i];
			VectorD th = Local_Coords(points->X(nb) - pos, lf);
			for (size_t j=0;j<d-1;j++)data[i*(size_t)d+j]=th[j];
			data[i*(size_t)d+(size_t)d-(size_t)1] = phi(nb);
		}		
		mls.Fit(&data[0], (int)n, pi);//// degenerate to LS
		return mls.c;
	}

	////helper function, fit local geometry on an arbitrary position using MLS
	VectorX Fit_Local_Geometry_MLS(const VectorD& pos, const MatrixD& lf){
		using namespace LeastSquares;
		Array<int> tang_nbs = Find_Tangential_Nbs(pos, lf);
		MLSPointSet<d-1> mls;
		int n=(int)tang_nbs.size();
		if(n==0) std::cout<<"Warning: no neighbor around"<<pos.transpose()<<std::endl;

		
		////TOFIX? discontinuous weight function
		Array<real> data(n*d);
		int pi=-1;
		for (size_t i = 0; i < n; i++) {
			int nb = tang_nbs[i];
			VectorD th = Local_Coords(points->X(nb) - pos, lf);
			for(int j=0;j<d;j++)data[i*d+j]=th[j];
		}		
		mls.Fit(&data[0], (int)n, pi);//// degenerate to LS
		return mls.c;

		////continuous weight function
		// Array<real> data(n*d);
		// Array<real> weights(n);
		// for(int i=0;i<n;i++){
		// 	int nb = tang_nbs[i];
		// 	VectorD th;Project_To_TPlane(points->X(nb)-pos,lf,th);
		// 	for(int j=0;j<d;j++)data[i*d+j]=th[j];
		// 	real frac=(points->X(nb)-pos).norm(); weights[i]=pow(1-frac,4)*(4*frac+1);}
		// mls.Fit(&data[0],&weights[0],n);
		// return mls.c;
	}

	//////////////////////////////////////////////////////////////////////////
	////MPS implementation
	template<typename F> VectorT Grad_TPlane_MPS(const int p,const F& phi)
	{
		VectorT grad=VectorT::Zero();
		int nbs_num = (int)tang_nbs[p].size();
		real nden=(real)0;
		for(int i=0;i<nbs_num;i++){
			int nb = tang_nbs[p][i];
			if(p==nb)continue;
			VectorD u=points->X(nb)-points->X(p);
			VectorT t = Project_To_TPlane(u, points->E(p));
			real r=t.norm();
			real w=W_MPS(abs(r));
			VectorT dir=t/r;
			grad+=(phi(nb)-phi(p))/r*dir*w;
			nden+=w;}
		if(nden!=(real)0)grad*=((real)d/nden);
		return grad;
	}

	template<typename F> real Laplacian_TPlane_MPS(const int p,const F& phi)
	{
		real lap=(real)0;
		int nbs_num = (int)tang_nbs[p].size();
		real nden=(real)0;
		for(int i=0;i<nbs_num;i++){
			int nb=tang_nbs[p][i];
			if(p==nb)continue;
			VectorD u=points->X(nb)-points->X(p);
			VectorT t = Project_To_TPlane(u, points->E(p));
			real r=t.norm();
			real w=W_MPS(abs(r));
			lap+=(phi(nb)-phi(p))*w;
			nden+=(r*r)*w;}
		if(nden!=(real)0)lap*=((real)d*(real)2/nden);
		return lap;	
	}

	template<typename F> VectorT Laplacian_VecT_TPlane_MPS(const int p,const F& u)
	{
		VectorT lap=VectorT::Zero();
		int nbs_num = tang_nbs[p].size();
		real nden=(real)0;
		for(int i=0;i<nbs_num;i++){
			int nb = tang_nbs[p][i];
			if(p==nb)continue;
			VectorD u=points->X(nb)-points->X(p);
			VectorT t;Project_To_TPlane(u,points->E(p),t);
			real r=t.norm();
			real w=W_MPS(abs(r));
			lap+=(u(nb)-u(p))*w;
			nden+=(r*r)*w;}
		if(nden!=(real)0)lap*=((real)d*(real)2/nden);
		return lap;	
	}

	//////////////////////////////////////////////////////////////////////////
	////Differential operators on surface
	////Operators on i-th point
	template<typename F> VectorD Grad(const int i,const F& phi) //Return the global gradient vector
	{
		VectorT grad=Grad_TPlane(i,phi);
		if constexpr (d==2){
			return points->E(i)*Vector2(grad[0],0);}
		else if constexpr (d==3){
			return points->E(i)*Vector3(grad[0],grad[1],0);}
	}

	template<typename F> real Div(const int i, const F& phi, int typ = -1)	//Here F is a global vector field
	{
		if constexpr (d == 3) {
			VectorT e1, e2; e1 = VectorT::Unit(0); e2 = VectorT::Unit(1);
			// std::function<real(const int)> sqrt_gu1=[&](const int idx)->real{return sqrt(Local_Metric_Tensor(i,idx).determinant())*phi(i)[0];};
			// std::function<real(const int)> sqrt_gu2=[&](const int idx)->real{return sqrt(Local_Metric_Tensor(i,idx).determinant())*phi(i)[1];};
			std::function<real(const int)> sqrt_gu1 = [&](const int idx)->real {
				VectorT tv; Project_To_TPlane(phi(idx), points->E(i), tv);
				return sqrt(Local_Metric_Tensor(i, idx).determinant()) * tv(0); };
			std::function<real(const int)> sqrt_gu2 = [&](const int idx)->real {
				VectorT tv; Project_To_TPlane(phi(idx), points->E(i), tv);
				return sqrt(Local_Metric_Tensor(i, idx).determinant()) * tv(1); };
			real sqrt_gu1_par1 = Grad_TPlane(i, sqrt_gu1, typ).dot(e1);
			real sqrt_gu2_par2 = Grad_TPlane(i, sqrt_gu2, typ).dot(e2);
			return (sqrt_gu1_par1 + sqrt_gu2_par2);
		}
		////TODO:what if det(g)=0?}
		else if constexpr (d == 2) {
			VectorT e1 = VectorT::Unit(0);
			// std::function<real(const int)> sqrt_gu1=[&](const int idx)->real{return sqrt(Local_Metric_Tensor(i,idx).determinant())*phi(i)[0];};
			std::function<real(const int)> sqrt_gu1 = [&](const int idx)->real {
				VectorT tv; Project_To_TPlane(phi(idx), points->E(i), tv);
				return sqrt(Local_Metric_Tensor(i, idx).determinant()) * tv(0); };
			real sqrt_gu1_par1 = Grad_TPlane(i, sqrt_gu1, typ).dot(e1);
			return sqrt_gu1_par1;
		}
	}

	template<typename F> real Laplacian(const int i,const F& phi)
	{
		if constexpr (d==3){
			//// Applying the product rule:
			//// par(sqrt(g)*g^kl*(par(s)/par(x_l)))/par xk
			////   =par(sqrt(g)*g^kl)/par(xk)*par(s)/par(x_l)
			////   +sqrt(g)*g^kl*par^2(s)/par(x_l)/par(x_k);
			VectorT grad=Grad_TPlane(i,phi);
			VectorT e0,e1;e0=VectorT::Unit(0);e1=VectorT::Unit(1);

			std::function<real(const int)> sqrtg_g00=[&](const int idx)->real{ 
					MatrixT g=Local_Metric_Tensor(i,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(0,0);};
            real sqrtg_g00_par0=Grad_TPlane(i,sqrtg_g00).dot(e0);

            std::function<real(const int)> sqrtg_g01=[&](const int idx)->real{
					MatrixT g=Local_Metric_Tensor(i,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(0,1);};
            real sqrtg_g01_par0=Grad_TPlane(i,sqrtg_g01).dot(e0);
            
			std::function<real(const int)> sqrtg_g10=[&](const int idx)->real{
					MatrixT g=Local_Metric_Tensor(i,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(1,0);};
			real sqrtg_g10_par1=Grad_TPlane(i,sqrtg_g10).dot(e1);
            
			std::function<real(const int)> sqrtg_g11=[&](const int idx)->real{
					MatrixT g=Local_Metric_Tensor(i,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(1,1);};
			real sqrtg_g11_par1=Grad_TPlane(i,sqrtg_g11).dot(e1);

			return grad(0)*(sqrtg_g00_par0+sqrtg_g10_par1)+
				grad(1)*(sqrtg_g01_par0+sqrtg_g11_par1)+Laplacian_TPlane(i,phi);
		}else if constexpr (d==2){
			VectorT grad=Grad_TPlane(i,phi);
			VectorT e0=VectorT::Unit(0);

			std::function<real(const int)> sqrtg_g00=[&](const int idx)->real{ 
					MatrixT g=Local_Metric_Tensor(i,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(0,0);};
            real sqrtg_g00_par0=Grad_TPlane(i,sqrtg_g00).dot(e0);
			
			return grad(0)*sqrtg_g00_par0+Laplacian_TPlane(i,phi);
		}
		return (real)0;
	}



	//////////////////////////////////////////////////////////////////////////
	////Differential operators on surface
	////Operators on an arbitrary position
	template<typename F> VectorD Grad(const VectorD& pos, const MatrixD& lf, const F& phi)
	{
		VectorT grad=Grad_TPlane(pos,lf,phi);
		if constexpr (d==2){
			return lf*Vector2(grad[0],0);}
		else if constexpr (d==3){
			return lf*Vector3(grad[0],grad[1],0);}
	}

	template<typename F> real Div(const VectorD& pos, const MatrixD& lf, const F& phi, int typ = -1)
	{
		if constexpr (d == 3) {
			VectorT e1, e2; e1 = VectorT::Unit(0); e2 = VectorT::Unit(1);
			// std::function<real(const int)> sqrt_gu1=[&](const int idx)->real{return sqrt(Local_Metric_Tensor(i,idx).determinant())*phi(i)[0];};
			// std::function<real(const int)> sqrt_gu2=[&](const int idx)->real{return sqrt(Local_Metric_Tensor(i,idx).determinant())*phi(i)[1];};
			std::function<real(const int)> sqrt_gu1 = [&](const int idx)->real {
				VectorT tv; Project_To_TPlane(phi(idx), lf, tv);
				return sqrt(Local_Metric_Tensor(lf, idx).determinant()) * tv(0); };
			std::function<real(const int)> sqrt_gu2 = [&](const int idx)->real {
				VectorT tv; Project_To_TPlane(phi(idx), lf, tv);
				return sqrt(Local_Metric_Tensor(lf, idx).determinant()) * tv(1); };
			real sqrt_gu1_par1 = Grad_TPlane(pos,lf,sqrt_gu1, typ).dot(e1);
			real sqrt_gu2_par2 = Grad_TPlane(pos,lf,sqrt_gu2, typ).dot(e2);
			return (sqrt_gu1_par1 + sqrt_gu2_par2);
		}
		else if constexpr (d == 2) {
			VectorT e1 = VectorT::Unit(0);
			// std::function<real(const int)> sqrt_gu1=[&](const int idx)->real{return sqrt(Local_Metric_Tensor(lf,idx).determinant())*phi(i)[0];};
			std::function<real(const int)> sqrt_gu1 = [&](const int idx)->real {
				VectorT tv; Project_To_TPlane(phi(idx), lf, tv);
				return sqrt(Local_Metric_Tensor(lf, idx).determinant()) * tv(0); };
			real sqrt_gu1_par1 = Grad_TPlane(pos,lf,sqrt_gu1, typ).dot(e1);
			return sqrt_gu1_par1;
		}
	}

	template<typename F> real Laplacian(const VectorD& pos,const MatrixD& lf,const F& phi)
	{
		if constexpr (d==3){
			//// Applying the product rule:
			//// par(sqrt(g)*g^kl*(par(s)/par(x_l)))/par xk
			////   =par(sqrt(g)*g^kl)/par(xk)*par(s)/par(x_l)
			////   +sqrt(g)*g^kl*par^2(s)/par(x_l)/par(x_k);
			VectorT grad=Grad_TPlane(pos,lf,phi);
			VectorT e0,e1;e0=VectorT::Unit(0);e1=VectorT::Unit(1);

			std::function<real(const int)> sqrtg_g00=[&](const int idx)->real{ 
					MatrixT g=Local_Metric_Tensor(lf,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(0,0);};
            real sqrtg_g00_par0=Grad_TPlane(pos,lf,sqrtg_g00).dot(e0);

            std::function<real(const int)> sqrtg_g01=[&](const int idx)->real{
					MatrixT g=Local_Metric_Tensor(lf,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(0,1);};
            real sqrtg_g01_par0=Grad_TPlane(pos,lf,sqrtg_g01).dot(e0);
            
			std::function<real(const int)> sqrtg_g10=[&](const int idx)->real{
					MatrixT g=Local_Metric_Tensor(lf,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(1,0);};
			real sqrtg_g10_par1=Grad_TPlane(pos,lf,sqrtg_g10).dot(e1);
            
			std::function<real(const int)> sqrtg_g11=[&](const int idx)->real{
					MatrixT g=Local_Metric_Tensor(lf,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(1,1);};
			real sqrtg_g11_par1=Grad_TPlane(pos,lf,sqrtg_g11).dot(e1);

			return grad(0)*(sqrtg_g00_par0+sqrtg_g10_par1)+
				grad(1)*(sqrtg_g01_par0+sqrtg_g11_par1)+Laplacian_TPlane(pos,lf,phi);
		}else if constexpr (d==2){
			VectorT grad=Grad_TPlane(pos,lf,phi);
			VectorT e0=VectorT::Unit(0);

			std::function<real(const int)> sqrtg_g00=[&](const int idx)->real{ 
					MatrixT g=Local_Metric_Tensor(lf,idx); MatrixT g_inv=g.inverse();
					return sqrt(g.determinant())*g_inv(0,0);};
            real sqrtg_g00_par0=Grad_TPlane(pos,lf,sqrtg_g00).dot(e0);
			
			return grad(0)*sqrtg_g00_par0+Laplacian_TPlane(pos,lf,phi);
		}
		return (real)0;
	}

	VectorD Curvature_Vector(const int i);
	VectorD Curvature_Vector(const VectorD& pos, const MatrixD& lf);


	//// helper function, return j's metric tensor in i's local frame
	MatrixT Local_Metric_Tensor(const int i, const int j) { return Metric_Tensor(Local_Dzdx(i, j)); }
	//// helper function, return j's metric tensor in local frame lf
	MatrixT Local_Metric_Tensor(const MatrixD& lf, const int j) { return Metric_Tensor(Local_Dzdx(lf, j)); }

	//// helper function, return j's dzdx in i's local frame
	VectorT Local_Dzdx(const int i, const int j);
	//// helper function, return j's dzdx in local frame lf
	VectorT Local_Dzdx(const MatrixD& lf, const int j);

	//////////////////////////////////////////////////////////////////////////
	////Kernel functions
	////MPS kernel function
	real W_MPS(const real r) const
	{
		if(r<t_r)return t_r/r-(real)1;
		else return (real)0;
	}

	////Kernel for PCA normal calculation
	real W_PCA(const real r) const
	{
		if(r<v_r)return (real)1-pow(r/v_r,3);
		else return (real)0;	
	}

	//////////////////////////////////////////////////////////////////////////
	////Point reseeding, deletion, and relaxation with SPH
	////Note: These functions should not change any point attributes
	void Update_Number_Density() { avg_nden = Calculate_Number_Density(points->XRef(), nden); }
	real Calculate_Number_Density(const Array<VectorD>& X, Array<real>& nden);
	real Calculate_Number_Density(const Array<VectorD>& X, const VectorD& pos) const;
	
	////Before reseeding, we assume nden is already updated.
	////Currently we do not delete points
	void Point_Reseeding();

	////This function redistribute the point position on the surface by applying SPH forces.
	////The SPH calculations assume a fixed neighboring relations. But it allows floating positions.
	////TOFIX: These parameters do not work for different dx yet.
	void Point_Relaxation();


	real Ray_Casting_MLS(VectorD ray_p, VectorD ray_d, real max_distance, MatrixD& lf_o, VectorX& c_o);

	//// Solve 
	void Solve_Ray_Sphere_Intersection(const VectorD& p, const VectorD& dir, const VectorD& s, const real r, real& norm_dist, real& tang_dist, real& delta_dist); 
	real Solve_Intersection(const VectorD& p, const VectorD& dir, const VectorX& coeff); 

	//////////////////////////////////////////////////////////////////////////
	////Statistics
	void Print_Statistics();
	

protected:
	//////////////////////////////////////////////////////////////////////////
	////particle system operations
	inline bool Valid_Particle(const int i) const {return points->I(i)!=-1;}
	inline int Append_Particles(int n) {return points->Add_Elements(n);}
	inline int Append_Particle() {return points->Add_Element();}
	int Add_Particle(){if(!invalid_point_heap.empty()){int p=invalid_point_heap.top();invalid_point_heap.pop();points->I(p)=0;return p;}else return Append_Particle();}
	void Remove_Particle(const int i){invalid_point_heap.push(i);points->I(i)=-1;points->V(i)=VectorD::Zero();points->X(i)=VectorD::Zero();}
};

#endif

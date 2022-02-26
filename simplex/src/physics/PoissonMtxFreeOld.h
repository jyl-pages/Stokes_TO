//////////////////////////////////////////////////////////////////////////
// Poisson solver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __PoissonMtxFreeOld_h__
#define __PoissonMtxFreeOld_h__
#include "Common.h"
#include "Field.h"
#include "MacGrid.h"
#include "FaceField.h"
#include "BitFaceField.h"
#include "Interpolation.h"
#include "KrylovSolver.h"
#include "Hashtable.h"
#include "Timer.h"
#include "MultiGridMtxFree.h"

namespace KrylovPoisson
{
using namespace MultiGrid;
template<int d,class T=real> class KrylovPoissonX
{using KXW=KrylovPoissonX<d,T>;using KX=Field<T,d>;Typedef_VectorDii(d);
public:
	Field<T,d>* kx=nullptr;
	int n=0;
	bool own_data=true;

	KrylovPoissonX(KX* _kx):kx(_kx)
	{own_data=false;n=kx->counts.prod();}
	
	KrylovPoissonX(const VectorDi& counts)
	{if(kx!=nullptr&&own_data)delete kx;kx=new KX();kx->Resize(counts,Zero<T>());n=kx->counts.prod();}
	~KrylovPoissonX(){if(own_data){delete kx;}}
	
	//////////////////////////////////////////////////////////////////////////
	////Krylov operations
	////*this=ax+y
	void Axpy(real a,const KXW& x,const KXW& y)
	{
		for(int i=0;i<n;i++){
			(*kx).array[i]=a*(*x.kx).array[i]+(*y.kx).array[i];}
	}

	////*this=x
	void Copy(const KXW& x)
	{*kx=*x.kx;}

	real Dot(const KXW& x) const
	{
		T dot=(T)0;
		for(int i=0;i<n;i++){
			dot+=(*kx).array[i]*(*x.kx).array[i];}
		return dot;
	}

	void Set_Zero()
	{kx->Fill(Zero<T>());}

	void Print()
	{
		std::cout<<"n: "<<n<<":  ";
		for(int i=0;i<n;i++)std::cout<<(*kx).array[i]<<", ";std::cout<<std::endl;
	}

	//////////////////////////////////////////////////////////////////////////
	////Multigrid functions
	void Restriction(KrylovPoissonX& coarse_kx,const Vector<int,d>& factor)
	{
		MultiGrid::Restriction<d,T>((*kx),(*coarse_kx.kx),factor);
	}

	void Prolongation(KrylovPoissonX& fine_kx,const Vector<int,d>& factor)
	{
		MultiGrid::Prolongation<d,T>((*fine_kx.kx),(*kx),factor);
	}
};

template<int d,class T=real> class KrylovPoissonA
{using KXW=KrylovPoissonX<d,T>;using KX=Field<T,d>;Typedef_VectorDii(d);
public:
	MacGrid<d> mac_grid;
	FaceField<T,d>* alpha=nullptr;
	BitField<d>* psi_D=nullptr;					////psi_D flags	
	Field<T,d>* psi_D_values=nullptr;				////psi_D values
	BitFaceField<d>* psi_N=nullptr;				////psi_N flags
	FaceField<T,d>* psi_N_values=nullptr;			////psi_N values
	bool own_data=false;

	KrylovPoissonA(const VectorDi& counts,FaceField<T,d>* _alpha,
		BitField<d>* _psi_D,Field<T,d>* _psi_D_values,
		BitFaceField<d>* _psi_N,FaceField<T,d>* _psi_N_values)
	{
		mac_grid.Initialize(counts,(real)1);
		alpha=_alpha;psi_D=_psi_D;psi_D_values=_psi_D_values;psi_N=_psi_N;psi_N_values=_psi_N_values;
		own_data=false;
	}

	KrylovPoissonA(const VectorDi& counts)
	{
		mac_grid.Initialize(counts,(real)1);
		alpha=new FaceField<T,d>(counts,(real)1);
		psi_D=new BitField<d>(counts,false);
		psi_D_values=new Field<T,d>(counts,(real)0);
		psi_N=new BitFaceField<d>(counts,false);
		psi_N_values=new FaceField<T,d>(counts,(real)0);

		own_data=true;
	}

	~KrylovPoissonA()
	{if(own_data){delete alpha;delete psi_D;delete psi_D_values;delete psi_N;delete psi_N_values;}}

	//////////////////////////////////////////////////////////////////////////
	////Krylov opreations
	////rst=Ax, A is negative Laplacian multiplying delta x squared
	void Multiply(const KXW& x,KXW& rst) const 
	{
		int n=x.kx->counts.prod();
		#pragma omp parallel for
		for(int i=0;i<n;i++){
			const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
			////psi_D cell
			if((*psi_D)(cell)){(*rst.kx).array[i]=(*x.kx).array[i];continue;}
			////non-psi_D cell
			T Ax=Zero<T>();
			for(int j=0;j<Grid<d>::Number_Of_Nb_C();j++){
				int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(j,axis,side);
				const VectorDi nb_cell=Grid<d>::Nb_C(cell,j);
				real nb_val=(real)0;real a=(real)1;
				if(mac_grid.grid.Valid_Cell(nb_cell)){
					VectorDi face=mac_grid.Cell_Incident_Face(axis,cell,side);
					a=(*alpha)(axis,face);
					if(!(*psi_N)(axis,face)&&!(*psi_D)(nb_cell))nb_val=(*x.kx)(nb_cell);}
				Ax+=a*((*x.kx)(cell)-nb_val);}
			(*rst.kx)(cell)=Ax;}
		//rst.Print();
	}

	////rst=b-Ax
	void Res(const KXW& b,const KXW& x,KXW& rst) const 
	{
		Multiply(x,rst);
		rst.Axpy((real)-1.,rst,b);
	}		

	//////////////////////////////////////////////////////////////////////////
	////Multigrid operations
	void Gauss_Seidel_Smoothing(KrylovPoissonX<d,T>& x,const KrylovPoissonX<d,T>& b,const int iter_num=1,const bool reversed_order=false)
	{
		//std::cout<<"gs: "<<x.kx->counts.transpose()<<std::endl;
		//KrylovPoissonX<d,T> r0(x.kx->counts);
		//Res(b,x,r0);
		//std::cout<<"res 0: "<<r0.Dot(r0)<<std::endl;

		int n=x.kx->counts.prod();
		for(int iter=0;iter<iter_num;iter++){
			for(int ii=0;ii<n;ii++){
				int i=(reversed_order?n-1-ii:ii);
				const VectorDi& cell=mac_grid.grid.Cell_Coord(i);
				if((*psi_D)(cell)){(*x.kx).array[i]=(*b.kx).array[i];continue;}

				T new_x=(*b.kx).array[i];
				real dia_coef=(real)0;
				for(int j=0;j<Grid<d>::Number_Of_Nb_C();j++){
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(j,axis,side);
					const VectorDi nb_cell=Grid<d>::Nb_C(cell,j);
					if(!mac_grid.grid.Valid_Cell(nb_cell)){
						dia_coef+=(real)1;
						continue;}
					VectorDi face=mac_grid.Cell_Incident_Face(axis,cell,side);
					if((*psi_N)(axis,face))continue;
					real nb_val=(real)0;
					if(!(*psi_D)(nb_cell))nb_val=(*x.kx)(nb_cell);
					new_x+=nb_val;
					dia_coef+=(real)1;}
				if(dia_coef>(real)0)new_x/=dia_coef;
				(*x.kx).array[i]=new_x;}}

		//Res(b,x,r0);
		//std::cout<<"res 1: "<<r0.Dot(r0)<<std::endl;
	}

	void Solve(KrylovPoissonX<d,T>& kx,const KrylovPoissonX<d,T>& kb)
	{
		using namespace KrylovPoisson;
		Vector<int,d> counts=mac_grid.grid.cell_counts;
		KrylovPoissonX<d,real> kr(counts),kd(counts),kq(counts),ks(counts);
		KrylovPoissonPrecIdentity<d,real> kP;

		KrylovSolver::Preconditioned_Conjugate_Gradient((*this),kx,kb,kr,kd,kq,ks,kP);		
	}
};

////Identity preconditioner
template<int d,class T=real> class KrylovPoissonPrecIdentity
{using KXW=KrylovPoissonX<d,T>;
public:
	KrylovPoissonPrecIdentity(){Initialize();}
	void Initialize(){}
	template<class KXW> void Precond(const KXW& r,KXW& rst){rst.Copy(r);}
};

using namespace MultiGrid;

////Geometric multigrid preconditioner
template<int d,class T,class KA,class KX> class KrylovPoissonPrecGMG
{public:
	const KA* ka=nullptr;
	MultiGrid::MultiGridSystemMtxFree<d,T,KA,KX> gmg;
	bool data_initialized=false;

	KrylovPoissonPrecGMG(const KA* _ka):ka(_ka){Initialize();}
	~KrylovPoissonPrecGMG()
	{
		if(data_initialized){
			for(int i=0;i<gmg.params.levels;i++){
				if(i>0)delete gmg.A[i];
				delete gmg.x[i];
				delete gmg.b[i];
				delete gmg.r[i];}}
	}

	template<int d> Vector<int,d> Coarsened(const Vector<int,d>& counts,const int coarsen_factor)
	{Vector<int,d> coarsened;for(int i=0;i<d;i++){coarsened[i]=counts[i]/coarsen_factor;if(coarsened[i]==0)coarsened[i]=1;}return coarsened;}

	void Initialize()
	{
		////TODO: one_over_intp=(T)1/params.relax_coef;
		const MacGrid<d>& mac_grid=ka->mac_grid;
		Vector<int,d> counts=mac_grid.grid.cell_counts;

		auto& A=gmg.A;
		auto& b=gmg.b;
		auto& x=gmg.x;
		auto& r=gmg.r;
		auto& params=gmg.params;
		auto& coarsen_factor=gmg.coarsen_factor;

		int levels=0;
		if(params.use_auto_calculated_levels){
			levels=(int)std::ceil(std::log2(AuxFunc::Max(counts)))-2;
			levels=std::max(1,levels);
			std::cout<<"Auto mg level, #levels: "<<levels<<std::endl;}
		else levels=params.levels;
		//levels=4;
		params.levels=levels;
		params.smooth_num=4;
		A.resize(levels);for(int i=1;i<levels;i++)A[i]=nullptr;
		b.resize(levels);for(int i=0;i<levels;i++)b[i]=nullptr;
		x.resize(levels);for(int i=0;i<levels;i++)x[i]=nullptr;
		r.resize(levels);for(int i=0;i<levels;i++)r[i]=nullptr;
		coarsen_factor.resize(levels);

		//{
		//	Field<real,d> fine(counts,(real)1);
		//	Field<real,d> coarse(counts/2,(real)0);
		//	MultiGrid::Restriction<d,real>(fine,coarse,Vector2i::Ones()*2);
		//	auto c=coarse.counts;
		//	for(int ii=0;ii<c[0];ii++){
		//		for(int jj=0;jj<c[1];jj++){
		//			std::cout<<coarse(Vector2i(ii,jj))<<", ";}
		//		std::cout<<std::endl;}
		//	std::cout<<"coarse end"<<std::endl;		
		//}

		//{
		//	Field<real,d> fine(counts,(real)0);
		//	Field<real,d> coarse(counts/2,(real)1);
		//	MultiGrid::Prolongation<d,real>(coarse,fine,Vector2i::Ones()*2);
		//	auto c=fine.counts;
		//	for(int ii=0;ii<c[0];ii++){
		//		for(int jj=0;jj<c[1];jj++){
		//			std::cout<<fine(Vector2i(ii,jj))<<", ";}
		//		std::cout<<std::endl;}
		//	std::cout<<"fine end"<<std::endl;
		//}

		for(int i=0;i<levels;i++){
			if(i==0){
				A[0]=const_cast<KA*>(ka);
			}
			else{
				Vector<int,d> old_counts=counts;
				counts=Coarsened<d>(old_counts,params.coarsen_factor);
				Vector<int,d> factor=old_counts.cwiseQuotient(counts);
				coarsen_factor[i-1]=factor;
				A[i]=new KA(counts);
				Restriction_Boundary_Condition<d,T>((*A[i-1]->psi_D),(*A[i]->psi_D),(*A[i-1]->psi_D_values),(*A[i]->psi_D_values),factor);
				////TODO: psi_N
				////may need alpha also
				
				//auto c=(*A[i]->psi_D).counts;
				//std::cout<<"counts: "<<c.transpose()<<std::endl;
				//for(int ii=0;ii<c[0];ii++){
				//	for(int jj=0;jj<c[1];jj++){
				//		std::cout<<(*A[i]->psi_D_values)(Vector2i(ii,jj))<<", ";}
				//	std::cout<<std::endl;}

				std::cout<<"counts: "<<counts.transpose()<<", factor: "<<factor.transpose()<<std::endl;
			}
				
			b[i]=new KX(counts);
			x[i]=new KX(counts);
			r[i]=new KX(counts);}

		data_initialized=true;
	}

	////solve Ax=b
	template<class KXW> void Precond(const KXW& _b,/*rst*/KXW& _x)
	{
		auto& x=gmg.x;
		auto& b=gmg.b;
		auto& A=gmg.A;
		(*b[0]).Copy(_b);
		(*x[0]).Copy(_x);
		gmg.V_Cycle();
		_x.Copy((*x[0]));

		//using namespace KrylovPoisson;
		//{
		//	KrylovPoissonX<d,T> r0((*x[0]).kx->counts);
		//	(*A[0]).Res((*b[0]),(*x[0]),r0);
		//	std::cout<<"res 0: "<<r0.Dot(r0)<<std::endl;		
		//}

		//for(int i=0;i<10;i++){
		//	gmg.V_Cycle();
		//	_x.Copy((*x[0]));

		//	{
		//		KrylovPoissonX<d,T> r0((*x[0]).kx->counts);
		//		(*A[0]).Res((*b[0]),(*x[0]),r0);
		//		std::cout<<"res 1: "<<r0.Dot(r0)<<std::endl;		
		//	}
		//}

		////_x.Copy(_b);

		//{
		//	KrylovPoissonX<d,T> r0((*x[0]).kx->counts);
		//	(*A[0]).Res((*b[0]),(*x[0]),r0);
		//	std::cout<<"res 1: "<<r0.Dot(r0)<<std::endl;		
		//}
	}
};
};

template<int d,class T=real> class PoissonMtxFreeOld
{Typedef_VectorDii(d);
public:
	MacGrid<d> mac_grid;
	BitField<d> psi_D;						////psi_D flags
	Field<T,d> psi_D_values;				////psi_D values
	BitFaceField<d> psi_N;					////psi_N flags
	FaceField<T,d> psi_N_values;			////psi_N values
	FaceField<T,d> alpha;					////alpha
	Field<T,d> rhs;							////rhs
	Field<T,d> p;							////solution

	//////////////////////////////////////////////////////////////////////////
	////initialization
	void Initialize(const VectorDi& cell_counts,const real dx,const VectorD& domain_min=VectorD::Zero())
	{
		mac_grid.Initialize(cell_counts,dx,domain_min);
		psi_D.Resize(mac_grid.grid.cell_counts,false);
		psi_D_values.Resize(mac_grid.grid.cell_counts,(T)0);
		psi_N.Resize(mac_grid.grid.cell_counts,false);
		psi_N_values.Resize(mac_grid.grid.cell_counts,(T)0);
		alpha.Resize(mac_grid.grid.cell_counts,(T)1);
		p.Resize(mac_grid.grid.cell_counts);p.Fill((T)0);
		rhs.Resize(mac_grid.grid.cell_counts);rhs.Fill((T)0);

		std::cout<<"Initialize "<<cell_counts.transpose()<<std::endl;
	}
	void Initialize(const Grid<d>& _grid){Initialize(_grid.cell_counts,_grid.dx,_grid.domain_min);}

	//////////////////////////////////////////////////////////////////////////
	////boundary conditions
	void Set_Psi_N(const int axis,const VectorDi& face,const real value){psi_N_values(axis,face)=value;psi_N.Set(axis,face);}
	void Set_Psi_D(const VectorDi& cell,const real value=(real)0){psi_D_values(cell)=value;psi_D.Set(cell);}

	//////////////////////////////////////////////////////////////////////////
	////matrix-free implementation

	void Update_Rhs_With_BC()
	{
		real dx=mac_grid.grid.dx;
		real dx2=(T)pow(mac_grid.grid.dx,2);

		////set b
		int n=(int)rhs.array.size();
		for(int i=0;i<n;i++){rhs.array[i]=-dx2*rhs.array[i];}

		////set Dirichlet boundary rhs
		int cn=mac_grid.grid.cell_counts.prod();
		for(int j=0;j<cn;j++){
			VectorDi cell=mac_grid.grid.Cell_Coord(j);
			if(!psi_D(cell))continue;
			real val=psi_D_values(cell);
			////rhs
			rhs(cell)=val;
			////nbs' rhs
			for(int i=0;i<Grid<d>::Number_Of_Nb_C();i++){VectorDi nb_cell=Grid<d>::Nb_C(cell,i);
				if(mac_grid.grid.Valid_Cell(nb_cell)&&!psi_D(nb_cell)){
					int axis=0;int side=0;mac_grid.grid.Nb_C_Axis_And_Side(i,axis,side);
					VectorDi face=cell+VectorDi::Unit(axis)*side;real a=alpha(axis,face);rhs(nb_cell)+=a*val;}}}

		////set Neumann boundary rhs
		for(int axis=0;axis<d;axis++){int face_num=mac_grid.Number_Of_Faces(axis);
			for(int i=0;i<face_num;i++){VectorDi face=mac_grid.Face_Coord(axis,i);
				if(!psi_N(axis,face))continue;
				real value=psi_N_values(axis,face);
				const VectorDi cell_0=MacGrid<d>::Face_Incident_Cell(axis,face,0);
				const VectorDi cell_1=MacGrid<d>::Face_Incident_Cell(axis,face,1);
				if(mac_grid.grid.Valid_Cell(cell_0)){
					rhs(cell_0)+=alpha(axis,face)*value*dx;}
				if(mac_grid.grid.Valid_Cell(cell_1)){
					rhs(cell_1)-=alpha(axis,face)*value*dx;}}}
	}

	//////////////////////////////////////////////////////////////////////////
	////solve
	void Build()
	{
		Update_Rhs_With_BC();
	}

	void Solve()
	{
		using namespace KrylovPoisson;
		KrylovPoissonA<d,real> kA(p.counts,&alpha,&psi_D,&psi_D_values,&psi_N,&psi_N_values);

		KrylovPoissonX<d,real> kx(&p),kb(&rhs);
		KrylovPoissonX<d,real> kr(p.counts),kd(p.counts),kq(p.counts),ks(p.counts);
		//KrylovPoissonPrecIdentity<d,real> kP;
		KrylovPoissonPrecGMG<d,real,KrylovPoissonA<d,real>,KrylovPoissonX<d,real> > kP(&kA);

		//kb.Print();
		KrylovSolver::Preconditioned_Conjugate_Gradient(kA,kx,kb,kr,kd,kq,ks,kP);
	}

	void Build_And_Solve()
	{
		Build();
		Solve();
	}
};

#endif
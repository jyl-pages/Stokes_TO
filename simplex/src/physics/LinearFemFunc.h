#ifndef __LinearFemFunc_h__
#define __LinearFemFunc_h__
#include "Common.h"
#include "Constants.h"
#include "SparseFunc.h"

namespace ElasticParamImpl 
{
    template<int dim> inline real Lambda(const real youngs,const real poisson){return youngs*poisson/(((real)1+poisson)*((real)1-(real)2*poisson));}
    template<> inline real Lambda<2>(const real youngs,const real poisson){return youngs*poisson/((real)1-pow(poisson,2));}

	template<int dim> inline real Bulk(const real youngs,const real poisson){return (real).5*youngs/((real)1-poisson);}
	template<> inline real Bulk<3>(const real youngs,const real poisson){return (real)one_third*youngs/((real)1-(real)2*poisson);}

	template<int dim>  inline real E_Hat_Coef(const real poisson){return (real)0;}	////E=E_hat*coef
	template<> inline real E_Hat_Coef<2>(const real poisson){return ((real)1-poisson*poisson);}
	template<> inline real E_Hat_Coef<3>(const real poisson){return ((real)1-2*poisson)*((real)1+poisson);}

	template<int dim> inline real G_Hat_Coef(const real poisson){return (real)0;}	////G=E_hat*coef
	template<> inline real G_Hat_Coef<2>(const real poisson){return (real).5*((real)1-poisson);}	
	template<> inline real G_Hat_Coef<3>(const real poisson){return (real).5*((real)1-(real)2*poisson);}
}

class ElasticParam
{public:
	real youngs_modulus=(real)1.;
	real poisson_ratio=(real).45;
	
	ElasticParam(const real _youngs,const real _poisson):youngs_modulus(_youngs),poisson_ratio(_poisson){}

	////Lame parameters
	template<int dim> real Lambda() const {return Lambda<dim>(youngs_modulus,poisson_ratio);}
	template<int dim> real Mu() const {return Mu<dim>(youngs_modulus,poisson_ratio);}
	
	template <int dim> static real Lambda(const real youngs,const real poisson) {return ElasticParamImpl::Lambda<dim>(youngs, poisson);}
	template<int dim> static real Mu(const real youngs,const real poisson){return (real).5*youngs/((real)1+poisson);}	////shear

	template<int dim> real Bulk() const {return Bulk<dim>(youngs_modulus,poisson_ratio);}
	template<int dim> static real Bulk(const real youngs,const real poisson){ElasticParamImpl::Bulk<dim>(youngs, poisson);}
	template<int dim> real Shear() const{return Shear<dim>(youngs_modulus,poisson_ratio);}
	template<int dim> static real Shear(const real youngs,const real poisson){return (real).5*youngs/((real)1+poisson);}

	////E_hat
	template<int dim> static real E_Hat_Coef(const real poisson){ElasticParamImpl::E_Hat_Coef<dim>(poisson);}
	////G_hat
	template<int dim> static real G_Hat_Coef(const real poisson){ElasticParamImpl::G_Hat_Coef<dim>(poisson);}
};

template<int d> class LinearFemFunc
{Typedef_VectorDii(d);Typedef_MatrixD(d);Typedef_VectorEi(d+1);
public:
	//////////////////////////////////////////////////////////////////////////
	////Material model
	static void Strain_Stress_Matrix_Linear(const real youngs,const real poisson,MatrixX& E);
	static void Strain_Stress_Matrix_Linear_From_Shear_And_Lame(const real shear,const real lame,MatrixX& E);

	//////////////////////////////////////////////////////////////////////////
	////Hex element
	static void dNde(const VectorD& natural_coord,MatrixX& dNde);					////dNde, return a dxn matrix, n is the number of basis functions
	static VectorD dNde(const VectorD& natural_coord,const int idx);				////dNde for the idx'th basis function, return a d-dim vector
	
	static void Cell_dNdX(const VectorD& natural_coord,const real dx,MatrixX& dNdX);		////dNdX
	static void Hex_dNdX(const VectorD& natural_coord,const MatrixD& J_invT,MatrixX& dNdX);

	static void Cell_dXde(const VectorD& natural_coord,const real dx,MatrixD& dXde);						////Jacobian
	static void Hex_dXde(const VectorD& natural_coord,const ArrayF2P<VectorD,d>& X,MatrixD& dXde);			

	static real Cell_Strain_Displacement_Matrix(const VectorD& natural_coord,const real dx,MatrixX& B);
	static real Hex_Strain_Displacement_Matrix(const VectorD& natural_coord,const ArrayF2P<VectorD,d>& X,MatrixX& B);

	static void Cell_Stiffness_Matrix(const real youngs,const real poisson,const real dx,MatrixX& K_e,const MatrixX* E0=nullptr);
	static void Hex_Stiffness_Matrix(const real youngs,const real poisson,const ArrayF2P<VectorD,d>& X,MatrixX& K_e,const MatrixX* E0=nullptr);

	////Hex helper functions
	static void Add_Cell_Stiffness_Matrix(SparseMatrixT& K,const MatrixX& K_e,const Array<int>& nodes);
	static void Add_Cell_Vector(VectorX& f,const VectorX& fi,const Array<int>& nodes);
	static void Set_Cell_B_Elements_Helper(const int r,const int c,const VectorX& dN,MatrixX& B,const real coef=(real)1);

	//////////////////////////////////////////////////////////////////////////
	////Tet element
	static void Tet_dNdX(const ArrayF<VectorD,d+1>& tet,ArrayF<ArrayF<real,d+1>,d>& dNdx);
	static void Tet_Strain_Displacement_Matrix(const ArrayF<VectorD,d+1>& tet,MatrixX& B);
	static void Tet_Stiffness_Matrix(const real youngs,const real poisson,const ArrayF<VectorD,d+1>& tet,MatrixX& K_e,const MatrixX* E0=nullptr);
	////Tet helper functions
	static void Add_Tet_Stiffness_Matrix(SparseMatrixT& K,const MatrixX& K_e,const VectorEi& nodes);
	static void Set_Tet_B_Elements_Helper(const int r,const int c,/*rst*/MatrixX& B,const ArrayF<real,d+1>& e,const real coef=(real)1);

	//////////////////////////////////////////////////////////////////////////
	////Beam element
	static void Add_Beam_Stiffness_Matrix(SparseMatrixT& K,const MatrixX& K_e,const Vector2i& nodes);

	//////////////////////////////////////////////////////////////////////////
	////Helper functions
	static void Set_Dirichlet_Boundary_Helper(SparseMatrixT& K,VectorX& b,const int i,const real psi_D_value);
	////Gaussian integration
	static void Initialize_Gaussian_Integration_Points(ArrayF2P<VectorD,d>& points,ArrayF2P<real,d>& weights);
	static void Initialize_Natural_Coord_Points(ArrayF2P<VectorD,d>& points,ArrayF2P<real,d>& weights);
	static VectorD Natural_Coord_Point(const int idx);

	static void Symmetric_Tensor_Vector_To_Matrix(const VectorX& v,MatrixD& m);
	static void Symmetric_Tensor_Matrix_To_Vector(const MatrixD& m,VectorX& v);
	static real Von_Mises_Stress(const Matrix<real,d>& e);

	//static void Cell_dXde(const VectorD& natural_coord,const ArrayF2P<VectorD,d>& X,MatrixD& dXde);

	//////Jacobian matrix dXde
	//void Compute_dXde(const VectorD& natural_coord,const real dx,MatrixD& dXde);
	//void Compute_dXde(const VectorD& natural_coord,const ArrayF2P<VectorD,d>& X,MatrixD& dXde);
	//real Build_Cell_Strain_Displacement_Matrix(const VectorD& natural_coord,const real dx,MatrixX& B,
	//	const bool use_jacobian,const ArrayF2P<VectorD,d>* X);
	//////hex stiffness matrix K
	//void Build_Cell_Stiffness_Matrix(const real youngs,const real poisson,const real dx,MatrixX& K_e,const MatrixX* E0=nullptr,
	//	const bool use_jacobian=false,const ArrayF2P<VectorD,d>* X=nullptr);

	////////////////////////////////////////////////////////////////////////////
	//////Tet Fem
	//////tet displacement-strain matrix B
	//void Compute_Tet_dNdx(const ArrayF<VectorD,d+1>& tet,ArrayF<ArrayF<real,d+1>,d>& dNdx);
	//void Build_Tet_Strain_Displacement_Matrix(const ArrayF<VectorD,d+1>& tet,MatrixX& B);
	//////tet stiffness matrix K
	//void Build_Tet_Stiffness_Matrix(const real youngs,const real poisson,const ArrayF<VectorD,d+1>& tet,MatrixX& K_e,const MatrixX* E0=nullptr);
	//////tet stiffness-thermal matrix Kh_e
	//void Build_Tet_Thermal_Traction_Matrix(const real youngs,const real poisson,const ArrayF<VectorD,d+1>& tet,const VectorD& alpha,MatrixX& Kh_e);
	//////derivative
	//void Compute_dJdX(const VectorD& natural_coord,const int idx,const int axis,MatrixD& dJdX);
	//real Compute_dJacdX(const MatrixD& dJdX,const MatrixD& J_inv,const real jac);	////d|J|/dX
	//void Build_Hex_dBdX(const VectorD& natural_coord,const int idx,const int axis,const ArrayF2P<VectorD,d>& X,MatrixX& B);
	//void Build_Hex_dKdX(const int idx,const int axis,const MatrixX& E0,const ArrayF2P<VectorD,d>& X,MatrixX& K_e);
	
	////////////////////////////////////////////////////////////////////////////
	//////Helper functions
	//////Gaussian integration points and weights
	//void Initialize_Gaussian_Integration_Points(ArrayF2P<VectorD,d>& points,ArrayF2P<real,d>& weights);
	//void Initialize_Natural_Coord_Points(ArrayF2P<VectorD,d>& points,ArrayF2P<real,d>& weights);
	//VectorD Natural_Coord_Point(const int idx);
	//void Symmetric_Tensor_Vector_To_Matrix(const VectorX& v,/*rst*/MatrixD& m);
	//void Symmetric_Tensor_Matrix_To_Vector(const MatrixD& m,/*rst*/VectorX& v);
	//real Von_Mises_Stress(const Matrix<real,d>& e);
};

#endif

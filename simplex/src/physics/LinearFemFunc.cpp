#include "Grid.h"
#include "AuxFunc.h"
#include "MeshFunc.h"
#include "LinearFemFunc.h"

//////////////////////////////////////////////////////////////////////////
////Material model
template<> void LinearFemFunc<2>::Strain_Stress_Matrix_Linear(const real youngs,const real poisson,MatrixX& E)
{
	////zero strain in z
	//E.resize(3,3);E.fill((real)0);real c=youngs/(((real)1-(real)2*poisson)*((real)1+poisson));real G=youngs/((real)2*((real)1+poisson));
	//E(0,0)=E(1,1)=c*((real)1-poisson);E(0,1)=E(1,0)=c*poisson;E(2,2)=G;

	////zero stress in z
	E.resize(3,3);E.fill((real)0);real e=youngs/((real)1-pow(poisson,2));real G=youngs/((real)2*((real)1+poisson));
	E(0,0)=E(1,1)=e;E(0,1)=E(1,0)=e*poisson;E(2,2)=G;
}

template<> void LinearFemFunc<3>::Strain_Stress_Matrix_Linear(const real youngs,const real poisson,MatrixX& E)
{
	E.resize(6,6);E.fill((real)0);real e=youngs/(((real)1-(real)2*poisson)*((real)1+poisson));real G=youngs/((real)2*((real)1+poisson));
	E(0,0)=E(1,1)=E(2,2)=e*((real)1-poisson);E(0,1)=E(1,0)=E(0,2)=E(2,0)=E(1,2)=E(2,1)=e*poisson;E(3,3)=E(4,4)=E(5,5)=G;
}

template<> void LinearFemFunc<2>::Strain_Stress_Matrix_Linear_From_Shear_And_Lame(const real shear,const real lame,MatrixX& E)
{
	E.resize(3,3);E.fill((real)0);
	E(0,0)=E(1,1)=(real)2*shear+lame;
	E(0,1)=E(1,0)=lame;
	E(2,2)=shear;
}

template<> void LinearFemFunc<3>::Strain_Stress_Matrix_Linear_From_Shear_And_Lame(const real shear,const real lame,MatrixX& E)
{
	E.resize(6,6);E.fill((real)0);
	E(0,0)=E(1,1)=E(2,2)=(real)2*shear+lame;
	E(0,1)=E(1,0)=E(0,2)=E(2,0)=E(1,2)=E(2,1)=lame;
	E(3,3)=E(4,4)=E(5,5)=shear;
}

//////////////////////////////////////////////////////////////////////////
////Hex element
template<> void LinearFemFunc<2>::dNde(const Vector2& natural_coord,MatrixX& dNde)
{
	const real s=natural_coord[0];const real t=natural_coord[1];
	dNde<<t-1,-t-1,-t+1,t+1,
		s-1,-s+1,-s-1,s+1;
	dNde*=(real).25;
}

template<> void LinearFemFunc<3>::dNde(const Vector3& natural_coord,MatrixX& dNde)
{
	const real s=natural_coord[0];const real t=natural_coord[1];const real q=natural_coord[2];
	dNde<<-(1-t)*(1-q),-(1-t)*(1+q),-(1+t)*(1-q),-(1+t)*(1+q),(1-t)*(1-q),(1-t)*(1+q),(1+t)*(1-q),(1+t)*(1+q),
		-(1-s)*(1-q),-(1-s)*(1+q),(1-s)*(1-q),(1-s)*(1+q),-(1+s)*(1-q),-(1+s)*(1+q),(1+s)*(1-q),(1+s)*(1+q),
		-(1-s)*(1-t),(1-s)*(1-t),-(1-s)*(1+t),(1-s)*(1+t),-(1+s)*(1-t),(1+s)*(1-t),-(1+s)*(1+t),(1+s)*(1+t);
	dNde*=(real).125;
}

template<> Vector2 LinearFemFunc<2>::dNde(const Vector2& natural_coord,const int idx)
{
	const real s=natural_coord[0];const real t=natural_coord[1];
	switch(idx){
	case 0:return (real).25*Vector2(t-1,s-1);
	case 1:return (real).25*Vector2(-t-1,-s+1);
	case 2:return (real).25*Vector2(-t+1,-s-1);
	case 3:return (real).25*Vector2(t+1,s+1);
	default:return Vector2::Zero();}
}

template<> Vector3 LinearFemFunc<3>::dNde(const Vector3& natural_coord,const int idx)
{
	const real s=natural_coord[0];const real t=natural_coord[1];const real q=natural_coord[2];
	switch(idx){
	case 0:return (real).125*Vector3(-(1-t)*(1-q),-(1-s)*(1-q),-(1-s)*(1-t));
	case 1:return (real).125*Vector3(-(1-t)*(1+q),-(1-s)*(1+q),(1-s)*(1-t));
	case 2:return (real).125*Vector3(-(1+t)*(1-q),(1-s)*(1-q),-(1-s)*(1+t));
	case 3:return (real).125*Vector3(-(1+t)*(1+q),(1-s)*(1+q),(1-s)*(1+t));
	case 4:return (real).125*Vector3((1-t)*(1-q),-(1+s)*(1-q),-(1+s)*(1-t));
	case 5:return (real).125*Vector3((1-t)*(1+q),-(1+s)*(1+q),(1+s)*(1-t));
	case 6:return (real).125*Vector3((1+t)*(1-q),(1+s)*(1-q),-(1+s)*(1+t));
	case 7:return (real).125*Vector3((1+t)*(1+q),(1+s)*(1+q),(1+s)*(1+t));
	default:return Vector3::Zero();}
}

template<int d> void LinearFemFunc<d>::Cell_dNdX(const VectorD& natural_coord,const real dx,MatrixX& dNdX)
{
	dNde(natural_coord,dNdX);dNdX*=((real).5/dx);
}

template<int d> void LinearFemFunc<d>::Hex_dNdX(const VectorD& natural_coord,const MatrixD& J_invT,MatrixX& dNdX)
{
	dNde(natural_coord,dNdX);
	dNdX=J_invT*dNdX;
}

template<int d> void LinearFemFunc<d>::Cell_dXde(const VectorD& natural_coord,const real dx,MatrixD& dXde)
{
	int n=Grid<d>::Number_Of_Cell_Incident_Nodes();dXde=MatrixD::Zero();VectorDi cell=VectorDi::Zero();
	for(int i=0;i<n;i++){
		VectorDi node=Grid<d>::Cell_Incident_Node(cell,i);	
		VectorD X=dx*node.template cast<real>();
		dXde+=X*dNde(natural_coord,i).transpose();}
}


template<int d> void LinearFemFunc<d>::Hex_dXde(const VectorD& natural_coord,const ArrayF2P<VectorD,d>& X,MatrixD& dXde)
{
	int n=Grid<d>::Number_Of_Cell_Incident_Nodes();dXde=MatrixD::Zero();VectorDi cell=VectorDi::Zero();
	for(int i=0;i<n;i++){
		VectorDi node=Grid<d>::Cell_Incident_Node(cell,i);	
		dXde+=X[i]*dNde(natural_coord,i).transpose();}
}

template<> real LinearFemFunc<2>::Cell_Strain_Displacement_Matrix(const Vector<real,2>& natural_coord,const real dx,MatrixX& B)
{
	const int d=2;int tensor_n=d*(d+1)/2;int r=tensor_n;int vtx_n=(int)pow(2,d);int c=d*vtx_n;
	B.resize(r,c);B.fill(0);
	MatrixX dNdX(d,vtx_n);Cell_dNdX(natural_coord,dx,dNdX);
	real J_det=pow((real).5*dx,d);
			
	const int x=0;const int y=1;
	Set_Cell_B_Elements_Helper(0,x,dNdX.row(0),B);
	Set_Cell_B_Elements_Helper(1,y,dNdX.row(1),B);
	Set_Cell_B_Elements_Helper(2,x,dNdX.row(1),B);
	Set_Cell_B_Elements_Helper(2,y,dNdX.row(0),B);
		
	return J_det;
}

template<> real LinearFemFunc<3>::Cell_Strain_Displacement_Matrix(const Vector<real,3>& natural_coord,const real dx,MatrixX& B)
{
	const int d=3;int tensor_n=d*(d+1)/2;int r=tensor_n;int vtx_n=(int)pow(2,d);int c=d*vtx_n;
	B.resize(r,c);B.fill(0);
	MatrixX dNdX(d,vtx_n);Cell_dNdX(natural_coord,dx,dNdX);
	real J_det=pow((real).5*dx,d);

	const int x=0;const int y=1;const int z=2;
	Set_Cell_B_Elements_Helper(0,x,dNdX.row(0),B);
	Set_Cell_B_Elements_Helper(1,y,dNdX.row(1),B);
	Set_Cell_B_Elements_Helper(2,z,dNdX.row(2),B);
	Set_Cell_B_Elements_Helper(3,x,dNdX.row(1),B);	
	Set_Cell_B_Elements_Helper(3,y,dNdX.row(0),B);	
	Set_Cell_B_Elements_Helper(4,y,dNdX.row(2),B);	
	Set_Cell_B_Elements_Helper(4,z,dNdX.row(1),B);	
	Set_Cell_B_Elements_Helper(5,x,dNdX.row(2),B);	
	Set_Cell_B_Elements_Helper(5,z,dNdX.row(0),B);	
	
	return J_det;
}

template<> real LinearFemFunc<2>::Hex_Strain_Displacement_Matrix(const Vector2& natural_coord,const ArrayF2P<Vector2,2>& X,MatrixX& B)
{
	const int d=2;int tensor_n=d*(d+1)/2;int r=tensor_n;int vtx_n=(int)pow(2,d);int c=d*vtx_n;
	B.resize(r,c);B.fill(0);
	
	real J_det=(real)0;
	Matrix2 J_invT;
	Hex_dXde(natural_coord,X,J_invT);
	J_det=J_invT.determinant();	////J_det=|J|
	J_invT=J_invT.inverse().transpose();
	MatrixX dNdX(d,vtx_n);Hex_dNdX(natural_coord,J_invT,dNdX);
			
	const int x=0;const int y=1;
	Set_Cell_B_Elements_Helper(0,x,dNdX.row(0),B);
	Set_Cell_B_Elements_Helper(1,y,dNdX.row(1),B);
	Set_Cell_B_Elements_Helper(2,x,dNdX.row(1),B);
	Set_Cell_B_Elements_Helper(2,y,dNdX.row(0),B);
		
	return J_det;
}

template<> real LinearFemFunc<3>::Hex_Strain_Displacement_Matrix(const Vector3& natural_coord,const ArrayF2P<Vector3,3>& X,MatrixX& B)
{
	const int d=3;int tensor_n=d*(d+1)/2;int r=tensor_n;int vtx_n=(int)pow(2,d);int c=d*vtx_n;
	B.resize(r,c);B.fill(0);
	
	real J_det=(real)0;
	MatrixD J_invT;
	Hex_dXde(natural_coord,X,J_invT);
	J_det=J_invT.determinant();	////J_det=|J|
	J_invT=J_invT.inverse().transpose();
	MatrixX dNdX(d,vtx_n);Hex_dNdX(natural_coord,J_invT,dNdX);

	const int x=0;const int y=1;const int z=2;
	Set_Cell_B_Elements_Helper(0,x,dNdX.row(0),B);
	Set_Cell_B_Elements_Helper(1,y,dNdX.row(1),B);
	Set_Cell_B_Elements_Helper(2,z,dNdX.row(2),B);
	Set_Cell_B_Elements_Helper(3,x,dNdX.row(1),B);	
	Set_Cell_B_Elements_Helper(3,y,dNdX.row(0),B);	
	Set_Cell_B_Elements_Helper(4,y,dNdX.row(2),B);	
	Set_Cell_B_Elements_Helper(4,z,dNdX.row(1),B);	
	Set_Cell_B_Elements_Helper(5,x,dNdX.row(2),B);	
	Set_Cell_B_Elements_Helper(5,z,dNdX.row(0),B);	
	
	return J_det;
}

template<int d> void LinearFemFunc<d>::Cell_Stiffness_Matrix(const real youngs,const real poisson,const real dx,MatrixX& K_e,const MatrixX* E0/*=nullptr*/)
{
	int n=d*(int)pow(2,d);K_e.resize(n,n);K_e.fill((real)0);
	ArrayF2P<VectorD,d> points;ArrayF2P<real,d> weights;Initialize_Gaussian_Integration_Points(points,weights);
	for(auto i=0;i<points.size();i++){
		MatrixX B;real J=Cell_Strain_Displacement_Matrix(points[i],dx,B);
		MatrixX K0;if(E0!=0){K0=B.transpose()*(*E0)*B*J*weights[i];}
		else{MatrixX E;Strain_Stress_Matrix_Linear(youngs,poisson,E);K0=B.transpose()*E*B*J*weights[i];}
		K_e+=K0;}
}

template<int d> void LinearFemFunc<d>::Hex_Stiffness_Matrix(const real youngs,const real poisson,const ArrayF2P<VectorD,d>& X,MatrixX& K_e,const MatrixX* E0/*=nullptr*/)
{
	int n=d*(int)pow(2,d);K_e.resize(n,n);K_e.fill((real)0);
	ArrayF2P<VectorD,d> points;ArrayF2P<real,d> weights;Initialize_Gaussian_Integration_Points(points,weights);
	for(auto i=0;i<points.size();i++){
		MatrixX B;real J=Hex_Strain_Displacement_Matrix(points[i],X,B);
		MatrixX K0;if(E0!=0){K0=B.transpose()*(*E0)*B*J*weights[i];}
		else{MatrixX E;Strain_Stress_Matrix_Linear(youngs,poisson,E);K0=B.transpose()*E*B*J*weights[i];}
		K_e+=K0;}
}

template<int d> void LinearFemFunc<d>::Add_Cell_Stiffness_Matrix(/*rst*/SparseMatrixT& K,const MatrixX& K_e,const Array<int>& nodes)
{const int node_n=(int)nodes.size();for(int Ke_i=0;Ke_i<node_n;Ke_i++){
	int K_i=nodes[Ke_i];for(int Ke_j=0;Ke_j<node_n;Ke_j++){int K_j=nodes[Ke_j];SparseFunc::Add_Block<d>(K,K_i,K_j,K_e,Ke_i,Ke_j);}}}


template<int d> void LinearFemFunc<d>::Add_Cell_Vector(/*rst*/VectorX& f,const VectorX& fi,const Array<int>& nodes)
{const int node_n=(int)nodes.size();for(int i=0;i<node_n;i++){int p=nodes[i];for(int j=0;j<d;j++)f[p*d+j]+=fi[i*d+j];}}

template<int d> void LinearFemFunc<d>::Set_Cell_B_Elements_Helper(const int r,const int c,const VectorX& dN,MatrixX& B,const real coef)
{for(int i=0;i<(int)dN.size();i++)B(r,c+i*d)=dN[i]*coef;}

//////////////////////////////////////////////////////////////////////////
////Tet element
template<> void LinearFemFunc<2>::Tet_dNdX(const ArrayF<Vector2,3>& tri,ArrayF<ArrayF<real,3>,2>& dNdx)
{
	real x1=tri[0].x();real y1=tri[0].y();real x2=tri[1].x();real y2=tri[1].y();real x3=tri[2].x();real y3=tri[2].y();
	dNdx[0][0]=y2-y3;dNdx[0][1]=y3-y1;dNdx[0][2]=y1-y2;
	dNdx[1][0]=x3-x2;dNdx[1][1]=x1-x3;dNdx[1][2]=x2-x1;
}

template<> void LinearFemFunc<3>::Tet_dNdX(const ArrayF<Vector3,4>& tet,ArrayF<ArrayF<real,4>,3>& dNdx)
{
	real x1=tet[0].x();real y1=tet[0].y();real z1=tet[0].z();
	real x2=tet[1].x();real y2=tet[1].y();real z2=tet[1].z();
	real x3=tet[2].x();real y3=tet[2].y();real z3=tet[2].z();
	real x4=tet[3].x();real y4=tet[3].y();real z4=tet[3].z();
	real x12=x1-x2;real x13=x1-x3;real x14=x1-x4;real x23=x2-x3;real x24=x2-x4;real x34=x3-x4;
	real x21=-x12;real x31=-x13;real x41=-x14;real x32=-x23;real x42=-x24;real x43=-x34; 
	real y12=y1-y2;real y13=y1-y3;real y14=y1-y4;real y23=y2-y3;real y24=y2-y4;real y34=y3-y4;
	real y21=-y12;real y31=-y13;real y41=-y14;real y32=-y23;real y42=-y24;real y43=-y34; 
	real z12=z1-z2;real z13=z1-z3;real z14=z1-z4;real z23=z2-z3;real z24=z2-z4;real z34=z3-z4;
	real z21=-z12;real z31=-z13;real z41=-z14;real z32=-z23;real z42=-z24;real z43=-z34;
	dNdx[0][0]=y42*z32-y32*z42;dNdx[1][0]=x32*z42-x42*z32;dNdx[2][0]=x42*y32-x32*y42;
	dNdx[0][1]=y31*z43-y34*z13;dNdx[1][1]=x43*z31-x13*z34;dNdx[2][1]=x31*y43-x34*y13;
	dNdx[0][2]=y24*z14-y14*z24;dNdx[1][2]=x14*z24-x24*z14;dNdx[2][2]=x24*y14-x14*y24;
	dNdx[0][3]=y13*z21-y12*z31;dNdx[1][3]=x21*z13-x31*z12;dNdx[2][3]=x13*y21-x12*y31;		
}

template<> void LinearFemFunc<2>::Tet_Strain_Displacement_Matrix(const ArrayF<Vector2,3>& tri,MatrixX& B)
{	
	const int d=2;int tensor_n=d*(d+1)/2;int row=tensor_n;int vtx_n=d+1;int col=d*vtx_n;
	B.resize(row,col);B.setZero();	////3*6 for 2D
	ArrayF<ArrayF<real,3>,2> dNdx;Tet_dNdX(tri,dNdx);
	real coef=(real)1/((real)2*MeshFunc::Triangle_Area(tri));
	Set_Tet_B_Elements_Helper(0,0,B,dNdx[0],coef);
	Set_Tet_B_Elements_Helper(1,1,B,dNdx[1],coef);
	Set_Tet_B_Elements_Helper(2,0,B,dNdx[1],coef);
	Set_Tet_B_Elements_Helper(2,1,B,dNdx[0],coef);
}

template<> void LinearFemFunc<3>::Tet_Strain_Displacement_Matrix(const ArrayF<Vector3,4>& tet,MatrixX& B)
{	
	const int d=3;int tensor_n=d*(d+1)/2;int row=tensor_n;int vtx_n=d+1;int col=d*vtx_n;
	B.resize(row,col);B.setZero();	////6*12 for 3D
	ArrayF<ArrayF<real,4>,3> dNdx;Tet_dNdX(tet,dNdx);
	real coef=(real)1/((real)6*MeshFunc::Simplex_Size(tet));	//Tetrahedron_Volume
	Set_Tet_B_Elements_Helper(0,0,B,dNdx[0],coef);
	Set_Tet_B_Elements_Helper(1,1,B,dNdx[1],coef);
	Set_Tet_B_Elements_Helper(2,2,B,dNdx[2],coef);
	Set_Tet_B_Elements_Helper(3,0,B,dNdx[1],coef);
	Set_Tet_B_Elements_Helper(3,1,B,dNdx[0],coef);
	Set_Tet_B_Elements_Helper(4,1,B,dNdx[2],coef);
	Set_Tet_B_Elements_Helper(4,2,B,dNdx[1],coef);
	Set_Tet_B_Elements_Helper(5,0,B,dNdx[2],coef);
	Set_Tet_B_Elements_Helper(5,2,B,dNdx[0],coef);
}

template<int d> void LinearFemFunc<d>::Tet_Stiffness_Matrix(const real youngs,const real poisson,const ArrayF<VectorD,d+1>& tet,MatrixX& K_e,const MatrixX* E0)
{
	MatrixX B;Tet_Strain_Displacement_Matrix(tet,B);real vol=MeshFunc::Simplex_Size(tet);
	if(E0!=nullptr){K_e=vol*B.transpose()*(*E0)*B;}
	else{MatrixX E;Strain_Stress_Matrix_Linear(youngs,poisson,E);K_e=vol*B.transpose()*E*B;}
}

template<int d> void LinearFemFunc<d>::Add_Tet_Stiffness_Matrix(SparseMatrixT& K,const MatrixX& K_e,const VectorEi& nodes)
{
	const int node_n=(int)nodes.size();
	for(int Ke_i=0;Ke_i<node_n;Ke_i++){int K_i=nodes[Ke_i];
		for(int Ke_j=0;Ke_j<node_n;Ke_j++){int K_j=nodes[Ke_j];
			SparseFunc::Add_Block<d>(K,K_i,K_j,K_e,Ke_i,Ke_j);}}
}

template<int d> void LinearFemFunc<d>::Set_Tet_B_Elements_Helper(const int r,const int c,MatrixX& B,const ArrayF<real,d+1>& e,const real coef/*=(real)1*/)
{for(int i=0;i<(int)e.size();i++)B(r,c+i*d)=e[i]*coef;}


template<int d> void LinearFemFunc<d>::Add_Beam_Stiffness_Matrix(SparseMatrixT& K,const MatrixX& K_e,const Vector2i& nodes)
{
	constexpr int dof_per_vtx=d*(d+1)/2;
	for(int Ke_i=0;Ke_i<2;Ke_i++){int K_i=nodes[Ke_i];
		for(int Ke_j=0;Ke_j<2;Ke_j++){int K_j=nodes[Ke_j];
			SparseFunc::Add_Block<dof_per_vtx>(K,K_i,K_j,K_e,Ke_i,Ke_j);}}
}

template<int d> void LinearFemFunc<d>::Set_Dirichlet_Boundary_Helper(SparseMatrixT& K,VectorX& b,const int i,const real psi_D_value)
{
	for(InnerIteratorT iter(K,i);iter;++iter){int j=(int)iter.col();
		if(i==j){b(j)=psi_D_value*K.coeff(i,j);}
		else{real K_ij=K.coeff(i,j);b(j)-=K_ij*psi_D_value;
			K.coeffRef(i,j)=(real)0;K.coeffRef(j,i)=(real)0;}}
}

//////////////////////////////////////////////////////////////////////////
////Gaussian integration
template<> void LinearFemFunc<2>::Initialize_Gaussian_Integration_Points(ArrayF2P<Vector2,2>& points,ArrayF2P<real,2>& weights)
{real c=one_over_sqrt_three;points[0]=Vector2(-c,-c);points[1]=Vector2(-c,c);points[2]=Vector2(c,-c);points[3]=Vector2(c,c);AuxFunc::Fill(weights,(real)1);}

template<> void LinearFemFunc<3>::Initialize_Gaussian_Integration_Points(ArrayF2P<Vector3,3>& points,ArrayF2P<real,3>& weights)
{real c=one_over_sqrt_three;points[0]=Vector3(-c,-c,-c);points[1]=Vector3(-c,-c,c);points[2]=Vector3(-c,c,-c);points[3]=Vector3(-c,c,c);
points[4]=Vector3(c,-c,-c);points[5]=Vector3(c,-c,c);points[6]=Vector3(c,c,-c);points[7]=Vector3(c,c,c);AuxFunc::Fill(weights,(real)1);}	

template<> void LinearFemFunc<2>::Initialize_Natural_Coord_Points(ArrayF2P<Vector2,2>& points,ArrayF2P<real,2>& weights)
{real c=(real)1;points[0]=Vector2(-c,-c);points[1]=Vector2(-c,c);points[2]=Vector2(c,-c);points[3]=Vector2(c,c);AuxFunc::Fill(weights,(real)1);}

template<> void LinearFemFunc<3>::Initialize_Natural_Coord_Points(ArrayF2P<Vector3,3>& points,ArrayF2P<real,3>& weights)
{real c=(real)1;points[0]=Vector3(-c,-c,-c);points[1]=Vector3(-c,-c,c);points[2]=Vector3(-c,c,-c);points[3]=Vector3(-c,c,c);
points[4]=Vector3(c,-c,-c);points[5]=Vector3(c,-c,c);points[6]=Vector3(c,c,-c);points[7]=Vector3(c,c,c);AuxFunc::Fill(weights,(real)1);}

template<> Vector2 LinearFemFunc<2>::Natural_Coord_Point(const int idx)
{static const real c=(real)1;static Vector2 points[4]={{-c,-c},{-c,c},{c,-c},{c,c}};return points[idx];}

template<> Vector3 LinearFemFunc<3>::Natural_Coord_Point(const int idx)
{static const real c=(real)1;static Vector3 points[8]={{-c,-c,-c},{-c,-c,c},{-c,c,-c},{-c,c,c},{c,-c,-c},{c,-c,c},{c,c,-c},{c,c,c}};return points[idx];}

template<> void LinearFemFunc<2>::Symmetric_Tensor_Vector_To_Matrix(const VectorX& v,Matrix2& m)
{assert(v.size()==3);m<<v[0],v[2],v[2],v[1];}

template<> void LinearFemFunc<3>::Symmetric_Tensor_Vector_To_Matrix(const VectorX& v,Matrix3& m)
{assert(v.size()==6);m<<v[0],v[3],v[5],v[3],v[1],v[4],v[5],v[4],v[2];}

template<> void LinearFemFunc<2>::Symmetric_Tensor_Matrix_To_Vector(const Matrix2& m,VectorX& v)
{v.resize(3);v<<m(0,0),m(1,1),m(0,1);}

template<> void LinearFemFunc<3>::Symmetric_Tensor_Matrix_To_Vector(const Matrix3& m,VectorX& v)
{v.resize(6);v<<m(0,0),m(1,1),m(2,2),m(0,1),m(1,2),m(0,2);}

template<> real LinearFemFunc<2>::Von_Mises_Stress(const Matrix2& e)
{return pow(e(0,0),2)+pow(e(1,1),2)-e(0,0)*e(1,1)+(real)3*pow(e(0,1),2);}
template<> real LinearFemFunc<3>::Von_Mises_Stress(const Matrix3& e)
{return sqrt((real).5*(pow(e(0,0)-e(1,1),2)+pow(e(1,1)-e(2,2),2)+pow(e(2,2)-e(0,0),2))+(real)3*(pow(e(0,1),2)+pow(e(0,2),2)+pow(e(1,2),2)));}

template class LinearFemFunc<2>;
template class LinearFemFunc<3>;
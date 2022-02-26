#include "NonlinearFemFunc.h"
#include "LinearFemFunc.h"

template<int d> void NonlinearFemFunc<d>::Deformation_Gradient(const real dx,const VectorX& cell_u,const VectorD& natural_coord,/*rst*/MatrixD& F)
{
	F=MatrixD::Identity();int number_of_cell_nodes=(int)pow(2,d);
	MatrixD J_invT;LinearFemFunc<d>::Cell_dXde(natural_coord,dx,J_invT);
	J_invT=J_invT.inverse().transpose();
	for(int i=0;i<number_of_cell_nodes;i++){
		VectorD u_i;for(int j=0;j<d;j++)u_i[j]=cell_u[d*i+j];
		VectorD dNdX=J_invT*LinearFemFunc<d>::dNde(natural_coord,i);
		F+=u_i*dNdX.transpose();}
}

template<int d> void NonlinearFemFunc<d>::dPdx_Neohookean(const real mu,const real lambda,const VectorD& natural_coord,const MatrixD& F,const VectorD& dNidX_T,/*rst*/Array<MatrixD>& dPdx)
{
	dPdx.resize(d);real J=F.determinant();
	MatrixD F_invT=F.inverse().transpose();
	for(int i=0;i<d;i++){
		MatrixD dF=MatrixD::Zero();dF.row(i)=dNidX_T.transpose();
		dPdx[i]=mu*dF+(mu-lambda*log(J))*F_invT*dF.transpose()*F_invT+lambda*(F_invT.transpose()*dF).trace()*F_invT;}
}

template<int d> void NonlinearFemFunc<d>::Cell_Stiffness_Matrix_And_f_Nonlinear(const real mu,const real lambda,const real dx,const VectorX& cell_u,/*rst*/MatrixX& K_e,/*rst*/VectorX* f_elas/*=0*/,/*rst*/real* e_elas/*=0*/)
{
	int number_of_cell_nodes=(int)pow(2,d);int n=d*number_of_cell_nodes;K_e.resize(n,n);K_e.fill((real)0);
	ArrayF2P<VectorD,d> points;ArrayF2P<real,d> weights;LinearFemFunc<d>::Initialize_Gaussian_Integration_Points(points,weights);
	real J=(real)pow((real).5*dx,d);
	if(e_elas){*e_elas=(real)0;}

	for(int p=0;p<(int)points.size();p++){
		MatrixD J_invT;LinearFemFunc<d>::Cell_dXde(points[p],dx,J_invT);
		J_invT=J_invT.inverse().transpose();
		MatrixD F;Deformation_Gradient(dx,cell_u,points[p],F);
		MatrixD P;if(f_elas){First_PK_Stress_Neohookean(mu,lambda,F,P);}
		////update e_elas
		if(e_elas){*e_elas+=Elastic_Energy_Neohookean(mu,lambda,F)*J;}

		for(int i=0;i<number_of_cell_nodes;i++){
			VectorD dNidX_T=J_invT*LinearFemFunc<d>::dNde(points[p],i);
			if(f_elas){VectorD fi=P*dNidX_T*J;for(int j=0;j<d;j++)(*f_elas)[i*d+j]+=fi[j];}

			////update K
			for(int j=i;j<number_of_cell_nodes;j++){
				////compute dP_i/dx_j
				VectorD dNjdX_T=J_invT*LinearFemFunc<d>::dNde(points[p],j);
				Array<MatrixD> dPdx;dPdx_Neohookean(mu,lambda,points[p],F,dNjdX_T,dPdx);
				////assemble K_ij
				MatrixD K_ij=MatrixD::Zero();
				for(int k=0;k<d;k++){K_ij.col(k)=dPdx[k]*dNidX_T*J*weights[p];}
				////copy to K
				for(int ii=0;ii<d;ii++)for(int jj=0;jj<d;jj++){K_e(i*d+ii,j*d+jj)+=K_ij(ii,jj);}
				if(i!=j)for(int ii=0;ii<d;ii++)for(int jj=0;jj<d;jj++){K_e(j*d+jj,i*d+ii)=K_e(i*d+ii,j*d+jj);}}}}	
}


template class NonlinearFemFunc<2>;
template class NonlinearFemFunc<3>;

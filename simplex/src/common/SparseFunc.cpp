//////////////////////////////////////////////////////////////////////////
// Sparse solver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include <iostream>
#include "SparseFunc.h"

namespace SparseSolver{
    using ConjugateGradient=Eigen::ConjugateGradient<SparseMatrixT,Eigen::Upper,Eigen::IdentityPreconditioner>;
    using ICPCG=Eigen::ConjugateGradient<SparseMatrixT,Eigen::Upper,Eigen::IncompleteCholesky<real,Eigen::Upper,Eigen::NaturalOrdering<int> > >;
    using BiCGSTAB=Eigen::BiCGSTAB<SparseMatrixT>;
    using SparseLU=Eigen::SparseLU<SparseMatrixT>;

    bool Conjugate_Gradient(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params)
    {
        ConjugateGradient cg;
        cg.setMaxIterations(params.max_iter_num);
        cg.setTolerance(params.tolerance);
        cg.compute(A);
        if(cg.info()!=Eigen::Success){std::cerr<<"ERROR: [Sparse] Eigen CG solver factorization failed."<<std::endl;return false;}
        x = cg.solve(b);
        if(cg.info()!=Eigen::Success){std::cerr<<"ERROR: [Sparse] Eigen CG solver failed."<<std::endl;return false;}
        std::cout<<"EigCg solver converge: "<<cg.iterations()<<", errorerror: "<<cg.error()<<std::endl;return true;
    }

	bool IC_PCG(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params)
    {
        ICPCG cg;
        cg.setMaxIterations(params.max_iter_num);
        cg.setTolerance(params.tolerance);
        cg.compute(A);
        if(cg.info()!=Eigen::Success){std::cerr<<"ERROR: [Sparse] Eigen MICCG solver factorization failed."<<std::endl;return false;}
        x = cg.solve(b);
        if(cg.info()!=Eigen::Success){std::cerr<<"ERROR: [Sparse] Eigen MICCG solver failed."<<std::endl;return false;}
        std::cout<<"EigCg solver converge: "<<cg.iterations()<<", error: "<<cg.error()<<std::endl;return true;
    }

    bool BiCG_STAB(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const Params params)
    {
        BiCGSTAB bicg;
        bicg.setMaxIterations(params.max_iter_num);
        bicg.setTolerance(params.tolerance);
        bicg.compute(A);
        if(bicg.info()!=Eigen::Success){std::cerr<<"ERROR: [Sparse] Eigen BiCGSTAB solver factorization failed."<<std::endl;return false;}
        x = bicg.solve(b);
        if(bicg.info()!=Eigen::Success){std::cerr<<"ERROR: [Sparse] EigBicg solver failed."<<std::endl;return false;}
        std::cout<<"Eigen BiCGSTAB solver converge: "<<bicg.iterations()<<", error: "<<bicg.error()<<std::endl;return true;
    }

    template<> bool LU<float>(const SparseMatrix<float>& A,VectorN<float>& x,const VectorN<float>& b,const Params)
    {MatrixXf Ad(A);x=Ad.fullPivLu().solve(b);return true;}

	template<> bool LU<double>(const SparseMatrix<double>& A,VectorN<double>& x,const VectorN<double>& b,const Params)
    {MatrixXd Ad(A);x=Ad.fullPivLu().solve(b);return true;}

	template<> bool Least_Squares<float>(const SparseMatrix<float>& A, VectorN<float>& x, const VectorN<float>& b, const Params params/*=Params()*/)
	{MatrixXf Ad(A);x=Ad.colPivHouseholderQr().solve(b);return true;}

	template<> bool Least_Squares<double>(const SparseMatrix<double>& A, VectorN<double>& x, const VectorN<double>& b, const Params params/*=Params()*/)
	{MatrixXd Ad(A);x=Ad.colPivHouseholderQr().solve(b);return true;}
};
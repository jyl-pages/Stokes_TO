//////////////////////////////////////////////////////////////////////////
// AMGCL Solver
// Copyright (c) (2018-), Bo Zhu
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __AmgclSolver_h__
#define __AmgclSolver_h__

#include <numeric>
#include <functional>
#include "SparseFunc.h"
#include "GeometricMultiGrid.h"
#include "Timer.h"

//// include AMGCL headers
#ifndef USE_CUDA
#include "assert.h"
#include <amgcl/amg.hpp>
#include <amgcl/make_solver.hpp>
#include <amgcl/adapter/eigen.hpp>
#include <amgcl/backend/builtin.hpp>
#include <amgcl/adapter/crs_builder.hpp>
#include <amgcl/coarsening/smoothed_aggregation.hpp>
#include <amgcl/relaxation/spai0.hpp>
#include <amgcl/solver/bicgstab.hpp>
#include <amgcl/profiler.hpp>

//// see https://github.com/ddemidov/amgcl/issues/103
//// to enable Eigen matrix with different backend
AMGCL_USE_EIGEN_VECTORS_WITH_BUILTIN_BACKEND()
#endif

namespace AmgclSolver{

template<int d> bool GMGPCG_AMGCL(const SparseMatrix<real>& A,VectorN<real>& x,const VectorN<real>& b,const GeometricMultiGrid::Params params/*=Params()*/)
{
#ifdef USE_CUDA
	std::cerr<<"Error: [GMGPCG_AMGCL] USE_AMGCL disabled"<<std::endl;
	return false;
#else
    Timer<real> timer;timer.Reset();
    // Setup the solver:
    typedef amgcl::make_solver<
        amgcl::amg<
            amgcl::backend::builtin<real>,
            amgcl::coarsening::smoothed_aggregation,
            amgcl::relaxation::spai0
            >,
        amgcl::solver::bicgstab<amgcl::backend::builtin<real> >
        > Solver;

    timer.Reset();
    Solver solve(A);
    timer.Elapse_And_Output_And_Reset("AMGCL Allocate A");
    std::cout << solve << std::endl;

    // Solve the system for the given RHS:
    int    iters;
    real error;
    timer.Reset();
    std::tie(iters, error) = solve(b, x);
    timer.Elapse_And_Output_And_Reset("AMGCL Solve");

    std::cout <<"#Iters="<< iters << ", error=" << error << std::endl;
    return 0;

#endif
}

};
#endif
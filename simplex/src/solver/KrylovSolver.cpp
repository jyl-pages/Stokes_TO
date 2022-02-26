//////////////////////////////////////////////////////////////////////////
// Krylov solver
// Copyright (c) (2018-), Bo Zhu, Jinyuan Liu, Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#include "KrylovSolver.h"

namespace KrylovSolver {
	inline void rotmat(const real& x, const real& y, real* c, real* s)
	{
		if (x == 0.)
		{
			*c = 0.;
			*s = 1.;
		}
		else if (fabs(y) > fabs(x))
		{
			real tmp = x / y;
			*s = 1.0 / sqrt(1.0 + tmp * tmp);
			*c = tmp * (*s);
		}
		else
		{
			real tmp = y / x;
			*c = 1.0 / sqrt(1.0 + tmp * tmp);
			*s = tmp * (*c);
		}
	}

	bool GMRES(const SparseMatrix<real>& A, VectorN<real>& x, const VectorN<real>& b, const Params params)
	{
		const int N = (int)b.size();
		int M = params.restart_iter_num;
		real b_n = b.norm();
		real err = (real)1;

		MatrixX V(MatrixX::Zero(N, M + 1));	// orthogonal vectors for Arnoldi iteration
		MatrixX H(MatrixX::Zero(M + 1, M));	// Hessenberg matrix for Arnoldi iteration

		ArrayX C(ArrayX::Zero(M + 1));	// cos of given rotations
		ArrayX S(ArrayX::Zero(M + 1));	// sin of given rotations

		for (int iter = 0; iter < params.max_iter_num; iter++)
		{
			// compute residual
			VectorX r = b - A * x;
			real r_n = r.norm();
			err = r_n / b_n;
			// return if already below tolerance
			if (err < params.tolerance) break;

			V.col(0) = r / r_n;	// obtain V_1
			VectorX g(VectorX::Zero(M + 1));	// vectorg = ||r|| * e1
			g(0) = r_n;

			// arnoldi iteration
			for (int i = 0; i < M; i++)
			{
				// form V_{i+1} and H(:,i)
				V.col(i + 1) = A * V.col(i);
				for (int j = 0; j <= i; j++)
				{
					H(j, i) = V.col(i + 1).dot(V.col(j));
					V.col(i + 1) -= H(j, i) * V.col(j);
				}
				H(i + 1, i) = V.col(i + 1).norm();
				if (H(i + 1, i) != 0) V.col(i + 1) /= H(i + 1, i);

				/*
				givens rotation
				| c, s| * |x| = |r|
				|-s, c|   |y|   |0|
				*/
				for (int j = 0; j < i; j++)
				{
					real tmp = C(j) * H(j, i) + S(j) * H(j + 1, i);
					H(j + 1, i) = -S(j) * H(j, i) + C(j) * H(j + 1, i);
					H(j, i) = tmp;
				}

				// rotate the last column
				rotmat(H(i, i), H(i + 1, i), &C(i), &S(i));
				H(i, i) = C(i) * H(i, i) + S(i) * H(i + 1, i);
				H(i + 1, i) = 0.;

				/* rotate g(i) and g(i+1)
				| c, s| * |g_i| = | c * g_i|
				|-s, c|   | 0 |   |-s * g_i|
				*/
				g(i + 1) = -S(i) * g(i);
				g(i) = C(i) * g(i);

				/*
				After all:
				H: [i+2, i+1]
				V: [N,   i+2]
				g: [i+2,   1]
				*/
				err = fabs(g(i + 1)) / b_n;
			}
			VectorX y = H.topLeftCorner(M, M).lu().solve(g.head(M));
			x += V.leftCols(M) * y;
			if (err < params.tolerance) break;
		}

		bool converged = err < params.tolerance;
		if (params.verbose) {
			std::cout << "KrylovSolve::GMRES: " << (converged ? "converge" : "fail")
				<< " with residual " << err << "." << std::endl;
		}
		return converged;
	}
}
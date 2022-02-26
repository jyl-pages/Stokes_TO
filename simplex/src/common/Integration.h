#pragma once
//////////////////////////////////////////////////////////////////////////
// Numerical Integration Method
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Common.h"
#include "AuxFunc.h"

namespace Integration {
	//https://en.wikipedia.org/wiki/Newton%E2%80%93Cotes_formulas
	//integrate f in [va,va+vh] with only one step
	template<int d, class T>
	T Newton_Cotes_Closed_Step(const Vector<real, d>& va, const Vector<real, d>& vh, std::function<T(const Vector<real, d>&)> f, const int degree) {
		real h = vh.norm();
		if (degree == 1) {//Trapezoidal rule
			T f0 = f(va), f1 = f(va + vh);
			return h / 2.0 * (f0 + f1);
		}
		else if (degree == 2) {//Simpson's rule
			T f0 = f(va), f1 = f(va + 0.5 * vh), f2 = f(va + vh);
			return h / 6.0 * (f0 + 4 * f1 + f2);
		}
		else if (degree == 3) {//Simpson's 3/8 rule
			T f0 = f(va), f1 = f(va + 1.0 / 3 * vh), f2 = f(va + 2.0 / 3 * vh), f3 = f(va + vh);
			return h / 8.0 * (f0 + 3 * f1 + 3 * f2 + f3);
		}
		else if (degree == 4) {
			T f0 = f(va), f1 = f(va + 0.25 * vh), f2 = f(va + 0.5 * vh), f3 = f(va + 0.75 * vh), f4 = f(va + vh);
			return h / 90.0 * (7 * f0 + 32 * f1 + 12 * f2 + 32 * f3 + 7 * f4);
		}
		else {
			AuxFunc::Crash_With_Info("Newton_Cotes: unsupported order");
			return Zero<T>();
		}
	}

	//Integrate f alone segment a-b. Splitting the segment into n pieces.
	template<int d, class T>
	T Newton_Cotes_Line(Vector<real, d> a, Vector<real, d> b, std::function<T(const Vector<real, d>&)> f, int n, const int degree = 2) {
		Vector<real, d> vh = (b - a) / (real)n;
		T sum = Zero<T>();
		for (int i = 0; i < n; i++) {
			sum += Newton_Cotes_Closed_Step<d, T>(a, vh, f, degree);
			a += vh;
		}
		return sum;
	}

	template<int d, class T>
	T Newton_Cotes(const Vector<real, d>& min_crn, const Vector<real, d> max_crn, std::function<T(const Vector<real, d>&)> f, const Vector<int, d>& res, const int degree = 2, int from_axis = 0) {
		std::function<T(const Vector1&)> int_f = nullptr;
		if (from_axis == d - 1) {
			int_f = [&](const Vector1& x)->real {Vector<real, d> pos = min_crn; pos[from_axis] = x[0]; return f(pos); };
		}
		else {
			int_f = [&](const Vector1& x)->real {
				Vector<real, d> _min = min_crn, _max = max_crn;
				_min[from_axis] = _max[from_axis] = x[0];
				return Newton_Cotes<d, T>(_min, _max, f, res, degree, from_axis + 1);
			};
		}
		Vector1 a(min_crn[from_axis]), b(max_crn[from_axis]);
		return Newton_Cotes_Line<1, T>(a, b, int_f, res[from_axis], degree);
	}
}

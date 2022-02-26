#include "RandomNumber.h"
#include <cmath>

namespace RandomFunc {
	RandomNumber unit_generator(0, 1);

	real Random_Real(real a, real b)
	{
		return unit_generator.Value() * (b - a) + a;
	}

	bool Sample_With_Probability(real prob)
	{
		real x = unit_generator.Value();
		return x <= prob;
	}

	int Stochastic_Coarsed_Real(real val)
	{
		int a = (int)floor(val);
		int b = a + 1;
		bool flg = Sample_With_Probability((real)b - val);
		return flg ? a : b;
	}

	template<int d>
	Vector<real, d> Random_Vector_Cartesian(real a, real b)
	{
		Vector<real, d> x;
		for (int i = 0; i < d; i++) x[i] = Random_Real(a, b);
		return x;
	}
	//note: these are just definitions
	template Vector<real, 1> Random_Vector_Cartesian<1>(real a, real b);
	template Vector<real, 2> Random_Vector_Cartesian<2>(real a, real b);
	template Vector<real, 3> Random_Vector_Cartesian<3>(real a, real b);

	//these are full implementations, so there's an empty <>
	template<> Vector<real, 1> Random_Vector_Spherical<1>(real max_r) {
		real r = Random_Real(0, max_r);
		bool flg = Sample_With_Probability(0.5);
		return flg ? Vector1(r) : Vector1(-r);
	}
	template<> Vector<real, 2> Random_Vector_Spherical<2>(real max_r) {
		real r = Random_Real(0, max_r);
		real theta = Random_Real(0, 2 * pi);
		return Vector2(r * cos(theta), r * sin(theta));
	}
	template<> Vector<real, 3> Random_Vector_Spherical<3>(real max_r) {
		real r = Random_Real(0, max_r);
		real theta = Random_Real(0, 2 * pi); 
		real phi = Random_Real(0, 2 * pi);
		return Vector3(r * sin(theta) * cos(phi), r * sin(theta) * sin(phi), r * cos(theta));
	}
}

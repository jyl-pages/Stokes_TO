#include "Common.h"

Vector1 Vec1(const real s) { Vector1 v; v[0] = s; return v; }
Vector1i Vec1i(const int s) { Vector1i v; v[0] = s; return v; }
Vector1d Vec1d(const double s) { Vector1d v; v[0] = s; return v; }
Vector1f Vec1f(const float s) { Vector1f v; v[0] = s; return v; }
Vector1c Vec1c(const C s) { Vector1c v; v[0] = s; return v; }

#define Define_Scalar_Valid(T) template<> bool Is_Valid_Number(const T &a){return (!std::isnan(a)) && std::isfinite(a);}
Define_Scalar_Valid(float); Define_Scalar_Valid(double); Define_Scalar_Valid(long double);
#undef Define_Scalar_Valid

#define Define_Vector_Valid(d) template<> bool Is_Valid_Number(const Vector<real,d> &a){for(int i=0;i<d;i++){if(!Is_Valid_Number(a[i])) return false;}return true;}
Define_Vector_Valid(1); Define_Vector_Valid(2); Define_Vector_Valid(3);
#undef Define_Vector_Valid

void Print(const std::string& msg)
{
	std::cout << msg << "\n";
}

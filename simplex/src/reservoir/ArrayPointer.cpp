//////////////////////////////////////////////////////////////////////////
// A class that can point to any Array<T>, to simplify class Particles<d>
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "ArrayPointer.h"

template class ArrayPointerDerived<int>;
template class ArrayPointerDerived<short>;
template class ArrayPointerDerived<bool>;
template class ArrayPointerDerived<float>;
template class ArrayPointerDerived<double>;
template class ArrayPointerDerived<C>;

template class ArrayPointerDerived<Vector1>;
template class ArrayPointerDerived<Vector2>;
template class ArrayPointerDerived<Vector3>;
template class ArrayPointerDerived<Vector4>;
template class ArrayPointerDerived<Vector5>;
template class ArrayPointerDerived<Vector6>;

template class ArrayPointerDerived<Matrix2>;
template class ArrayPointerDerived<Matrix3>;
template class ArrayPointerDerived<Matrix4>;
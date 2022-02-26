//////////////////////////////////////////////////////////////////////////
// ArrayIO, here provides more convenient I/O, based on basic functions provided by File.h
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ArrayIO_h__
#define __ArrayIO_h__


#include "File.h"
#include "Common.h"

namespace BinaryDataIO {
	//// Declaration of functions

	//Dynamic array of scalars. Memory layout:
	//Bytes [0,4): an integer n, indicating number of scalars (NOT number of bytes)
	//Bytes [4,4+n*sizeof(T)): binary value of n scalars
	bool Write_Scalar_Array(const std::string& file_name, const Array<bool>& arr);
	template<class T> bool Write_Scalar_Array(const std::string& file_name, const Array<T>& arr);
	template<class T, class F> bool Write_Scalar_Array(const std::string& file_name, F& f, std::uint32_t n);//i-th scalar is f(i)
	bool Read_Scalar_Array(const std::string& file_name, Array<bool>& arr);
	template<class T> bool Read_Scalar_Array(const std::string& file_name, Array<T>& arr);

	//Dynamic array of vectors. Memory layout:
	//Bytes [0,4): an integer n, indicating number of vectors (NOT number of bytes)
	//Bytes [4,4+n*d*sizeof(T)): binary value of n vectors
	template<class T, int d> bool Write_Vector_Array(const std::string& file_name, const Array<Vector<T, d> >& arr);
	template<class T, int d> bool Read_Vector_Array(const std::string& file_name, Array<Vector<T, d> >& arr);


	//Dynamic array of vectors. Memory layout:
	//Bytes [0,4): an integer n, indicating number of vectors (NOT number of bytes)
	//Bytes [4,4+n*3*sizeof(T)): binary value of n 3D vectors.
	//Note: always store 3D vectors in binary file. If the original vectors are less than 3D, we fill it with 0.
	template<class T, int d> bool Write_Vector_Array_3D(const std::string& file_name, const Array<Vector<T, d> >& arr);
	template<class T, int d, class F> bool Write_Vector_Array_3D(const std::string& file_name, F& f, std::uint32_t n);//i-th vector is f(i)
	template<class T, int d> bool Read_Vector_Array_3D(const std::string& file_name, Array<Vector<T, d> >& arr);

	//Dynamic array of Matrices. Memory layout:
	//Bytes [0,4): an integer n, indicating number of matrices (NOT number of bytes)
	//Bytes [4,4+n*d*d*sizeof(T)): binary value of matrices with row-major
	template<class T, int d> bool Write_Matrix_Array(const std::string& file_name, const Array<Matrix<T, d> >& arr);
	template<class T, int d> bool Read_Matrix_Array(const std::string& file_name, Array<Matrix<T, d> >& arr);

	//Dynamic array of different types.
	template<class T> bool Write_Array(const std::string& file_name, const Array<T>& arr) {	return Write_Scalar_Array(file_name, arr);	}
	template<class T> bool Write_Array(const std::string& file_name, const Array<Vector<T,1>> &arr){ return Write_Vector_Array(file_name, arr); }
	template<class T, int d> bool Write_Array(const std::string& file_name, const Array<Vector<T, d>>& arr) { return Write_Vector_Array(file_name, arr); }
	template<class T, int d> bool Write_Array(const std::string& file_name, const Array<Matrix<T, d>>& arr) { return Write_Matrix_Array(file_name, arr); }
	template<class T> bool Read_Array(const std::string& file_name, Array<T>& arr) { return Read_Scalar_Array(file_name, arr); }
	template<class T> bool Read_Array(const std::string& file_name, Array<Vector<T, 1>>& arr) { return Read_Vector_Array(file_name, arr); }
	template<class T, int d> bool Read_Array(const std::string& file_name, Array<Vector<T, d>>& arr) { return Read_Vector_Array(file_name, arr); }
	template<class T, int d> bool Read_Array(const std::string& file_name, Array<Matrix<T, d>>& arr) { return Read_Matrix_Array(file_name, arr); }


	//// Implementation of functions

	template<class T>
	bool Write_Scalar_Array(const std::string& file_name, const Array<T>& arr)
	{
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		std::uint32_t n = (std::uint32_t)arr.size();
		File::Write_Binary<std::uint32_t>(output, n);
		const T* data=(const T*)(arr.data());
		File::Write_Binary_Array<T>(output, data, (int)n);
		output.close();
		return true;
	}

	template<class T, class F>
	bool Write_Scalar_Array(const std::string& file_name, F& f, std::uint32_t n)
	{
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		File::Write_Binary<std::uint32_t>(output, n);
		T* data = new T[n];
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) data[i] = f(i);
		File::Write_Binary_Array<T>(output, data, n);
		delete[] data;
		output.close();
		return true;
	}
	template<class T>
	bool Read_Scalar_Array(const std::string& file_name, Array<T>& arr)
	{
		std::ifstream input(file_name, std::ios::binary);
		if (!input) return false;
		std::uint32_t n;
		File::Read_Binary<std::uint32_t>(input, n);
		arr.resize(n);
		T* data=(T*)(arr.data());
		File::Read_Binary_Array<T>(input, data, (int)n);
		input.close();
		return true;
	}

	template<class T, int d>
	bool Write_Vector_Array(const std::string& file_name, const Array<Vector<T, d>>& arr)
	{
		assert(1 <= d && d <= 3);
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		std::uint32_t n = (std::uint32_t)arr.size();
		File::Write_Binary<std::uint32_t>(output, n);
		T* data = new T[n * d];
		memset(data, 0, (int)n * d * sizeof(T));
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) {
			for (int axis = 0; axis < d; axis++) data[i * d + axis] = arr[i](axis);
		}
		File::Write_Binary_Array(output, data, n * d);
		delete[] data;
		output.close();
		return true;
	}
	template<class T, int d>
	bool Read_Vector_Array(const std::string& file_name, Array<Vector<T, d>>& arr)
	{
		assert(1 <= d && d <= 3);
		std::ifstream input(file_name, std::ios::binary);
		if (!input) return false;
		std::uint32_t n;
		File::Read_Binary(input, n);
		arr.resize(n);
		T* data = new T[n * d];
		File::Read_Binary_Array(input, data, n * d);
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) {
			for (int axis = 0; axis < d; axis++) arr[i](axis) = data[i * d + axis];
		}
		delete[] data;
		input.close();
		return true;
	}
	template<class T, int d>
	bool Write_Vector_Array_3D(const std::string& file_name, const Array<Vector<T, d>>& arr)
	{
		assert(1 <= d && d <= 3);
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		std::uint32_t n = (std::uint32_t)arr.size();
		File::Write_Binary(output, n);
		T* data = new T[n * 3];
		memset(data, 0, n * 3 * sizeof(T));
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) {
			for (int axis = 0; axis < d; axis++) data[i * 3 + axis] = arr[i](axis);
		}
		File::Write_Binary_Array(output, data, n * 3);
		delete[] data;
		output.close();
		return true;
	}
	template<class T, int d, class F>
	bool Write_Vector_Array_3D(const std::string& file_name, F& f, std::uint32_t n)
	{
		assert(1 <= d && d <= 3);
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		File::Write_Binary(output, n);
		T* data = new T[n * 3];
		memset(data, 0, n * 3 * sizeof(T));
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			Vector<T, d> vec = f(i);
			for (int axis = 0; axis < d; axis++) data[i * 3 + axis] = vec(axis);
		}
		File::Write_Binary_Array(output, data, n * 3);
		delete[] data;
		output.close();
		return true;
	}
	template<class T, int d>
	bool Read_Vector_Array_3D(const std::string& file_name, Array<Vector<T, d>>& arr)
	{
		assert(1 <= d && d <= 3);
		std::ifstream input(file_name, std::ios::binary);
		if (!input) return false;
		std::uint32_t n;
		File::Read_Binary(input, n);
		arr.resize(n);
		T* data = new T[n * 3];
		File::Read_Binary_Array(input, data, n * 3);
#pragma omp parallel for
		for (int i = 0; i < (int)n; i++) {
			for (int axis = 0; axis < d; axis++) arr[i](axis) = data[i * 3 + axis];
		}
		delete[] data;
		input.close();
		return true;
	}
	template<class T, int d>
	bool Write_Matrix_Array(const std::string& file_name, const Array<Matrix<T, d>>& arr)
	{
		assert(d > 0);
		std::ofstream output(file_name, std::ios::binary);
		if (!output) return false;
		std::uint32_t n = arr.size();
		File::Write_Binary(output, n);
		T* data = new T[n * d * d];
		memset(data, 0, n * d * d * sizeof(T));
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			std::uint32_t offset = i * d * d;
			memcpy(data + offset, arr[i].data(), d * d * sizeof(T));
		}
		File::Write_Binary_Array(output, data, n * d * d);
		delete[] data;
		output.close();
		return true;
	}
	template<class T, int d>
	bool Read_Matrix_Array(const std::string& file_name, Array<Matrix<T, d>>& arr)
	{
		assert(d > 0);
		std::ifstream input(file_name, std::ios::binary);
		if (!input) return false;
		std::uint32_t n;
		File::Read_Binary(input, n);
		arr.resize(n);
		T* data = new T[n * d * d];
		File::Read_Binary_Array(input, data, n * d * d);
#pragma omp parallel for
		for (int i = 0; i < n; i++) {
			std::uint32_t offset = i * d * d;
			memcpy(arr[i].data(), data + offset, d * d * sizeof(T));
		}
		delete[] data;
		input.close();
		return true;
	}
}

#endif

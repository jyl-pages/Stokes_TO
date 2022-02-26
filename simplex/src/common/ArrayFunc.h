//////////////////////////////////////////////////////////////////////////
// Auxillary functions relevent to Array(i.e., std::vector) operators
// Copyright (c) (2018-), Mengdi Wang
// This file is part of SimpleX, whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#ifndef __ArrayFunc_h__
#define __ArrayFunc_h__
#include "Common.h"
#include "AuxFunc.h"
#include <iostream>

namespace ArrayFunc {
	////Definitions
	//check if an Array contains invalid values, like nan, inf
	template<class T> bool Numerical_Check(const Array<T>& arr, const std::string& name = "", bool crash_on_fail = true);
	template<class T> int Largest_Norm_Element(const Array<T>& arr);//return -1 on empty
	template<class T> real Largest_Norm(const Array<T>& arr);

	////Implementations of template functions
	template<class T>
	bool Numerical_Check(const Array<T>& arr, const std::string& name, bool crash_on_fail)
	{
		int invalid_cnt = 0;
		for (int i = 0; i < arr.size(); i++) {
			//std::cout << "i: " << i << std::endl;
			//std::cout << "thing: " << arr[i] << std::endl;
			if (!Is_Valid_Number(arr[i])) {
				std::cerr << "ArrayFunc::Numerical_Check" << "for [" << name << "] failed at index " << i << " : " << arr[i] << "\n";
				invalid_cnt++;
				if (invalid_cnt >= 10) {
					std::cerr << "......\n";
					break;
				}
			}
		}
		if (invalid_cnt) {
			if (crash_on_fail) exit(1);
			else return false;
		}
		return true;
	}
	template<class T>
	int Largest_Norm_Element(const Array<T>& arr)
	{
		int idx = -1; real max_norm = -1.0;
		for (int i = 0; i < arr.size(); i++) {
			real v = AuxFunc::Norm<T>(arr[i]);
			if (v >= max_norm) {
				idx = i;
				max_norm = v;
			}
		}
		return idx;
	}
	template<class T>
	real Largest_Norm(const Array<T>& arr)
	{
		int idx = Largest_Norm_Element(arr);
		if (idx < 0) return 0;
		else return AuxFunc::Norm<T>(arr[idx]);
	}
}

#endif

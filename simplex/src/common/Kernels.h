//////////////////////////////////////////////////////////////////////////
// SPH Kernels
// Copyright (c) (2018-),Xiangxin Kong, Mengdi Wang
// Please see simplex/docs/kernels-math-en.md for documentation.
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////

#ifndef __Kernel_h__
#define __Kernel_h__
#include <iostream>
#include <cmath>
#include "Common.h"
#include "Constants.h"
#include "AuxFunc.h"

enum class KernelType { NONE, POLY6 = 0, SPIKY, CUBIC, QUINTIC, GAUSSIAN };

class UnitKernel {
public:
	Vector<real, 3> alpha;//0,1,2, -> d=1,2,3
	virtual real Weight(int d, real r)const = 0;
	virtual real Grad(int d, real r)const = 0;
};

class UnitPOLY6 : public UnitKernel {
public:
	UnitPOLY6() { alpha = Vector3(35.0 / 32.0, 4.0 / pi, 315.0 / (64.0 * pi)); }
	virtual real Weight(int d, real r)const { return r < 1 ? alpha[d - 1] * AuxFunc::Quick_Pow(1 - r * r, 3) : 0; }
	virtual real Grad(int d, real r)const { return r < 1 ? alpha[d - 1] * 6 * AuxFunc::Quick_Pow(r * r - 1, 2) * r : 0; }
};

class UnitSPIKY : public UnitKernel {
public:
	UnitSPIKY() { alpha = Vector3(2.0, 10.0 / pi, 15.0 / pi); }
	virtual real Weight(int d, real r)const { return r < 1 ? alpha[d - 1] * AuxFunc::Quick_Pow(1 - r, 3) : 0; }
	virtual real Grad(int d, real r)const { return r < 1 ? alpha[d - 1] * (-3) * AuxFunc::Quick_Pow(1 - r, 2) : 0; }
};

class UnitCUBIC :public UnitKernel {
public:
	UnitCUBIC() { alpha = Vector3(4.0 / 3.0, 40.0 / (7.0 * pi), 8.0 / pi); }
	virtual real Weight(int d, real r)const {
		if (0 <= r && r < 0.5) return alpha[d - 1] * ((real)6 * AuxFunc::Quick_Pow(r, 3) - (real)6 * AuxFunc::Quick_Pow(r, 2) + 1);
		else if (0.5 <= r && r < 1) return alpha[d - 1] * 2 * AuxFunc::Quick_Pow(1 - r, 3);
		else return 0;
	}
	virtual real Grad(int d, real r)const {
		if (0 <= r && r < 0.5) return alpha[d - 1] * 6.0 * r * (3.0 * r - 2.0);
		else if (0.5 <= r && r < 1) return alpha[d - 1] * (-6.0) * AuxFunc::Quick_Pow(1.0 - r, 2);
		else return 0;
	}
};

class UnitQUINTIC :public UnitKernel {
public:
	UnitQUINTIC() { alpha = Vector3(1.0 / 40.0, 63.0 / (478.0 * pi), 81.0 / (359.0 * pi)); }
	virtual real Weight(int d, real r)const {
		if (0 <= r && r < 1.0 / 3) return alpha[d - 1] * (AuxFunc::Quick_Pow(3 - 3 * r, 5) - 6 * AuxFunc::Quick_Pow(2 - 3 * r, 5) + 15 * AuxFunc::Quick_Pow(1 - 3 * r, 5));
		else if (1.0 / 3 <= r && r < 2.0 / 3) return alpha[d - 1] * (AuxFunc::Quick_Pow(3 - 3 * r, 5) - 6 * AuxFunc::Quick_Pow(2 - 3 * r, 5));
		else if (2.0 / 3 < r && r < 1.0) return alpha[d - 1] * (AuxFunc::Quick_Pow(3 - 3 * r, 5));
		else return 0;
	}
	virtual real Grad(int d, real r)const {
		if (0 <= r && r < 1.0 / 3)
			return alpha[d - 1] * (-(real)15 * AuxFunc::Quick_Pow(3 - 3 * r, 4) + (real)90 * AuxFunc::Quick_Pow(2 - 3 * r, 4) - (real)225 * AuxFunc::Quick_Pow(1 - 3 * r, 4));
		else if (1.0 / 3 <= r && r < 2.0 / 3)
			return alpha[d - 1] * (-(real)15 * AuxFunc::Quick_Pow(3 - 3 * r, 4) + (real)90 * AuxFunc::Quick_Pow(2 - 3 * r, 4));
		else if (2.0 / 3 < r && r < 1.0)
			return alpha[d - 1] * (-(real)15 * AuxFunc::Quick_Pow(3 - 3 * r, 4));
		else return 0;
	}
};

class UnitGAUSSIAN :public UnitKernel {
public:
	Vector3 beta;
	real trunc_num;//take sqrt(2)*sigma=1.0/trunc_num
	UnitGAUSSIAN(real _trunc_num = 3) :trunc_num(_trunc_num) {
		alpha = Vector3(trunc_num / sqrt(pi), AuxFunc::Quick_Pow(trunc_num, 2) / pi, AuxFunc::Quick_Pow(trunc_num, 3) / pow(pi, 1.5));
		beta = alpha * (-2) * AuxFunc::Quick_Pow(trunc_num, 2);
	}
	virtual real Weight(int d, real r)const { return alpha[d - 1] * exp(-AuxFunc::Quick_Pow(trunc_num * r, 2)); }
	virtual real Grad(int d, real r)const { return beta[d - 1] * r * exp(-AuxFunc::Quick_Pow(trunc_num * r, 2)); }
};

class KernelSPH {
public:
	static UnitPOLY6 poly6;
	static UnitSPIKY spiky;
	static UnitCUBIC cubic;
	static UnitQUINTIC quintic;
	static UnitGAUSSIAN gaussian;
	bool registered = false;
	std::array<UnitKernel*, 5> kernels;
	
	KernelType ref_type;
	//Truncated at h. Only non-zero in [0,h)
	real h;
	real h_pows[5];//3d in maximum, so we must have h_pows[4]
	KernelSPH(const real _h = 1.0, const KernelType _type = KernelType::SPIKY) :h(_h),ref_type(_type) {
		h_pows[0] = 1;
		for (int i = 1; i < 5; i++) { h_pows[i] = h_pows[i - 1] * h; }
		if (!registered) {
			kernels[(int)KernelType::POLY6] = &KernelSPH::poly6;
			kernels[(int)KernelType::SPIKY] = &KernelSPH::spiky;
			kernels[(int)KernelType::CUBIC] = &KernelSPH::cubic;
			kernels[(int)KernelType::QUINTIC] = &KernelSPH::quintic;
			kernels[(int)KernelType::GAUSSIAN] = &KernelSPH::gaussian;
			registered = true;
		}
	}

	real Weight(int d, real r, real h1, KernelType kernel_type = KernelType::NONE)const;
	real Weight(int d, real r, KernelType kernel_type = KernelType::NONE)const;// const { return Weight(d, r, h, kernel_type); }
	template<int d> real Weight(const Vector<real,d> &r, real h1, const KernelType& kernel_type = KernelType::NONE)const{return Weight(d, r.norm(), h1, kernel_type); }
	template<int d> real Weight(const Vector<real, d>& r, const KernelType& kernel_type = KernelType::NONE)const { return Weight<d>(r, h, kernel_type); }

	template<int d> Vector<real, d> Grad(const Vector<real, d>& r, real h1, KernelType kernel_type = KernelType::NONE)const{
		if (kernel_type == KernelType::NONE) kernel_type = ref_type;
		return kernels[(int)kernel_type]->Grad(d, r.norm() / h1) * r.normalized() / AuxFunc::Quick_Pow(h1, d + 1);
	}
	template<int d> Vector<real, d> Grad(const Vector<real, d>& r, KernelType kernel_type = KernelType::NONE)const {
		if (kernel_type == KernelType::NONE) kernel_type = ref_type;
		return kernels[(int)kernel_type]->Grad(d, r.norm() / h) * r.normalized() / h_pows[d + 1];
	}
};
#endif

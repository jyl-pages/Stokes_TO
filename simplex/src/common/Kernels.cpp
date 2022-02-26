//////////////////////////////////////////////////////////////////////////
// SPH Kernels
// Copyright (c) (2018-),Xiangxin Kong, Mengdi Wang
// Please see simplex/docs/kernels-math-en.md for documentation.
// This file is part of SimpleX,whose distribution is governed by the LICENSE file.
//////////////////////////////////////////////////////////////////////////
#include "Kernels.h"

real KernelSPH::Weight(int d, real r, real h1, KernelType kernel_type)const {
	if (kernel_type == KernelType::NONE) kernel_type = ref_type;
	return kernels[(int)kernel_type]->Weight(d, r / h1) / AuxFunc::Quick_Pow(h1, d);
}

real KernelSPH::Weight(int d, real r, KernelType kernel_type)const {
	if (kernel_type == KernelType::NONE) kernel_type = ref_type;
	return kernels[(int)kernel_type]->Weight(d, r / h) / h_pows[d];
}

UnitPOLY6 KernelSPH::poly6;
UnitSPIKY KernelSPH::spiky;
UnitCUBIC KernelSPH::cubic;
UnitQUINTIC KernelSPH::quintic;
UnitGAUSSIAN KernelSPH::gaussian;

#pragma once
#include "LinearMapping.h"
#include "cuda_runtime_api.h"
#include "gpuUtils.h"

template<class Attr>
class PartialMapping :public LinearMapping
{
public:
	Attr *attr_x, *attr_y;
	LinearMapping *mapping;
	int xdof, ydof;
	Scalar *temp;

	PartialMapping(Attr *_attr_y, LinearMapping *_mapping, Attr *_attr_x, Scalar *buffer=nullptr)
	{
		attr_x = _attr_x;
		attr_y = _attr_y;
		mapping = _mapping;
		xdof = mapping->xDoF();
		ydof = mapping->yDoF();
		if (buffer == nullptr)
			cudaMalloc(&temp, sizeof(Scalar)*xdof);
		else
			temp = buffer;
	}

	int xDoF() override
	{
		return xdof;
	}

	int yDoF() override
	{
		return ydof;
	}

	void applyMapping(Scalar *Ap, Scalar *p) override
	{
		cudaMemset(Ap, 0, sizeof(Scalar)*ydof);
		auto fix = [=] __device__(Scalar& tv, bool tmask) {if (!tmask) tv = (Scalar)0; };
		cudaMemcpy(temp, p, sizeof(Scalar)*xdof, cudaMemcpyDeviceToDevice);
		cwise_mapping_wrapper(temp, *attr_x, fix, xdof);
		mapping->applyMapping(Ap, temp);
		cwise_mapping_wrapper(Ap, *attr_y, fix, ydof);
	}
};
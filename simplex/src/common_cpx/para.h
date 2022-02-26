#pragma once
#include "cuda_runtime_api.h"
#include <cstdio>

#define USE_DOUBLE

#ifdef USE_DOUBLE
typedef double Scalar;
#else
typedef float Scalar;
#endif

#define CUDA_DEBUG

static const char* _cudaGetErrorEnum(cudaError_t error) {
	return cudaGetErrorName(error);
}

template <typename T>
void check(T result, char const* const func, const char* const file,
	int const line) {
	if (result) {
		fprintf(stderr, "CUDA error at %s:%d code=%d(%s) \"%s\" \n", file, line,
			static_cast<unsigned int>(result), _cudaGetErrorEnum(result), func);
		cudaDeviceReset();
		// Make sure we call CUDA Device Reset before exiting
		exit(1);
	}
}

#define checkCudaErrors(val) check((val), #val, __FILE__, __LINE__)

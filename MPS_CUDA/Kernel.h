#ifndef __KERNEL_H__
#define __KERNEL_H__

#include <cuda_runtime.h>

class kernel
{
public:
	__host__ __device__ double weight(double r, double re);

	__host__ __device__ kernel();
	__host__ __device__ ~kernel();
	__host__ __device__ kernel(const kernel &k);

};

#endif
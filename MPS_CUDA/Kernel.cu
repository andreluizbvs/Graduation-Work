#include "Kernel.h"

__host__ __device__ double kernel::weight(double r, double re)
{
	double ORIG_rr;
	double ORIG_ww;

	ORIG_rr = (r / re);
	if (ORIG_rr<1.0) ORIG_ww = 1.0 / ORIG_rr - 1.0;
	else ORIG_ww = 0.0;

	return(ORIG_ww);
}

__host__ __device__ kernel::kernel()
{

}

__host__ __device__ kernel::~kernel()
{

}

__host__ __device__ kernel::kernel(const kernel &k)
{

}
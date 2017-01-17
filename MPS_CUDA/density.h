#ifndef __DENSITY_H__
#define __DENSITY_H__

#include "Particles.h"
#include "Neighborhood.h"
#include <cmath>
#include "Kernel.h"
#include"cuda_runtime.h"

#include <vector>

using namespace std;

__global__ void set_n_kernel2D(Particle2D* particles, double* n0p);
__global__ void cal_n_kernel2D(int offset, int * nei_1D, Particle2D* particles, double r, int nump);
__global__ void set_n_kernel(Particle3D* particles, double* n0p);
__global__ void cal_n_kernel(int offset, int * nei_1D, Particle3D* particles, double r, int nump);

class density
{
public:
//protected:
	kernel k;
	
//public:
	__host__ __device__ void cal_n2D(neighbor nei, Particle2D* particles, double r, int nump);
	__host__ __device__ void cal_n(neighbor nei, Particle3D* particles, double r, int nump);

	__host__ __device__ density();
	__host__ __device__ ~density();
	__host__ __device__ density(const density &d);
};

#endif
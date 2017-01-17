#ifndef __NEIGHBORHOOD_H__
#define __NEIGHBORHOOD_H__

#include "Particles.h"
#include <cuda_runtime.h>
#include <vector>

using namespace std;

__global__ void set_nei_kernel2D(int offset, Particle2D* particles, double r2, int nump, int *nei);
__global__ void set_nei_kernel(int offset, Particle3D* particles, double r2, int nump, int *nei);

class neighbor
{
public:
//protected:
	int** neighbors;

//public:
	__host__ __device__ neighbor(int cont, int max_nei);
	__host__ __device__ ~neighbor();
	__host__ __device__ neighbor(const neighbor&n);

	__host__ __device__ void set2D(Particle2D* particles, double r2, int nump);
	__host__ __device__ void set(Particle3D* particles, double r2, int nump);
	//__host__ __device__ void set_one_part(Particle2D* particles, double r2, int nump, int i);
	//void set_kernel(int offset, Particle2D* particles, double r2, int nump);
	__host__ __device__ int **get();
	__host__ __device__ void clear(int nump);
};


#endif
#ifndef __PARTICLESACTIONS_H__
#define __PARTICLESACTIONS_H__

#include "Particles.h"
#include "Neighborhood.h"
#include <vector>
#include <cmath>
#include <cuda_runtime.h>
#include "SystemResolution.h"


using namespace std;

__global__ void external_force_kernel2D(int offset, Particle2D* particles, double *dt, double g);
__global__ void mov_part_kernel2D(int offset, Particle2D* particles, double *dt);
__global__ void collision_kernel2D(int offset, Particle2D* particles, double *dt, int *nei_1D, double radius, double Density, int nump);
__global__ void set_pr_kernel2D(int offset, Particle2D* particles, double *dp);
__global__ void external_force_kernel(int offset, Particle3D* particles, double *dt, double g);
__global__ void mov_part_kernel(int offset, Particle3D* particles, double *dt);
__global__ void collision_kernel(int offset, Particle3D* particles, double *dt, int *nei_1D, double radius, double Density, int nump);
__global__ void set_pr_kernel(int offset, Particle3D* particles, double *dp);

class ParticlesActions
{
public:
	void external_force(Particle2D* particles, double dt, double g, int nump);
	//void external_force_cu(Particle2D* particles, double dt, double g, int nump);
	//void external_force_HL(Particle2D* particles, double dt, double g, neighbor neiICCG, double radius_ICCG, double n0pICCG, double rho, int nump);
	void mov_part(Particle2D* particles, double dt, int nump);
	void collision(Particle2D* particles, double dt, neighbor nei, double radius, double Density, int nump);
	
	__host__ __device__ ParticlesActions();
	__host__ __device__ ~ParticlesActions();
	__host__ __device__ ParticlesActions(const ParticlesActions& p);
	
};

#endif
#ifndef __SYSTEMRESOLOTION_H__
#define __SYSTEMRESOLOTION_H__

#include "Particles.h"
#include "Neighborhood.h"
#include "Kernel.h"
#include "Neighborhood.h"
#include <vector>
#include"cuda_runtime.h"

using namespace std;

__global__ void conf_dp_kernel(double* dp, double* dpArray, int size, int *bcon, int *contB);
__global__ void set_bcon_kernel2D(int offset, Particle2D* particles, double n0p, double dirichlet, int nump, int *bcon);
__global__ void cal_B_HS_ECS_kernel2D(int offset, Particle2D* particles, double n0p, double *dt, int *nei, double radius, int nump, int *bcon, int *contTemp, double *temp_b);

__global__ void set_bcon_kernel(int offset, Particle3D* particles, double n0p, double dirichlet, int nump, int *bcon);
__global__ void cal_B_HS_ECS_kernel(int offset, Particle3D* particles, double n0p, double *dt, int *nei, double radius, int nump, int *bcon, int *contTemp, double *temp_b);
__global__ void set_B_HS_ECS_kernel(int offset, int *bcon, int *contB, double *scrB, double *temp_b, int *contTemp);
__global__ void cal_B_kernel(int offset, Particle3D* particles, double n0p, double *dt, int *contTemp, double *temp_b);

__global__ void cal_A_HL_kernel2D(int offset, Particle2D* particles, int *neiICCG, double n0pICCG, double radius_ICCG, double rho, int num_d, int nump, int *bcon, double *temp_a);
__global__ void cal_A_HL_kernel(int offset, Particle3D* particles, int *neiICCG, double n0pICCG, double radius_ICCG, double rho, int num_d, int nump, int *bcon, double *temp_a);
__global__ void set_A_HL_kernel1(int offset, int *bcon, double *temp_a2, double *temp_a, int nump, int *bcon_cont);
__global__ void set_A_HL_kernel2(int offsetX, int offsetY, double *srcA, int *contB, double *temp_a2, int nump);

class SystemResolution
{
//protected:
public:
	//vector<double> B;
	double *B;
	int* Bcon;
	double* dp;
	vector<vector<double>> A;
	//double *srcA;
	//kernel k;

//public:
	__host__ __device__ SystemResolution(Particle2D* particles, int nump);
	__host__ __device__ ~SystemResolution();
	__host__ __device__ SystemResolution(const SystemResolution &s);

	void clear(int nump);

	void set_bcon(Particle2D* particles, double n0p, double dirichlet, int nump);
	void set_B(Particle2D* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho, int nump);
	void set_B_HS_ECS(Particle2D* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho, int *nei, double radius, int nump, int *bcon);
	void set_A(Particle2D* particles, neighbor nei, double n0pICCG, double radius, double lambda, double rho, double dt, int num_d, int nump);
	void set_A_HL(Particle2D* particles, int *neiICCG, double n0pICCG, double radius_ICCG, double lambda, double rho, double dt, int num_d, int nump, int* bcon/*, double *srcA*/);
	__host__ __device__ int* GetBcon();
	__host__ __device__ double* get_dp();
	
	void system_resolution(int MaxIterPress, int MinIterPress, double Convergence, int allsize, int *bcon, int contB, double *srcB, double *srcA);

//protected:
	__host__ __device__ void conf_dp(double *dp, int size, int *bcon, int contB);


};

#endif
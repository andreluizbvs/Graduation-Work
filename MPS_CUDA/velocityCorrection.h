#ifndef __velocityCorrection_h__
#define __velocityCorrection_h__

#include "Particles.h"
#include "Neighborhood.h"
#include "Kernel.h"
#include <vector>
#include "cuda_runtime.h"

__global__ void set_pmin_kernel2D(int offset, Particle2D* particles, int *nei, int nump, double *pmin);
__global__ void set_ddv_kernel2D(int offset, Particle2D* particles, int *nei, double rho, double *dt, double n0p, double radius, double num_D, int nump, double *pmin, Point2D *DdvCorrector);
__global__ void updateVelocity_kernel2D(int offset, Particle2D* particles, double* dt, int nump, Point2D *DdvCorrector);

__global__ void set_pmin_kernel(int offset, Particle3D* particles, int *nei, int nump, double *pmin);
__global__ void set_ddv_kernel(int offset, Particle3D* particles, int *nei, double rho, double *dt, double n0p, double radius, double num_D, int nump, double *pmin, Point3D *DdvCorrector);
__global__ void updateVelocity_kernel(int offset, Particle3D* particles, double* dt, int nump, Point3D *DdvCorrector);

class velocityCorrection
{
protected:
	//double *pmin;
	kernel k;

public:
	velocityCorrection();
	~velocityCorrection();
	velocityCorrection(velocityCorrection &v);

	void UpdateVelocity(Particle2D* particles, neighbor nei,double rho,double dt, double n0p, double radius, double num_D, int nump);

};

#endif
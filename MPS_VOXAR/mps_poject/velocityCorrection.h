#ifndef __velocityCorrection_h__
#define __velocityCorrection_h__

#include "Particles.h"
#include "Neighborhood.h"
#include "Kernel.h"
#include <vector>


class velocitycorrection
{
protected:
	vector <double> pmin;
	kernel k;

public:
	velocitycorrection();
	~velocitycorrection();
	velocitycorrection(velocitycorrection &v);

	void UpdateVelolicty(vector<Particle2D>* particles, neighbor nei,double rho,double dt, double n0p, double radius, double num_D);

};

#endif
#ifndef __DENSITY_H__
#define __DENSITY_H__

#include "Particles.h"
#include "Neighborhood.h"
#include <cmath>
#include "Kernel.h"

#include <vector>

using namespace std;

class density
{
protected:
	kernel k;
public:
	void cal_n(neighbor nei, vector<Particle2D>* particles,double r);
	
	density();
	~density();
	density(const density &d);
};

#endif
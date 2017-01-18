#ifndef __PARTICLESACTIONS_H__
#define __PARTICLESACTIONS_H__

#include "Particles.h"
#include "Neighborhood.h"
#include <vector>
#include <cmath>

using namespace std;

class ParticlesActions
{
public:
	void external_force(vector<Particle2D>* particles, double dt, double g);
	void mov_part(vector<Particle2D>* particles, double dt);
	void collision(vector<Particle2D>* particles, double dt, neighbor nei, double radius,double Density);
	
	ParticlesActions();
	~ParticlesActions();
	ParticlesActions(const ParticlesActions& p);
	
};

#endif
#ifndef __NEIGHBORHOOD_H__
#define __NEIGHBORHOOD_H__

#include "Particles.h"
#include <vector>

using namespace std;

class neighbor
{
protected:
	vector<vector<int>> neighbors;

public:
	neighbor(int cont);
	~neighbor();
	neighbor(const neighbor&n);

	void set(vector<Particle2D>* particles,double r2);
	vector<vector<int>> get();
	void clear();
};


#endif
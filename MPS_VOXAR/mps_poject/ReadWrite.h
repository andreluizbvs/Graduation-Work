#ifndef __READWRITE_H__
#define __READWRITE_H__

#include "Particles.h"
#include "Particles.h"
#include "data_in.h"
#include <fstream>
#include <vector>


class ReadWrite
{
public:
	ReadWrite();
	ReadWrite(const ReadWrite&r);
	~ReadWrite();

	void ReadParticles2D(vector<Particle2D>* particles,char* name);
	
	void ReadData(char *name, data_in *dados);

	void WriteOut(vector<Particle2D>* particles,int number);
};

#endif
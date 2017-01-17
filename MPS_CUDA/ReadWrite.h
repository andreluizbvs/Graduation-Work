#ifndef __READWRITE_H__
#define __READWRITE_H__

#include "Particles.h"
#include "data_in.h"
#include <fstream>
#include <vector>
#include <iostream>

class ReadWrite
{
public:
	ReadWrite();
	ReadWrite(const ReadWrite&r);
	~ReadWrite();

	void ReadParticles2D(Particle2D* particles,char* name, int nump);
	void ReadParticles3D(Particle3D* particles, char* name, int nump);
	
	void ReadData(char *name, data_in *dados);

	void WriteOut2D(Particle2D* particles, int number, int nump);
	void WriteOut(Particle3D* particles, int number, int nump);
};

#endif
#include "Neighborhood.h"

// Cálculo dos vizinhos: A estrutura usada para armazenar é matriz simples.
// A organização na matriz é: Primeiro elemento da linha "i" indica o número de vizinhos da partícula "i".
// O resto da linha "i" contém o índice do vizinho da partícula "i"
// [n_vizinhos,v1,v2,v3,v4,v5,...,v8]

__global__ void set_nei_kernel(int offset, Particle3D* particles, double r2, int nump, int *nei)
{
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	//unsigned int i = offset + (/*blockDim.x **/ blockIdx.x /*+ threadIdx.x*/);
	//unsigned int j = offset + (/*blockDim.x * blockIdx.x +*/ threadIdx.x);

	Point3D Part_j, Part_i;
	double x, y, z;

	for (int j = 0; j < nump; j++){
		if (!((particles[j].is_wall() && particles[i].is_wall()) || (i == j))){
			Part_j = (particles)[j].get_po();
			Part_i = (particles)[i].get_po();

			x = Part_j.x - Part_i.x;
			y = Part_j.y - Part_i.y;
			z = Part_j.z - Part_i.z;

			if (((x*x) + (y*y) + (z*z)) <= r2)
			{
				nei[(i * nump) + 0]++;
				nei[(i * nump) + nei[(i * nump) + 0]] = j;
			}
		}
	}
}

__host__ __device__ neighbor::neighbor(const neighbor &n)
{
	this->neighbors = n.neighbors;
}

__host__ __device__ neighbor::neighbor(int cont, int max_nei)
{
	neighbors = new int*[cont];
	for (int i = 0; i < cont; ++i) neighbors[i] = new int[max_nei];
	for (int i = 0; i < cont; i++) {
		memset(neighbors[i], 0, sizeof(neighbors[i][0]) * max_nei);
	}
}

__host__ __device__ neighbor::~neighbor()
{

}

__host__ __device__ void neighbor::clear(int nump)
{
	for (int i = 0; i < nump; i++) {
		memset(neighbors[i], 0, sizeof(neighbors[i][0]) * nump);
	}
}

__host__ __device__ int** neighbor::get()
{
	return neighbors;
}

__host__ __device__ void neighbor::set(Particle3D* particles, double r2, int nump)
{
	Point3D Part_j, Part_i;
	double x, y, z;

	for (int i = 0; i<nump; i++)
	{
		for (int j = i + 1; j<nump; j++)
		{
			if ((particles)[j].is_wall() && (particles)[i].is_wall())
				continue;

			Part_j = (particles)[j].get_po();
			Part_i = (particles)[i].get_po();

			x = Part_j.x - Part_i.x;
			y = Part_j.y - Part_i.y;
			z = Part_j.z - Part_i.z;

			if (((x*x) + (y*y) + (z*z)) <= r2)
			{
				neighbors[i][0]++;
				neighbors[j][0]++;

				neighbors[i][neighbors[i][0]] = j;
				neighbors[j][neighbors[j][0]] = i;
			}
		}
	}
}

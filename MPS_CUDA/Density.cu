#include "density.h"

__global__ void set_n_kernel(Particle3D* particles, double* n0p)
{
	n0p[0] = particles[0].get_n();
}

__global__ void cal_n_kernel(int offset, int* nei_1D, Particle3D* particles, double r, int nump)
{
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	int j = 0;
	double dx = 0.0, dy = 0.0, dz = 0.0, dist = 0.0, kern = 0.0, n = 0.0;
	Point3D Part_j, Part_i;
	kernel k;

	n = 0.0;

	for (int a = 1; a <= nei_1D[(i*nump) + 0]; a++)
	{
		j = nei_1D[(i *nump) + a];

		Part_j = particles[j].get_po();
		Part_i = particles[i].get_po();

		dx = Part_j.x - Part_i.x;
		dy = Part_j.y - Part_i.y;
		dz = Part_j.z - Part_i.z;

		dist = sqrt(dx*dx + dy*dy + dz*dz);

		kern = k.weight(dist, r);

		n += kern;
	}
	particles[i].set_n(n);
}

__host__ __device__ void density::cal_n2D(neighbor nei, Particle2D* particles, double r, int nump)
{
	int j = 0;
	double dx = 0.0, dy = 0.0, dist = 0.0, k = 0.0, n = 0.0;
	Point2D Part_j, Part_i;
	for (int i = 0; i < nump; i++)
	{
		n = 0.0;

		for (int a = 1; a <= nei.get()[i][0]; a++)
		{
			j = nei.get()[i][a];

			Part_j = particles[j].get_po();
			Part_i = particles[i].get_po();

			dx = Part_j.x - Part_i.x;
			dy = Part_j.y - Part_i.y;

			dist = sqrt(dx*dx + dy*dy);

			k = this->k.weight(dist, r);

			n += k;
		}
		particles[i].set_n(n);
	}
}

__host__ __device__ density::density()
{

}

__host__ __device__ density::~density()
{

}

__host__ __device__ density::density(const density &d)
{
	this->k = d.k;
}
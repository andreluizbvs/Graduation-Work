#include "ParticlesActions.h"

__global__ void set_pr_kernel(int offset, Particle3D* particles, double *dp)
{
	unsigned int a = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	particles[a].set_pr(dp[a]);
}

__global__ void external_force_kernel(int offset, Particle3D* particles, double *dt, double g)
{
	unsigned int a = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	if (particles[a].is_fluid()){
		particles[a].set_v(particles[a].get_v().x, particles[a].get_v().y - g*dt[0], (particles[a].get_v().z));
	}
}

__global__ void mov_part_kernel(int offset, Particle3D* particles, double *dt)
{
	unsigned int a = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	if ((particles)[a].is_fluid()){
		(particles)[a].set_po((particles)[a].get_po().x + (particles)[a].get_v().x*dt[0], (particles)[a].get_po().y + (particles)[a].get_v().y*dt[0], (particles)[a].get_po().z + (particles)[a].get_v().z*dt[0]);
	}
}

__global__ void collision_kernel(int offset, Particle3D* particles, double *dt, int *nei_1D, double radius, double Density, int nump)
{
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	double dx = 0.0, dy = 0.0, dz = 0.0, dist = 0.0;
	double vg[3], vr[3], vbas, vm[3];
	double m1 = 0.0, m2 = 0.0;
	Point3D Part_j, Part_i, Velocity_j, Velocity_i;

	int j = 0;;

	m1 = Density;

	for (int a = 1; a <= nei_1D[(i*nump) + 0]; a++)
	{
		j = nei_1D[(i*nump) + a];

		if (j <= i)
			continue;

		Part_j = (particles)[j].get_po();
		Part_i = (particles)[i].get_po();

		dx = Part_j.x - Part_i.x;
		dy = Part_j.y - Part_i.y;
		dz = Part_j.z - Part_i.z;

		dist = dx*dx + dy*dy + dz*dz;
		if (dist<radius)
		{
			Velocity_j = (particles)[j].get_v();
			Velocity_i = (particles)[i].get_v();

			dist = sqrt(dist);

			m2 = Density;
			vg[0] = (m1*Velocity_i.x + m2*Velocity_j.x) / (m1 + m2);
			vg[1] = (m1*Velocity_i.y + m2*Velocity_j.y) / (m1 + m2);
			vg[2] = (m1*Velocity_i.z + m2*Velocity_j.z) / (m1 + m2);

			vr[0] = m1*(Velocity_i.x - vg[0]);
			vr[1] = m1*(Velocity_i.y - vg[1]);
			vr[2] = m1*(Velocity_i.z - vg[2]);

			vbas = (vr[0] * dx + vr[1] * dy + vr[2] * dz) / dist;

			if (vbas<0.0)continue;

			vm[0] = (1.2)*vbas*dx / dist;
			vm[1] = (1.2)*vbas*dy / dist;
			vm[2] = (1.2)*vbas*dz / dist;

			if ((particles)[i].is_fluid())
			{
				(particles)[i].set_v(Velocity_i.x - vm[0] / m1, Velocity_i.y - vm[1] / m1, Velocity_i.z - vm[2] / m1);
				(particles)[i].set_po(Part_i.x - dt[0] * vm[0] / m1, Part_i.y - dt[0] * vm[1] / m1, Part_i.z - dt[0] * vm[2] / m1);
			}

			if ((particles)[j].is_fluid())
			{
				(particles)[j].set_v(Velocity_j.x + vm[0] / m2, Velocity_j.y + vm[1] / m2, Velocity_j.z + vm[2] / m2);
				(particles)[j].set_po(Part_j.x + dt[0] * vm[0] / m2, Part_j.y + dt[0] * vm[1] / m2, Part_j.z + dt[0] * vm[2] / m2);
			}
		}
	}
}

//void ParticlesActions::external_force_cu(Particle2D* particles, double dt, double g, int nump){
//
//	Particle2D *particles_d = NULL;
//	cudaError_t err = cudaSuccess;
//	err = cudaMalloc((void **)&particles_d, nump*sizeof(Particle2D));
//	
//	if (err != cudaSuccess)
//	{
//		fprintf(stderr, "Failed to allocate device particles (error code %s)!\n", cudaGetErrorString(err));
//		exit(EXIT_FAILURE);
//	}
//	
//	int division = nump / 1024;
//	int mod = nump - (division * 1024);
//	err = cudaMemcpy(particles_d, particles, nump*sizeof(Particle2D), cudaMemcpyHostToDevice);
//	
//	if (err != cudaSuccess)
//	{
//		fprintf(stderr, "Failed to copy particles from host to device (error code %s)!\n", cudaGetErrorString(err));
//		exit(EXIT_FAILURE);
//	}
//	
//	external_force_kernel << <division, 1024 >> >(0, particles_d, dt, g);
//	external_force_kernel << <1, mod >> >(division * 1024, particles_d, dt, g);
//	
//	err = cudaGetLastError();
//	if (err != cudaSuccess)
//	{
//		fprintf(stderr, "Failed to launch external force kernel (error code %s)!\n", cudaGetErrorString(err));
//		system("pause");
//		exit(EXIT_FAILURE);
//	}
//	
//	err = cudaMemcpy(particles, particles_d, nump*sizeof(Particle2D), cudaMemcpyDeviceToHost);
//	if (err != cudaSuccess)
//	{
//		fprintf(stderr, "Failed to copy particles from device to host (error code %s)!\n", cudaGetErrorString(err));
//		exit(EXIT_FAILURE);
//	}
//	
//	err = cudaFree(particles_d);
//
//	if (err != cudaSuccess)
//	{
//		fprintf(stderr, "Failed to free device vector A (error code %s)!\n", cudaGetErrorString(err));
//		exit(EXIT_FAILURE);
//	}
//	
//	err = cudaDeviceReset();
//
//	if (err != cudaSuccess)
//	{
//		fprintf(stderr, "Failed to deinitialize the device! error=%s\n", cudaGetErrorString(err));
//		exit(EXIT_FAILURE);
//	}
//
//	//return particles;
//}
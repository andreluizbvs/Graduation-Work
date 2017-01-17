#include "velocityCorrection.h"

// Correção das velocidades e das posições das partículas de acordo com os valores de pressão calculados.

__global__ void set_pmin_kernel(int offset, Particle3D* particles, int *nei, int nump, double *pmin){
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	
	//double *pmin = new double[nump];
	double AuxPmin = 0.0;
	int j = 0;
	double p;

	if (particles[i].is_fluid()){
		AuxPmin = (particles)[i].get_pr();
		for (int a = 1; a <= nei[(i*nump) + 0]; a++){
			j = nei[(i*nump) + a];

			if ((particles)[j].is_wall()){
			}
			else{
				p = (particles)[j].get_pr();

				if (p < AuxPmin)
					AuxPmin = p;
			}
		}
		pmin[i] = (AuxPmin);		
	}
	else{
		pmin[i] = 0.0;
	}
}

__global__ void set_ddv_kernel(int offset, Particle3D* particles, int *nei, double rho, double *dt, double n0p, double radius, double num_D, int nump, double *pmin, Point3D *DdvCorrector){
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	Point3D Part_j, Part_i, new_point;
	int j = 0;
	double dx = 0.0, dy = 0.0, dz = 0.0, dist = 0.0, ddv = 0.0;
	//Point2D *DdvCorrector = new Point2D[nump];
	kernel k;

	new_point.x = 0.0;
	new_point.y = 0.0;
	new_point.z = 0.0;

	DdvCorrector[i] = (new_point);

	if ((particles)[i].is_fluid()){
		for (int a = 1; a <= nei[(i*nump) + 0]; a++){
			j = nei[(i*nump) + a];

			if (!particles[j].is_wall()){
				dx = 0.0;
				dy = 0.0;
				dz = 0.0;
				dist = 0.0;
				ddv = 0.0;

				Part_j = (particles)[j].get_po();
				Part_i = (particles)[i].get_po();

				dx = Part_i.x - Part_j.x;
				dy = Part_i.y - Part_j.y;
				dz = Part_i.z - Part_j.z;

				dist = sqrt(dx*dx + dy*dy + dz*dz);

				ddv = dt[0];

				ddv = ddv*(((particles)[i].get_pr() + (particles)[j].get_pr()) - (pmin[i] + pmin[j]));

				ddv = ddv / dist *k.weight(dist, radius);
				ddv = ddv / rho;
				ddv = ddv / n0p;
				ddv = ddv* num_D;

				DdvCorrector[i].x += ((ddv * dx) / dist);
				DdvCorrector[i].y += ((ddv * dy) / dist);
				DdvCorrector[i].z += ((ddv * dz) / dist);
			}
		}
	}
}

__global__ void updateVelocity_kernel(int offset, Particle3D* particles, double* dt, int nump, Point3D *DdvCorrector){
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	
	Point3D v_i, Part_i;

	Part_i = (particles)[i].get_po();
	v_i = (particles)[i].get_v();

	if ((particles)[i].is_fluid())
	{
		(particles)[i].set_po(Part_i.x + DdvCorrector[i].x*dt[0], Part_i.y + DdvCorrector[i].y*dt[0], Part_i.z + DdvCorrector[i].z*dt[0]);
		(particles)[i].set_v(v_i.x + DdvCorrector[i].x, v_i.y + DdvCorrector[i].y, v_i.z + DdvCorrector[i].z);
	}
}

velocityCorrection::velocityCorrection()
{

}

velocityCorrection::~velocityCorrection()
{

}

velocityCorrection::velocityCorrection(velocityCorrection &v)
{
	//this->pmin = v.pmin;
	this->k = k;
}

void velocityCorrection::UpdateVelocity(Particle2D* particles, neighbor nei, double rho, double dt, double n0p, double radius, double num_D, int nump)
{
	double *pmin = new double[nump];
	double AuxPmin = 0.0;
	int j = 0;

	Point2D* DdvCorrector = new Point2D[nump];

	for (int i = 0; i < nump; i++)
	{
		if (!(particles)[i].is_fluid()){
			pmin[i] = 0.0;
			continue;
		}

		AuxPmin = (particles)[i].get_pr();
		for (int a = 1; a <= nei.get()[i][0]; a++)
		{
			j = nei.get()[i][a];

			if ((particles)[j].is_wall())
				continue;

			double p = (particles)[j].get_pr();

			if (p < AuxPmin)
				AuxPmin = p;
		}
		pmin[i] = (AuxPmin);
	}

	Point2D Part_j, Part_i;
	Point2D new_point;
	double dx = 0.0, dy = 0.0, dist = 0.0,/* kernel = 0.0,*/ ddv = 0.0;

	for (int i = 0; i < nump; i++)
	{
		new_point.x = 0.0;
		new_point.y = 0.0;

		DdvCorrector[i] = (new_point);

		if (!(particles)[i].is_fluid())
			continue;

		for (int a = 1; a <= nei.get()[i][0]; a++)
		{
			j = nei.get()[i][a];

			if ((particles)[j].is_wall())
				continue;

			dx = 0.0;
			dy = 0.0;
			dist = 0.0;
			//kernel = 0.0;
			ddv = 0.0;

			Part_j = (particles)[j].get_po();
			Part_i = (particles)[i].get_po();

			dx = Part_i.x - Part_j.x;
			dy = Part_i.y - Part_j.y;

			dist = sqrt(dx*dx + dy*dy);

			ddv = dt;

			ddv = ddv*(((particles)[i].get_pr() + (particles)[j].get_pr()) - (pmin[i] + pmin[j]));

			ddv = ddv / dist *k.weight(dist, radius);
			ddv = ddv / rho;
			ddv = ddv / n0p;
			ddv = ddv* num_D;

			DdvCorrector[i].x += ((ddv * dx) / dist);
			DdvCorrector[i].y += ((ddv * dy) / dist);
		}
	}

	Point2D v_i;

	for (int i = 0; i < nump; i++)
	{
		Part_i = (particles)[i].get_po();
		v_i = (particles)[i].get_v();

		if ((particles)[i].is_fluid())
		{
			(particles)[i].set_po(Part_i.x + DdvCorrector[i].x*dt, Part_i.y + DdvCorrector[i].y*dt);
			(particles)[i].set_v(v_i.x + DdvCorrector[i].x, v_i.y + DdvCorrector[i].y);
		}
	}
	delete[] pmin;
	delete[] DdvCorrector;
}

#include "ParticlesActions.h"
#include <iostream>


//funções de interações fisicas basicas
// acrecimo de força externa (gravidade)
// deslocamento de particulas
//tratamendo de colisões

ParticlesActions::ParticlesActions()
{

}

ParticlesActions::ParticlesActions(const ParticlesActions&p)
{

}

ParticlesActions::~ParticlesActions()
{

}

void ParticlesActions::external_force(vector<Particle2D>* particles, double dt,double g)
{
	Point2D velocity;
	for(int a=0; a<particles->size(); a++)
	{
		if((*particles)[a].is_fluid())
		{
			velocity = (*particles)[a].get_v();
			(*particles)[a].set_v(velocity.x,velocity.y - g*dt);
		}
	}
}

void ParticlesActions::mov_part(vector<Particle2D>* particles, double dt)
{
	Point2D position;
	Point2D velocity;
	for(int a=0; a<particles->size(); a++)
	{
		if((*particles)[a].is_fluid())
		{
			position = (*particles)[a].get_po();
			velocity = (*particles)[a].get_v();
			(*particles)[a].set_po(position.x + velocity.x*dt,position.y + velocity.y*dt);
		}
	}
}

void ParticlesActions::collision(vector<Particle2D>* particles, double dt, neighbor nei, double radius, double Density)
{
	double dx = 0.0, dy = 0.0, dist = 0.0;
	double vg[2],vr[2],vbas,vm[2];
	double m1 = 0.0, m2 = 0.0;
	vector<vector<int>> nei_in = nei.get();
	int j = 0;

	for(int i=0; i<particles->size(); i++)
	{
		m1 = Density;

		for(int a=1; a<=nei_in[i][0]; a++)
		{
			j = nei_in[i][a];

			if(j<=i)
				continue;

			Point2D Part_j = (*particles)[j].get_po();
			Point2D	Part_i = (*particles)[i].get_po();

			dx = Part_j.x - Part_i.x;
			dy = Part_j.y - Part_i.y;

			dist = dx*dx + dy*dy;
			if(dist<radius)
			{
				//printf("prepare yourself %d %d\n", i, j);
				//printf("prepare yourself %d %d\n", i, j);
				Point2D Velocity_j = (*particles)[j].get_v();
				Point2D Velocity_i = (*particles)[i].get_v();

				dist = sqrt(dist);

				m2 = Density;
				vg[0] = (m1*Velocity_i.x + m2*Velocity_j.x ) / (m1 + m2);
				vg[1] = (m1*Velocity_i.y + m2*Velocity_j.y ) / (m1 + m2);

				vr[0] = m1*(Velocity_i.x - vg[0]);
				vr[1] = m1*(Velocity_i.y - vg[1]);

				vbas = (vr[0]*dx + vr[1]*dy)/dist;

				if(vbas<0.0)continue;

				vm[0] = (1.2)*vbas*dx/dist;
				vm[1] = (1.2)*vbas*dy/dist;

				if((*particles)[i].is_fluid())
				{
					(*particles)[i].set_v( Velocity_i.x - vm[0]/m1, Velocity_i.y - vm[1]/m1 );
					(*particles)[i].set_po( Part_i.x - dt*vm[0]/m1, Part_i.y - dt*vm[1]/m1 );
				}

				if((*particles)[j].is_fluid())
				{
					(*particles)[j].set_v( Velocity_j.x + vm[0]/m2, Velocity_j.y + vm[1]/m2 );
					(*particles)[j].set_po( Part_j.x + dt*vm[0]/m2, Part_j.y + dt*vm[1]/m2 );
				}

			}
		}
	}
}
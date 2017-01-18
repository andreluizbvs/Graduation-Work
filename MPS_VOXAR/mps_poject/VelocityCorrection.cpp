#include "velocityCorrection.h"


//correção da velocidade e da posição das particulas de acordo com os valores de pressão calculados.


velocitycorrection::velocitycorrection()
{

}

velocitycorrection::~velocitycorrection()
{

}

velocitycorrection::velocitycorrection(velocitycorrection &v)
{
	this->pmin = v.pmin;
	this->k = k;
}

void velocitycorrection::UpdateVelolicty(vector<Particle2D>* particles, neighbor nei, double rho, double dt, double n0p, double radius, double num_D)
{
	vector<Point2D> DdvCorrector;

	vector<vector<int>> nei_in = nei.get();
	double AuxPmin = 0.0;
	int j = 0;
	//FILE* debug;

	//fopen_s(&debug, "dv.vems", "w");
	for (int i = 0; i < particles->size(); i++)
	{
		if (!(*particles)[i].is_fluid()){
			pmin.push_back(0.0);
			continue;
		}
			

		AuxPmin = (*particles)[i].get_pr();
		for (int a = 1; a <= nei_in[i][0]; a++)
		{
			j = nei_in[i][a];

			if ((*particles)[j].is_wall())
				continue;

			double p = (*particles)[j].get_pr();

			if (p < AuxPmin)
				AuxPmin = p;
		}
		pmin.push_back(AuxPmin);
	}

	//FILE*pmindebug;
	//fopen_s(&pmindebug, "debug_pmin.vems", "w");
	//for (int a = 0; a < pmin.size(); a++)
	//	fprintf(pmindebug, "%lf\n", pmin[a]);
	//fclose(pmindebug);

	Point2D Part_j, Part_i;
	Point2D new_point;
	double dx = 0.0, dy = 0.0, dist = 0.0, kernel = 0.0, ddv = 0.0;

	for (int i = 0; i < particles->size(); i++)
	{
		new_point.x = 0.0;
		new_point.y = 0.0;

		DdvCorrector.push_back(new_point);

		if (!(*particles)[i].is_fluid())
			continue;

		for (int a = 1; a <= nei_in[i][0]; a++)
		{
			j = nei_in[i][a];

			if ((*particles)[j].is_wall())
				continue;

			dx = 0.0;
			dy = 0.0;
			dist = 0.0;
			kernel = 0.0;
			ddv = 0.0;

			Part_j = (*particles)[j].get_po();
			Part_i = (*particles)[i].get_po();

			dx = Part_i.x - Part_j.x;
			dy = Part_i.y - Part_j.y;

			dist = sqrt(dx*dx + dy*dy);

			ddv = dt;
			//ddv = ddv*((*particles)[j].get_pr() - pmin[i]);
			//if (j > pmin.size() || j > particles->size()) 
			//{
				//ddv = ddv*((*particles)[j].get_pr() - pmin[i]);
				//printf("bla %d\n", a);
			//}
			//else{
				//ddv = ddv*(((*particles)[i].get_pr() + (*particles)[j].get_pr()) - (pmin[i] + pmin[j]));
				//printf("%d\n", a);
			//}

			ddv = ddv*(((*particles)[i].get_pr() + (*particles)[j].get_pr()) - (pmin[i] + pmin[j]));

			ddv = ddv / dist *k.weight(dist, radius);
			ddv = ddv / rho;
			ddv = ddv / n0p;
			ddv = ddv* num_D;

			DdvCorrector[i].x += ((ddv * dx) / dist);
			DdvCorrector[i].y += ((ddv * dy) / dist);
		}

		//fprintf_s(debug, "%lf %lf\n", DdvCorrector[i].x, DdvCorrector[i].y);
	}

	Point2D v_i;

	for (int i = 0; i < particles->size(); i++)
	{
		Part_i = (*particles)[i].get_po();
		v_i = (*particles)[i].get_v();

		if ((*particles)[i].is_fluid())
		{
			(*particles)[i].set_po(Part_i.x + DdvCorrector[i].x*dt, Part_i.y + DdvCorrector[i].y*dt);
			(*particles)[i].set_v(v_i.x + DdvCorrector[i].x, v_i.y + DdvCorrector[i].y);
		}
	}

	pmin.clear();
	DdvCorrector.clear();
	//fclose(debug);
}

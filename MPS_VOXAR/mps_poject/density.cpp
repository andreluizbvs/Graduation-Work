#include "density.h"
#include <iostream>

void density::cal_n(neighbor nei, vector<Particle2D>* particles,double r)
{
	//FILE*debug_kerenl;
	//fopen_s(&debug_kerenl, "kernel.vems", "w");
	vector<vector<int>> nei_in = nei.get();
	vector<vector<int>> n_fluid;
	int j = 0/*, kj = 0*/;
	double dx = 0.0, dy = 0.0, dist = 0.0, k = 0.0;
	double n = 0.0;
	Point2D Part_j, Part_i;
	for(int i=0;i<particles->size();i++)
	{
		n = 0.0;
		//for (j = 0; j<fluid; j++)n_fluid[j][i] = 0.0;
		for(int a=1;a<=nei_in[i][0];a++)
		{
			j = nei_in[i][a];
			//kj = type[j];
			
			Part_j = (*particles)[j].get_po();
			Part_i = (*particles)[i].get_po();

			dx = Part_j.x - Part_i.x;
			dy = Part_j.y - Part_i.y;

			dist = sqrt(dx*dx + dy*dy);

			k = this->k.weight(dist,r);
			//fprintf(debug_kerenl,"%.15lf\n", k);
			n += k;
		}
		(*particles)[i].set_n(n);
	}

	//FILE* debug;

	//fopen_s(&debug, "saida_n.vems", "w");

	//for (int a = 0; a < particles->size(); a++)
	//	fprintf(debug, "%.15lf\n", (*particles)[a].get_n());
	//fclose(debug);
	//fclose(debug_kerenl);
}

density::density()
{

}

density::~density()
{

}

density::density(const density &d)
{
	this->k = d.k;
}
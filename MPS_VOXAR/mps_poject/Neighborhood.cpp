#include "Neighborhood.h"

// calculo dos vizinhos, a estrutura usada para armazenar é matriz simples
// a organização na matriz é: primeiro elementos da linha "i" indica o numero de vizinhos da particula "i"
// o resto da linha "i" contem o indice do vizinho de "i"
// [n_vizinhos,v1,v2,v3,4,v5...v8]



neighbor::neighbor(const neighbor &n)
{
	this->neighbors = n.neighbors;
}

neighbor::neighbor(int cont)
{
	vector<int>aux;
	for(int a=0;a<cont;a++)
	{
		for(int b=0;b<60;b++) //@albvs 60, na verdade 
		{
			aux.push_back(0);
		}
		neighbors.push_back(aux);
		aux.clear();
	}
}

neighbor::~neighbor()
{

}

void neighbor::clear()
{
	for(int a=0;a<neighbors.size();a++)
		neighbors[a][0]=0;
}

vector<vector<int>> neighbor::get()
{
	return neighbors;
}

void neighbor::set(vector<Particle2D>* particles,double r2)
{
	Point2D Part_j, Part_i;
	double x, y;
	for(int i=0; i<particles->size(); i++)
	{
		for(int j=i+1; j<particles->size(); j++)
		{
			if((*particles)[j].is_wall() && (*particles)[i].is_wall())
				continue;

			Part_j = (*particles)[j].get_po();
			Part_i = (*particles)[i].get_po();

			x = Part_j.x - Part_i.x;
			y = Part_j.y - Part_i.y;

			if(((x*x)+(y*y)) <= r2)
			{
				//cout<<x2+y2<<" "<<r2<<endl;

				neighbors[i][0]++;
				neighbors[j][0]++;

				neighbors[i][neighbors[i][0]] = j;
				neighbors[j][neighbors[j][0]] = i;
			}
		}

	}

}


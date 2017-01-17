//#include "Neighborhood.h"
//
//// C�lculo dos vizinhos: A estrutura usada para armazenar � matriz simples.
//// A organiza��o na matriz �: Primeiro elemento da linha "i" indica o n�mero de vizinhos da part�cula "i".
//// O resto da linha "i" cont�m o �ndice do vizinho da part�cula "i"
//// [n_vizinhos,v1,v2,v3,v4,v5,...,v8]
//
//
//neighbor::neighbor(const neighbor &n)
//{
//	this->neighbors = n.neighbors;
//}
//
//neighbor::neighbor(int cont)
//{
//	neighbors = new int*[cont];
//	for (int i = 0; i < cont; ++i) neighbors[i] = new int[cont];
//	for (int g = 0; g < cont; g++)
//	{
//		for (int h = 0; h < cont; h++)
//		{
//			neighbors[g][h] = 0;
//		}
//	}
//}
//
//neighbor::~neighbor()
//{
//
//}
//
//void neighbor::clear(int nump)
//{
//	for (int g = 0; g < nump; g++)
//	{
//		for (int h = 0; h < nump; h++)
//		{
//			neighbors[g][h] = 0;
//		}
//	}
//}
//
//int** neighbor::get()
//{
//	return neighbors;
//}
//
//void neighbor::set(Particle2D* particles,double r2, int nump)
//{
//	Point2D Part_j, Part_i;
//	double x, y;
//	
//	for(int i=0; i<nump; i++)
//	{
//		for(int j=i+1; j<nump; j++)
//		{
//			if((particles)[j].is_wall() && (particles)[i].is_wall())
//				continue;
//
//			Part_j = (particles)[j].get_po();
//			Part_i = (particles)[i].get_po();
//
//			x = Part_j.x - Part_i.x;
//			y = Part_j.y - Part_i.y;
//
//			if(((x*x)+(y*y)) <= r2)
//			{
//				neighbors[i][0]++;
//				neighbors[j][0]++;
//
//				neighbors[i][neighbors[i][0]] = j;
//				neighbors[j][neighbors[j][0]] = i;
//			}
//		}
//	}
//}
//

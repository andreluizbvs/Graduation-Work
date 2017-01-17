//#include "density.h"
//
//void density::cal_n(neighbor nei, Particle2D* particles,double r, int nump)
//{
//	int j = 0;
//	double dx = 0.0, dy = 0.0, dist = 0.0, k = 0.0, n = 0.0;
//	Point2D Part_j, Part_i;
//	for (int i = 0; i < nump; i++)
//	{
//		n = 0.0;
//		
//		for(int a=1;a<=nei.get()[i][0];a++)
//		{
//			j = nei.get()[i][a];
//			
//			Part_j = particles[j].get_po();
//			Part_i = particles[i].get_po();
//
//			dx = Part_j.x - Part_i.x;
//			dy = Part_j.y - Part_i.y;
//
//			dist = sqrt(dx*dx + dy*dy);
//
//			k = this->k.weight(dist,r);
//
//			n += k;
//		}
//		particles[i].set_n(n);
//	}
//}
//
//density::density()
//{
//
//}
//
//density::~density()
//{
//
//}
//
//density::density(const density &d)
//{
//	this->k = d.k;
//}
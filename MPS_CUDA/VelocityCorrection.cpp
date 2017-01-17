//#include "velocityCorrection.h"
//
//// Correção das velocidades e das posições das partículas de acordo com os valores de pressão calculados.
//
//
//velocityCorrection::velocityCorrection()
//{
//
//}
//
//velocityCorrection::~velocityCorrection()
//{
//
//}
//
//velocityCorrection::velocityCorrection(velocityCorrection &v)
//{
//	this->pmin = v.pmin;
//	this->k = k;
//}
//
//void velocityCorrection::UpdateVelocity(Particle2D* particles, neighbor nei, double rho, double dt, double n0p, double radius, double num_D, int nump)
//{
//	double *pmin = new double[nump];
//	double AuxPmin = 0.0;
//	int j = 0;
//
//	Point2D* DdvCorrector = new Point2D[nump];
//
//	for (int i = 0; i < nump; i++)
//	{
//		if (!(particles)[i].is_fluid()){
//			pmin[i] = 0.0;
//			continue;
//		}
//			
//		AuxPmin = (particles)[i].get_pr();
//		for (int a = 1; a <= nei.get()[i][0]; a++)
//		{
//			j = nei.get()[i][a];
//
//			if ((particles)[j].is_wall())
//				continue;
//
//			double p = (particles)[j].get_pr();
//
//			if (p < AuxPmin)
//				AuxPmin = p;
//		}
//		pmin[i] = (AuxPmin);
//	}
//
//	Point2D Part_j, Part_i;
//	Point2D new_point;
//	double dx = 0.0, dy = 0.0, dist = 0.0, kernel = 0.0, ddv = 0.0;
//
//	for (int i = 0; i < nump; i++)
//	{
//		new_point.x = 0.0;
//		new_point.y = 0.0;
//
//		DdvCorrector[i] = (new_point);
//
//		if (!(particles)[i].is_fluid())
//			continue;
//
//		for (int a = 1; a <= nei.get()[i][0]; a++)
//		{
//			j = nei.get()[i][a];
//
//			if ((particles)[j].is_wall())
//				continue;
//
//			dx = 0.0;
//			dy = 0.0;
//			dist = 0.0;
//			kernel = 0.0;
//			ddv = 0.0;
//
//			Part_j = (particles)[j].get_po();
//			Part_i = (particles)[i].get_po();
//
//			dx = Part_i.x - Part_j.x;
//			dy = Part_i.y - Part_j.y;
//
//			dist = sqrt(dx*dx + dy*dy);
//
//			ddv = dt;
//
//			ddv = ddv*(((particles)[i].get_pr() + (particles)[j].get_pr()) - (pmin[i] + pmin[j]));
//
//			ddv = ddv / dist *k.weight(dist, radius);
//			ddv = ddv / rho;
//			ddv = ddv / n0p;
//			ddv = ddv* num_D;
//
//			DdvCorrector[i].x += ((ddv * dx) / dist);
//			DdvCorrector[i].y += ((ddv * dy) / dist);
//		}
//	}
//
//	Point2D v_i;
//
//	for (int i = 0; i < nump; i++)
//	{
//		Part_i = (particles)[i].get_po();
//		v_i = (particles)[i].get_v();
//
//		if ((particles)[i].is_fluid())
//		{
//			(particles)[i].set_po(Part_i.x + DdvCorrector[i].x*dt, Part_i.y + DdvCorrector[i].y*dt);
//			(particles)[i].set_v(v_i.x + DdvCorrector[i].x, v_i.y + DdvCorrector[i].y);
//		}
//	}
//	delete [] pmin;
//	delete [] DdvCorrector;
//}

#include "TimerControl.h"

// Controle de tempo de cada interação: O tempo varia de acordo com a função abaixo

__global__ void set_cal_dt(int offset, Particle3D* particles, double *dt, double DtRatio, double AveParticleDis, double MaxDt, int nump)
{
	double vmax = 0;
	double v = 0;
	double dt_limt = 0;
	//double dt = 0;

	for (int a = 0; a<nump; a++)
	{
		if ((particles)[a].is_fluid())
		{
			v = pow((particles)[a].get_v().x, 2);
			v += pow((particles)[a].get_v().y, 2);
			v += pow((particles)[a].get_v().z, 2);

			if (v>vmax)
				vmax = v;
		}
	}

	vmax = sqrt(vmax);

	dt_limt = dt[0]*INCREASE_DT;

	if (vmax == 0.0){
		dt[0] = dt_limt;
	}else{

		dt[0] = DtRatio*AveParticleDis / vmax;

		if (dt[0] > dt_limt)
			dt[0] = dt_limt;

		if (dt[0] > MaxDt)
			dt[0] = MaxDt;
	}
}

double TimerControl::cal_dt(Particle2D* particles, double OldDt, double DtRatio, double AveParticleDis, double MaxDt, int nump)
{
	double vmax = 0;
	double v = 0;
	double dt_limt = 0;
	double dt = 0;

	for (int a = 0; a<nump; a++)
	{
		if ((particles)[a].is_fluid())
		{
			v = pow((particles)[a].get_v().x, 2);
			v += pow((particles)[a].get_v().y, 2);

			if (v>vmax)
				vmax = v;
		}
	}

	vmax = sqrt(vmax);

	dt_limt = OldDt*INCREASE_DT;

	if (vmax == 0.0)
		return dt_limt;

	dt = DtRatio*AveParticleDis / vmax;

	if (dt>dt_limt)
		dt = dt_limt;

	if (dt>MaxDt)
		dt = MaxDt;

	return dt;
}

TimerControl::TimerControl()
{

}

TimerControl::~TimerControl()
{

}

TimerControl::TimerControl(const TimerControl &t)
{

}
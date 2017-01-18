#include "TimerControl.h"


//controle de tempo de cada interação, tempo varia de acordo com a função abaixo


double TimerControl::cal_dt(vector<Particle2D>* particles,double OldDt, double DtRatio, double AveParticleDis, double MaxDt)
{
	double vmax = 0;
	double v = 0;
	double dt_limt = 0;
	double dt = 0;


	for(int a=0; a<particles->size(); a++)
	{
		if((*particles)[a].is_fluid())
		{
			v = pow((*particles)[a].get_v().x,2);
			v += pow((*particles)[a].get_v().y,2);
			
			if(v>vmax)
				vmax = v;
		}
	}


	vmax = sqrt(vmax);

	dt_limt = OldDt*INCREASE_DT;

	if(vmax==0.0)
		return dt_limt;

	dt = DtRatio*AveParticleDis/vmax;

	if(dt>dt_limt)
		dt = dt_limt;

	if(dt>MaxDt)
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
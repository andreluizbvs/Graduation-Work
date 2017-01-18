#ifndef __TIMERCONTROL_H__
#define __TIMERCONTROL_H__

#include "Particles.h"
#include <vector>
#include <cmath>

#define  INCREASE_DT 1.1

using namespace std;

class TimerControl
{
public:
	double cal_dt(vector<Particle2D>* particles,double OldDt, double DtRatio, double AveParticleDis, double MaxDt);
	
	TimerControl();
	~TimerControl();
	TimerControl(const TimerControl &c);
};


#endif
#ifndef __TIMERCONTROL_H__
#define __TIMERCONTROL_H__

#include "Particles.h"
#include <vector>
#include <cmath>
#include "cuda_runtime.h"

#define  INCREASE_DT 1.1

using namespace std;

__global__ void set_cal_dt2D(int offset, Particle2D* particles, double *dt, double DtRatio, double AveParticleDis, double MaxDt, int nump);
__global__ void set_cal_dt(int offset, Particle3D* particles, double *dt, double DtRatio, double AveParticleDis, double MaxDt, int nump);

class TimerControl
{
public:
	double cal_dt(Particle2D* particles, double OldDt, double DtRatio, double AveParticleDis, double MaxDt, int nump);

	TimerControl();
	~TimerControl();
	TimerControl(const TimerControl &c);
};


#endif
#include "ReadWrite.h"
#include "Particles.h"
#include "Neighborhood.h"
#include "density.h"
#include "TimerControl.h"
#include "ParticlesActions.h"
#include "data_in.h"
#include "SystemResolution.h"
#include "velocityCorrection.h"

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_svd.h>

#include <iostream>

using namespace std;

int main()
{
	double n0p = 0.0;
	double n0pICCG = 0.0;

	data_in*input_data = new data_in();
	vector<Particle2D>*particles;
	ReadWrite FileControl;
	TimerControl ControlTime;
	ParticlesActions Actions;
	velocitycorrection VCorrections;

	vector<double> dp;

	SystemResolution *sys;

	density op;

	particles = new vector < Particle2D > ;

	/********************************************************/
	FILE *outPointsDamBreak;
	errno_t errorCode = fopen_s(&outPointsDamBreak, "outPointsMPS.txt", "w");
	double maximum = 0.0;
	int location = 0;
	/********************************************************/

	FileControl.ReadData("mps.data", input_data);

	FileControl.ReadParticles2D(particles, "mps.grd");

	cout << particles->size() << endl;

	neighbor nei(particles->size());
	neighbor neiICCG(particles->size());

	nei.set(particles, input_data->radius2);
	neiICCG.set(particles, input_data->radius_ICCG2);

	op.cal_n(neiICCG, particles, input_data->radius_ICCG);

	n0pICCG = (*particles)[0].get_n();

	op.cal_n(nei, particles, input_data->radius);

	n0p = (*particles)[0].get_n();

	input_data->dt = input_data->MaxDt;

	FileControl.WriteOut(particles, 0);

	sys = new SystemResolution(particles);

	for (int it = 1; it <= input_data->MaxIteration; it++)
	{
		cout << it << endl;

		input_data->dt = ControlTime.cal_dt(particles, input_data->dt, input_data->DtRatio, input_data->AveParticleDis, input_data->MaxDt);

		Actions.external_force(particles, input_data->dt, input_data->Gravity);
		Actions.mov_part(particles, input_data->dt);
		Actions.collision(particles, input_data->dt, nei, input_data->radius2lim, input_data->DensityFluid0);

		nei.clear();
		neiICCG.clear();

		nei.set(particles, input_data->radius2);
		neiICCG.set(particles, input_data->radius_ICCG2);

		op.cal_n(nei, particles, input_data->radius);

		sys->set_bcon(particles, n0p, input_data->dirichlet);
		//condição de borda

		//sys->set_B(particles, n0p, input_data->dt, 2, input_data->lambda, n0pICCG, input_data->DensityFluid0);
		sys->set_B_HS_ECS(particles, n0p, input_data->dt, 2, input_data->lambda, n0pICCG, input_data->DensityFluid0, nei, input_data->radius);
		//vetor b
		//sys->set_A(particles, neiICCG, n0pICCG, input_data->radius_ICCG, input_data->lambda, input_data->DensityFluid0, input_data->dt, 2);
		sys->set_A_HL(particles, neiICCG, n0pICCG, input_data->radius_ICCG, input_data->lambda, input_data->DensityFluid0, input_data->dt, 2);
		//matriz a

		sys->system_resolution(neiICCG, input_data->dt, input_data->MaxIterPress, input_data->MinIterPress, input_data->Convergence, particles->size());
		//resultado do sistema

		dp = sys->get_dp();

		for (int a = 0; a < particles->size(); a++)
			(*particles)[a].set_pr(dp[a]);

		VCorrections.UpdateVelolicty(particles, neiICCG, input_data->DensityFluid0, input_data->dt, n0p, input_data->radius, 2);

		FileControl.WriteOut(particles, it);

		/*************************************************************************************///ALBVS begin

		//for (int k = 1; k < particles->size(); k++)
		//{
		//	if ((particles->at(k).get_po().x > maximum) && (abs(particles->at(k).get_v().x) + abs(particles->at(k).get_v().y) != 0.0))
		//	{
		//		maximum = particles->at(k).get_po().x;
		//		location = k + 1;
		//	}
		//}

		//fprintf(outPointsDamBreak, "%lf %lf %lf\n", maximum, ((it-1)*0.001)*sqrt(2 * 9.8 / 0.14389344), maximum / (0.14389344)); // L = 0.14389344, g = 9.8, t = it*10^-3

		/************************************************************************************///ALBVS end

		sys->clear();
		system("pause");
	}
	return 0;
}
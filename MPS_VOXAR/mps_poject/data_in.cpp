#include "data_in.h"


// todos os dados de entrada da simulação do arquivo com extenção .data
// e mais alguns que tem valores constantes durante toda a simulação.

data_in::data_in()
{
	int NumberOfFluid = 0;
	int wall_type = 0;
	int wall_p_type = 0;

	double DensityFluid0 = 0;
	double DensityFluid1 = 0;
	double DensityWall_P = 0;
	double DensityWall = 0;

	double ExpansionFluid0 = 0;
	double ExpansionFluid1 = 0;
	double ExpansionWall_P = 0;
	double ExpansionWall = 0;

	double Gravity = 0;
	double MaxTime = 0;

	int MaxIteration = 0;
	
	double MaxDt = 0;
	double DtRatio = 0;
	
	double Convergence = 0;

	int MaxIterPress = 0;
 	int MinIterPress = 0;

	double AveParticleDis = 0;

	double rePND = 0;
	double reICCG = 0;
	double minlen = 0;
	double dirichlet = 0;
	double lam_rat = 0;

	double col_rat = 0;

	double radius = 0, radius_ICCG = 0;
	double radius2 = 0, radius_ICCG2 = 0;
	double radius2lim = 0;
	double ra = 0;

	double lambda = 0;

	double dt = 0;
	double currentTime = 0;
}

void data_in::set_after_write()
{
	radius = AveParticleDis * rePND;
	radius_ICCG = AveParticleDis * reICCG;

	radius2 = radius*radius;
	radius_ICCG2 = radius_ICCG*radius_ICCG;

	radius2lim = AveParticleDis*AveParticleDis*minlen*minlen;

	ra = AveParticleDis/1.7724539;

	lambda = lam_rat*(radius_ICCG2*radius_ICCG2 / 12.0 - radius_ICCG*ra*ra*ra / 3.0 + ra*ra*ra*ra / 4.0)
		/ (radius_ICCG2 / 2.0 - radius_ICCG*ra + ra*ra/2.0);
}


data_in::~data_in()
{

}

data_in::data_in(const data_in &d)
{
	*(this) = d;
}
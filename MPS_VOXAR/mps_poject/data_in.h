#ifndef __DATA_H__
#define __DATA_H__


class data_in
{
public:
	data_in();
	data_in(const data_in&d);
	~data_in();
	void set_after_write();

	int NumberOfFluid;
	int wall_type;
	int wall_p_type;

	double DensityFluid0;
	double DensityFluid1;
	double DensityWall_P;
	double DensityWall;

	double ExpansionFluid0;
	double ExpansionFluid1;
	double ExpansionWall_P;
	double ExpansionWall;

	double Gravity;
	double MaxTime;

	int MaxIteration;
	
	double MaxDt;
	double DtRatio;
	
	double Convergence;

	int MaxIterPress;
 	int MinIterPress;

	double AveParticleDis;

	double rePND;
	double reICCG;
	double minlen;
	double dirichlet;
	double lam_rat;

	double col_rat;

	double radius, radius_ICCG;
	double radius2, radius_ICCG2;
	double radius2lim;
	double ra;

	double lambda;

	double dt;
	double currentTime;
};

#endif
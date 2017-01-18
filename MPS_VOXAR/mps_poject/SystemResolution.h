#ifndef __SYSTEMRESOLOTION_H__
#define __SYSTEMRESOLOTION_H__

#include "Particles.h"
#include "Neighborhood.h"
#include "Kernel.h"
#include "Neighborhood.h"
#include <vector>

#include <vnl/vnl_matrix.h>
#include <vnl/algo/vnl_matrix_inverse.h>
#include <vnl/vnl_inverse.h>
#include <vnl/algo/vnl_svd.h>
#include <vnl/algo/vnl_cholesky.h>
#include <vnl/vnl_inverse.h>


using namespace std;

class SystemResolution
{
public:
//protected:
	vector<double> B;
	vector<int> Bcon;
	vector<double> dp;
	vector<vector<double>> A;
	//vector<vector<double>> Anew;
	vector<vector<double>> IC;

	int cont_size_matrix;

	kernel k;

//public:
	SystemResolution(vector<Particle2D>* particles);
	~SystemResolution();
	SystemResolution(const SystemResolution &s);

	void clear();

	void set_bcon(vector<Particle2D>* particles, double n0p, double dirichlet);
	void set_B(vector<Particle2D>* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho);
	void set_B_HS_ECS(vector<Particle2D>* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho, neighbor nei, double radius);
	void set_A(vector<Particle2D>* particles, neighbor nei, double n0pICCG, double radius, double lambda, double rho, double dt, int num_d);
	void set_A_HL(vector<Particle2D>* particles, neighbor neiICCG, double n0pICCG, double radius_ICCG, double lambda, double rho, double dt, int num_d);
	vector<int> GetBcon();
	vector<double> get_dp();
	
	void system_resolution(neighbor nei, double dt, int MaxIterPress, int MinIterPress, double Convergence,int allsize);

protected:
	void conf_dp(double *dp,int size);


};

#endif
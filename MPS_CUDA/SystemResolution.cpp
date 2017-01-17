#include "SystemResolution.h"
#include "gmm.h"
#include "gmm\gmm_condition_number.h"
#include "AuxiliarFunctions.h"

// Resolução do sistema (resolução da PPE)
// Matriz A, vetor B, vetor bcon (condição de contorno)

SystemResolution::SystemResolution(Particle2D* particles, int nump)
{
	Bcon = new int[nump];
	dp = new double[nump];
}


SystemResolution::SystemResolution(const SystemResolution &s)
{
	this->A = s.A;
	this->B = s.B;
	this->Bcon = s.Bcon;
}


void SystemResolution::clear(int nump)
{
	for (int a = 0; a < A.size(); a++){ 
		A[a].clear(); 
	}
	A.clear();
	//B.clear();
	//delete[] B;
}


void SystemResolution::set_bcon(Particle2D* particles, double n0p, double dirichlet, int nump)
{
	for (int a = 0; a < nump; a++)
	{
		if ((particles)[a].is_wall()){
			Bcon[a] = -1;
		}
		else{
			if ((particles)[a].get_n() / n0p < dirichlet){
				Bcon[a] = 1;
			}
			else{
				Bcon[a] = 0;
			}
		}
	}
}

//ALBVS begin
//void SystemResolution::set_B_HS_ECS(Particle2D* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho, int *nei, double radius, int nump, int *bcon)
//{
//	//vector<double> temp_b;
//	double *temp_b = new double[nump];
//	B = new double[nump];
//	int contB = 0;
//
//	double aux = 0.0;
//	int j, l;
//	double vx, vy, sum = 0.0;
//	double dx = 0.0, dy = 0.0, dist = 0.0;
//	Point2D Part_j, Part_i, Part_jv, Part_iv;
//	double ECS = 0.0;
//	double alpha = 0.0;
//	double beta = 0.0;
//	int contTemp = 0;
//	for (int i = 0; i < nump; i++)
//	{
//		if (i > 0){
//			alpha = ((particles)[i].get_n() - n0p) / n0p; // (n[i]-n0p)/n0p
//			beta = (dt / n0p)*(temp_b[i - 1] / (-1.0 / (n0p*dt))); // (dt/n0)*(Dn/Dt)
//			ECS = fabs(alpha)*(beta / dt) + fabs(beta)*(alpha / dt);
//		}
//		sum = 0.0;
//		aux = 0.0;
//		for (l = 1; l <= nei[(i*nump) + 0]; l++){
//			j = nei[(i*nump) + l];
//
//			if ((particles)[j].is_wall()){
//				continue;
//			}
//			Part_j = (particles)[j].get_po();
//			Part_i = (particles)[i].get_po();
//
//			dx = ((Part_i.x) - (Part_j.x));
//			dy = ((Part_i.y) - (Part_j.y));
//
//			dist = sqrt(dx*dx + dy*dy);
//
//			Part_jv = (particles)[j].get_v();
//			Part_iv = (particles)[i].get_v();
//
//			vx = ((Part_iv.x) - (Part_jv.x));
//			vy = ((Part_iv.y) - (Part_jv.y));
//
//			sum += (radius / (dist*dist*dist))*(dx*vx + dy*vy); /////////////// MPS-HS
//		}
//		aux = (-1.0 / (n0p*dt))*sum + ECS; ///////////// MPS-ECS
//		//temp_b.push_back(aux);
//		temp_b[contTemp] = aux;
//		contTemp++;
//		aux = 0.0;
//	}
//
//	//int contB = 0;
//	for (int a = 0; a < /*temp_b.size()*/contTemp; a++){
//		if (bcon[a] == 0){
//			//B.push_back(temp_b[a]);
//			B[contB] = temp_b[a];
//			/*if (a < (contTemp - 1))*/contB++;
//		}
//	}
//
//	//temp_b.clear();
//	delete[] temp_b;
//}

//void SystemResolution::set_A_HL(Particle2D* particles, int *neiICCG, double n0pICCG, double radius_ICCG, double lambda, double rho, double dt, int num_d, int nump, int *bcon/*, double *srcA*/)
//{
//	vector<double> temp(nump);
//	//double *temp = new double[nump];
//	vector<vector<double>> temp_a;
//	//double *temp_a = new double[nump*nump];
//
//	int cont_line = 0;
//	int j;
//	Point2D Part_j, Part_i;
//	double dx = 0.0, dy = 0.0, dist = 0.0, kernel = 0.0, val = 1.0;
//
//	for (int i = 0; i < nump; i++)
//	{
//		for (int a = 1; a <= neiICCG[(i*nump) + 0]; a++)
//		{
//			j = neiICCG[(i*nump) + a];
//
//
//			if (bcon[j] == -1)
//			{
//				temp[j] = 0.0;
//			}
//			else
//			{
//				dx = 0.0;
//				dy = 0.0;
//				dist = 0.0;
//				kernel = 0.0;
//				val = 1.0;
//
//				Part_j = (particles)[j].get_po();
//				Part_i = (particles)[i].get_po();
//
//				dx = Part_j.x - Part_i.x;
//				dy = Part_j.y - Part_i.y;
//
//				dist = sqrt(dx*dx + dy*dy);
//
//				val = 3.0 * radius_ICCG;
//				val = val / (dist*dist*dist);
//				val = val / n0pICCG;
//				val = val / rho;
//
//				temp[j] = -val;
//				temp[cont_line] += val;
//			}
//		}
//
//		temp_a.push_back(temp);
//		cont_line++;
//
//		for (int v = 0; v < temp.size(); v++)
//			temp[v] = 0;
//	}
//
//	temp.clear();
//
//	for (int a = 0; a < nump; a++)
//	{
//		if (bcon[a] != 0)
//			continue;
//
//		for (int b = 0; b < nump; b++)
//		{
//			if (bcon[b] != 0)
//				continue;
//
//
//			temp.push_back(temp_a[a][b]);
//		}
//
//		A.push_back(temp);
//
//		temp.clear();
//
//		temp_a[a].clear();
//
//	}
//	temp.clear();
//	temp_a.clear();
//}
//ALBVS end

void SystemResolution::system_resolution(int MaxIterPress, int MinIterPress, double Convergence, int allsize, int *bcon, int contB, double *srcB, double *srcA)
{
	//size_t nbdof = B.size();
	size_t nbdof = contB;

	std::vector<double> Xgmm(nbdof), Bgmm(nbdof);

	gmm::row_matrix< gmm::rsvector<double> > Mgmm(nbdof, nbdof);
	
	//for (int a = 0; a < contB*contB; a++){
	//	//Mgmm[a][b] = A[a][b];
	//	if (a < contB)Bgmm[a] = srcB[a];
	//	Mgmm[a/contB][a%contB] = srcA[a/contB + a%contB];
	//}
	for (int a = 0; a < contB; a++){
		Bgmm[a] = srcB[a];
		for (int b = 0; b < contB; b++){
			Mgmm[a][b] = srcA[a*contB + b];
		}
	}

	gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > P(Mgmm, 10, 1e-6);

	gmm::iteration iter(1E-10);

	gmm::csc_matrix<double> M2;

	gmm::clean(Mgmm, 1E-16);

	//gmm::copy(Mgmm, M2);
	//printf("resolvendo\n");
	gmm::bicgstab(Mgmm, Xgmm, Bgmm, P, iter);

	double *dpArray = (double*)calloc(contB/*B.size()*/, sizeof(double));

	for (int a = 0; a < contB/*B.size()*/; a++)
		dpArray[a] = Xgmm[a];

	conf_dp(dpArray, allsize, bcon, contB);
	free(dpArray);

	Xgmm.clear();
	Bgmm.clear();
	Mgmm.clear_mat();
}


void SystemResolution::conf_dp(double* dpArray, int size, int *bcon, int contB)
{
	int desloc = 0;
	int nbdof = size;
	for (int a = 0; a < size; a++)
	{
		if (a - desloc < contB)
		{
			if (bcon[a] == 0)
			{
				if (dpArray[a - desloc]>0)
					dp[a] = dpArray[a - desloc];
				else
					dp[a] = 0;
			}
			else
			{
				dp[a] = 0;
				desloc++;
			}
		}
		else
			dp[a] = 0;
	}
}

//double* SystemResolution::get_dp()
//{
//	return this->dp;
//}

int* SystemResolution::GetBcon()
{
	return Bcon;
}


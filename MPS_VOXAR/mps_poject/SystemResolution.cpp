#include "SystemResolution.h"
#include "AuxiliarFunctions.h"

#include "gmm.h"
#include "gmm\gmm_condition_number.h"


//resolução de sistema
//matriz a, vetor b,  vetor bcon (condição de contorno)

SystemResolution::SystemResolution(vector<Particle2D>* particles)
{

	int n = particles->size();

	for (int a = 0; a < n; a++)
	{
		Bcon.push_back(0);
		dp.push_back(0.0);
	}

}



SystemResolution::SystemResolution(const SystemResolution &s)
{
	this->A = s.A;
	this->B = s.B;
	this->Bcon = s.Bcon;
}



void SystemResolution::clear()
{
	for (int a = 0; a < IC.size(); a++)
		IC[a].clear();

	IC.clear();
	B.clear();

	for (int a = 0; a < Bcon.size(); a++)
		Bcon[a] = 0;

	for (int a = 0; a < A.size(); a++)
		A[a].clear();
	A.clear();
}



void SystemResolution::set_bcon(vector<Particle2D>* particles, double n0p, double dirichlet)
{
	for (int a = 0; a < particles->size(); a++)
	{
		if ((*particles)[a].is_wall()){
			Bcon[a] = -1;
		}
		else{
			if ((*particles)[a].get_n() / n0p < dirichlet){
				Bcon[a] = 1;
			}
			/*albvs*/else{
				Bcon[a] = 0;
			}
		}
	}
}



void SystemResolution::set_B(vector<Particle2D>* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho)
{
	vector<double> temp_b;

	double aux = 0.0;
	double dif = 0.0;

	for (int a = 0; a < particles->size(); a++)
	{
		dif = ((*particles)[a].get_n() - n0p);
		aux = 1.0 / dt;
		aux = aux / dt;
		aux = aux * dif;
		aux = aux / n0p;

		/******************************/

		//aux = aux / (2 * num_d);
		//aux = aux * lambda;
		//aux = aux * n0pICCG;
		//aux = aux * rho;

		temp_b.push_back(aux);
		aux = 0.0;
	}

	//system("pause");

	for (int a = 0; a < temp_b.size(); a++)
	{
		if (Bcon[a] == 0)
			B.push_back(temp_b[a]);
	}

	temp_b.clear();

	//FILE*debug;
	//fopen_s(&debug, "saida_B.vems", "w");
	//for (int a = 0; a < B.size(); a++)
	//{
	//	fprintf(debug, "%.15lf\n", B[a]);
	//}

	//fclose(debug);
}


//ALBVS begin
void SystemResolution::set_B_HS_ECS(vector<Particle2D>* particles, double n0p, double dt, int num_d, double lambda, double n0pICCG, double rho, neighbor nei, double radius)
{
	vector<vector<int>> nei_in = nei.get();
	vector<double> temp_b;


	double aux = 0.0;
	int j, l;
	double vx, vy, sum = 0.0;
	double dx = 0.0, dy = 0.0, dist = 0.0;
	Point2D Part_j, Part_i, Part_jv, Part_iv;
	double ECS = 0.0;
	double alpha = 0.0;
	double beta = 0.0;

	for (int i = 0; i < particles->size(); i++)
	{
		if(i > 0){
			alpha = ((*particles)[i].get_n() - n0p) / n0p; // (n[i]-n0p)/n0p
			beta = (dt / n0p)*(temp_b[i-1] / (-1.0 / (n0p*dt))); // (dt/n0)*(Dn/Dt)
			ECS = fabs(alpha)*(beta / dt) + fabs(beta)*(alpha / dt);
		}
		sum = 0.0;
		aux = 0.0;
		for (l = 1; l <= nei_in[i][0]; l++){
			j = nei_in[i][l];
			
			if ((*particles)[j].is_wall()){ 
				//sum = 0.0;
				continue; 
			}
			Part_j = (*particles)[j].get_po();
			Part_i = (*particles)[i].get_po();

			dx = ((Part_i.x) - (Part_j.x));
			dy = ((Part_i.y) - (Part_j.y));

			dist = sqrt(dx*dx + dy*dy);

			Part_jv = (*particles)[j].get_v();
			Part_iv = (*particles)[i].get_v();

			vx = ((Part_iv.x) - (Part_jv.x));
			vy = ((Part_iv.y) - (Part_jv.y));
			//if (vx != 0.0 || vy != 0.0) {
			//	cout << vx << " " << vy << endl;
			//	//system("pause");
			//}
			sum += (radius / (dist*dist*dist))*(dx*vx + dy*vy); /////////////// MPS-HS
		}
		aux = (-1.0 / (n0p*dt))*sum + ECS; ///////////// MPS-ECS
		temp_b.push_back(aux);
		aux = 0.0;
	}

	//system("pause");

	for (int a = 0; a < temp_b.size(); a++)
	{
		if (Bcon[a] == 0)
			B.push_back(temp_b[a]);
	}

	//ofstream vectorB;
	//vectorB.open("vectorB.txt");
	//for (int i = 0; i < B.size(); i++) {
	//	vectorB << B[i] << endl;
	//}

	temp_b.clear();

	//FILE*debug;
	//fopen_s(&debug, "saida_B.vems", "w");
	//for (int a = 0; a < B.size(); a++)
	//{
	//	fprintf(debug, "%.15lf\n", B[a]);
	//}

	//fclose(debug);
}
//ALBVS end

void SystemResolution::set_A(vector<Particle2D>* particles, neighbor nei, double n0pICCG, double radius, double lambda, double rho, double dt, int num_d)
{
	vector<vector<int>> nei_in = nei.get();

	vector<vector<double>> temp_a;
	int cont_line = 0;

	vector<double> temp(particles->size());

	for (int i = 0; i < particles->size(); i++)
	{
		if (Bcon[i] != 0)return;

		for (int a = 1; a <= nei_in[i][0]; a++)
		{
			int j = nei_in[i][a];


			if (Bcon[j] == -1)
			{
				temp[j] = 0.0;
			}
			else
			{
				double dx = 0.0, dy = 0.0, dist = 0.0, kernel = 0.0, val = 1.0;
				Point2D Part_j, Part_i;

				Part_j = (*particles)[j].get_po();
				Part_i = (*particles)[i].get_po();

				dx = Part_j.x - Part_i.x;
				dy = Part_j.y - Part_i.y;

				dist = sqrt(dx*dx + dy*dy);

				kernel = this->k.weight(dist, radius);

				val = 2.0 * num_d;
				val = val / lambda;
				val = val*kernel;
				val = val / n0pICCG;
				val = val / rho;

				temp[j] = -val;
				temp[cont_line] += val;
			}
		}

		temp_a.push_back(temp);
		cont_line++;

		for (int v = 0; v < temp.size(); v++)
			temp[v] = 0.0;
	}

	temp.clear();

	for (int a = 0; a < particles->size(); a++)
	{
		if (Bcon[a] != 0)
			continue;

		for (int b = 0; b < particles->size(); b++)
		{
			if (Bcon[b] != 0)
				continue;


			temp.push_back(temp_a[a][b]);
		}

		A.push_back(temp);

		temp.clear();

		temp_a[a].clear();
	}

	//FILE*debug_;
	//fopen_s(&debug_, "matriz2.vems", "w");

	//for (int a = 0; a < Anew.size(); a++)
	//{
	//	for (int b = 0; b < Anew.size(); b++)
	//	{
	//		fprintf(debug_, "%.15lf ", Anew[a][b]);
	//	}
	//	fprintf(debug_, "\n");
	//}
	//fclose(debug_);
	//system("pause");


	temp.clear();
	temp_a.clear();
}

//ALBVS begin
void SystemResolution::set_A_HL(vector<Particle2D>* particles, neighbor neiICCG, double n0pICCG, double radius_ICCG, double lambda, double rho, double dt, int num_d)
{
	vector<vector<int>> nei_in = neiICCG.get();

	vector<vector<double>> temp_a;
	int cont_line = 0;
	int j;
	Point2D Part_j, Part_i;
	double dx = 0.0, dy = 0.0, dist = 0.0, kernel = 0.0, val = 1.0;

	vector<double> temp(particles->size());

	for (int i = 0; i < particles->size(); i++)
	{
		for (int a = 1; a <= nei_in[i][0]; a++)
		{
			j = nei_in[i][a];


			if (Bcon[j] == -1)
			{
				temp[j] = 0.0;
			}
			else
			{
				dx = 0.0;
				dy = 0.0;
				dist = 0.0;
				kernel = 0.0;
				val = 1.0;
				
				Part_j = (*particles)[j].get_po();
				Part_i = (*particles)[i].get_po();

				dx = Part_j.x - Part_i.x;
				dy = Part_j.y - Part_i.y;

				dist = sqrt(dx*dx + dy*dy);

				kernel = this->k.weight(dist, radius_ICCG);

				val = 3.0 * radius_ICCG;
				val = val / (dist*dist*dist);
				val = val / n0pICCG;
				val = val / rho;

				temp[j] = -val;
				temp[cont_line] += val;
			}
		}

		temp_a.push_back(temp);
		cont_line++;

		for (int v = 0; v < temp.size(); v++)
			temp[v] = 0;
	}

	temp.clear();

	for (int a = 0; a < particles->size(); a++)
	{
		if (Bcon[a] != 0)
			continue;

		for (int b = 0; b < particles->size(); b++)
		{
			if (Bcon[b] != 0)
				continue;


			temp.push_back(temp_a[a][b]);
		}

		A.push_back(temp);

		temp.clear();

		temp_a[a].clear();
	}

	//FILE*debug_;
	//fopen_s(&debug_, "matriz2.vems", "w");

	//for (int a = 0; a < Anew.size(); a++)
	//{
	//	for (int b = 0; b < Anew.size(); b++)
	//	{
	//		fprintf(debug_, "%.15lf ", Anew[a][b]);
	//	}
	//	fprintf(debug_, "\n");
	//}
	//fclose(debug_);
	//system("pause");


	temp.clear();
	temp_a.clear();
}
//ALBVS end

void SystemResolution::system_resolution(neighbor nei, double dt, int MaxIterPress, int MinIterPress, double Convergence, int allsize){

	int nbdof = B.size();

	std::vector<double> Xgmm(nbdof), Bgmm(nbdof);

	gmm::row_matrix< gmm::rsvector<double> > Mgmm(nbdof, nbdof);

	for (int a = 0; a < B.size(); a++)
	{
		Bgmm[a] = B[a];

		for (int b = 0; b < B.size(); b++)
		{
			Mgmm[a][b] = A[a][b];
		}
	}

	gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > P(Mgmm, 10, 1e-6);

	gmm::iteration iter(1E-8);

	gmm::csc_matrix<double> M2;

	gmm::clean(Mgmm, 1E-12);

	gmm::copy(Mgmm, M2);
	printf("resolvendo\n");
	gmm::bicgstab(M2, Xgmm, Bgmm, P, iter);
	

	//FILE* debug_;
	//fopen_s(&debug_, "p_puro.vems", "w");

	//for (int a = 0; a < Xgmm.size(); a++)
	//{
	//	fprintf(debug_, "%lf\n", Xgmm[a]);
	//}
	//fclose(debug_);
	//system("pause");

	double *dpArray = (double*)calloc(B.size(), sizeof(double));

	for (int a = 0; a < B.size(); a++)
		dpArray[a] = Xgmm[a];

	conf_dp(dpArray,allsize);
	free(dpArray);

	Xgmm.clear();
	Bgmm.clear();

	Mgmm.clear_mat();
}


void SystemResolution::conf_dp(double* dpArray, int size)
{
	int desloc=0;
	for (int a = 0; a < size; a++)
	{
		if (a - desloc < B.size())
		{
			if (Bcon[a] == 0)
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

vector<double> SystemResolution::get_dp()
{
	return this->dp;
}

vector <int> SystemResolution::GetBcon()
{
	return Bcon;
}


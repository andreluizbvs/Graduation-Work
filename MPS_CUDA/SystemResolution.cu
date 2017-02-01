#include "SystemResolution.h"

__global__ void conf_dp_kernel(double* dp, double* dpArray, int size, int *bcon, int *contB)
{
	int desloc = 0;
	int nbdof = size;
	for (int a = 0; a < size; a++)
	{
		if (a - desloc < contB[0])
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

__global__ void set_bcon_kernel(int offset, Particle3D* particles, double n0p, double dirichlet, int nump, int *bcon)
{
	int a = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	if (particles[a].is_wall()){
		bcon[a] = -1;
	}
	else{
		if ((particles[a].get_n() / n0p) < dirichlet){
			bcon[a] = 1;
		}
		else{
			bcon[a] = 0;
		}
	}
}

__global__ void cal_B_HS_ECS_kernel(int offset, Particle3D* particles, double n0p, double *dt, int *nei, double radius, int nump, int *bcon, int *contTemp, double *temp_b)
{
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	int j, l;
	double vx, vy, vz, sum = 0.0;
	double dx = 0.0, dy = 0.0, dz = 0.0, dist = 0.0;
	Point3D Part_j, Part_i, Part_jv, Part_iv;
	double ECS = 0.0;
	double alpha = 0.0;
	double beta = 0.0;

	if (i > 0){
		alpha = ((particles)[i].get_n() - n0p) / n0p; // (n[i]-n0p)/n0p
		beta = (dt[0] / n0p)*(temp_b[i - 1] / (-1.0 / (n0p*dt[0]))); // (dt/n0)*(Dn/Dt)
		ECS = fabs(alpha)*(beta / dt[0]) + fabs(beta)*(alpha / dt[0]);
	}
	sum = 0.0;

	for (l = 1; l <= nei[(i*nump) + 0]; l++){
		j = nei[(i*nump) + l];

		if ((particles)[j].is_wall()){
			continue;
		}
		Part_j = (particles)[i].get_po();
		Part_i = (particles)[j].get_po();

		dx = ((Part_i.x) - (Part_j.x));
		dy = ((Part_i.y) - (Part_j.y));
		dz = ((Part_i.z) - (Part_j.z));

		dist = sqrt(dx*dx + dy*dy + dz*dz);

		Part_jv = (particles)[i].get_v();
		Part_iv = (particles)[j].get_v();

		vx = ((Part_iv.x) - (Part_jv.x));
		vy = ((Part_iv.y) - (Part_jv.y));
		vz = ((Part_iv.z) - (Part_jv.z));

		sum += (radius / (dist*dist*dist))*(dx*vx + dy*vy + dz*vz); /////////////// MPS-HS
	}

	temp_b[i] = (-1.25 / (n0p*dt[0]))*sum + ECS; ///////////// MPS-ECS
	if (!(particles)[i].is_wall()){
		atomicAdd(&contTemp[0], 1);
	}


}

__global__ void cal_B_kernel(int offset, Particle3D* particles, double n0p, double *dt, int *contTemp, double *temp_b)
{
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);

	if (!(particles)[i].is_wall()){
		temp_b[i] = 1.0 / dt[0] / dt[0] * ((particles)[i].get_n() - n0p) / n0p;
		atomicAdd(&contTemp[0], 1);
	}

}

__global__ void set_B_HS_ECS_kernel(int offset, int *bcon, int *contB, double *srcB, double *temp_b, int *contTemp){
	
	//unsigned int a = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	for (int a = 0; a < contTemp[0]; a++){
		if (bcon[a] == 0){
			srcB[contB[0]] = temp_b[a];
			atomicAdd(&contB[0], 1);
		}
	}
}

__global__ void cal_A_HL_kernel(int offset, Particle3D* particles, int *neiICCG, double n0pICCG, double radius_ICCG, double rho, int num_d, int nump, int *bcon, double *temp_a)
{
	unsigned int i = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	int j;
	Point3D Part_j, Part_i;
	double dx = 0.0, dy = 0.0, dz = 0.0, dist = 0.0, val = 1.0;

	for (int a = 1; a <= neiICCG[(i*nump) + 0]; a++)
	{
		j = neiICCG[(i*nump) + a];

		if (bcon[j] == -1){
			temp_a[nump*i + j] = 0.0;
		}
		else{
			dx = 0.0;
			dy = 0.0;
			dz = 0.0;
			dist = 0.0;
			val = 1.0;

			Part_j = (particles)[j].get_po();
			Part_i = (particles)[i].get_po();

			dx = Part_j.x - Part_i.x;
			dy = Part_j.y - Part_i.y;
			dz = Part_j.z - Part_i.z;

			dist = sqrt(dx*dx + dy*dy + dz*dz);

			if (num_d == 2){
				val = 3.0 * radius_ICCG;
			}
			else{
				val = 2.0 * radius_ICCG;
			}
			val = val / (dist*dist*dist);
			val = val / n0pICCG;
			val = val / rho;

			temp_a[nump*i + j] = -val;
			temp_a[nump*i + i] += val;
		}
	}
}

__global__ void set_A_HL_kernel1(int offset, int *bcon, double *temp_a2, double *temp_a, int nump, int *bcon_cont){

	unsigned int a = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	int contRows = 0, contCols = 0;

	//for (int a = 0; a < nump; a++){
		//if (bcon[a] != 0){continue;}
	if (bcon[a] != 0){}else{
		for (int b = 0; b < nump; b++){
			if (bcon[b] != 0){}else{
				//temp_a[nump*contRows + contCols] = temp_a[a*nump + b];
				temp_a2[nump*(a - bcon_cont[a]) + (b - bcon_cont[b])] = temp_a[a*nump + b];
				//contCols++;
			}
		}
		//contRows++;
		//contCols = 0;
		//}
	}
}

__global__ void set_A_HL_kernel2(int offsetX, int offsetY, double *srcA, int *contB, double *temp_a, int nump){
	
	//unsigned int g = offset + (blockDim.x * blockIdx.x + threadIdx.x);
	unsigned int g = offsetX*1024 + blockIdx.x;
	//unsigned int h = offset + (blockDim.y * blockIdx.y + threadIdx.y);
	unsigned int h = offsetY*1024 + threadIdx.x;
	int idxOrig = (nump * g) + h;
	int idxDest = (contB[0] * g) + h;

	srcA[idxDest] = temp_a[idxOrig];
	
}

void SystemResolution::set_A_HL(Particle2D* particles, int *neiICCG, double n0pICCG, double radius_ICCG, double lambda, double rho, double dt, int num_d, int nump, int *bcon)
{
	//vector<double> temp(nump);
	double *temp = new double[nump];
	memset(temp, 0, sizeof(double)*nump);
	//vector<vector<double>> temp_a;
	double *temp_a = new double[nump*nump];
	memset(temp_a, 0, sizeof(double)*nump*nump);
	double *srcA = new double[nump*nump];
	memset(srcA, 0, sizeof(double)*nump*nump);

	int cont_line = 0;
	int j;
	Point2D Part_j, Part_i;
	double dx = 0.0, dy = 0.0, dist = 0.0, val = 1.0;

	for (int i = 0; i < nump; i++)
	{
		for (int a = 1; a <= neiICCG[(i*nump) + 0]; a++)
		{
			j = neiICCG[(i*nump) + a];


			if (bcon[j] == -1)
			{
				temp[j] = 0.0;
			}
			else
			{
				dx = 0.0;
				dy = 0.0;
				dist = 0.0;
				//kernel = 0.0;
				val = 1.0;

				Part_j = (particles)[j].get_po();
				Part_i = (particles)[i].get_po();

				dx = Part_j.x - Part_i.x;
				dy = Part_j.y - Part_i.y;

				dist = sqrt(dx*dx + dy*dy);
				
				if (num_d == 2){
					val = 3.0 * radius_ICCG;
				}
				else{
					val = 2.0 * radius_ICCG;
				}
				val = val / (dist*dist*dist);
				val = val / n0pICCG;
				val = val / rho;

				temp[j] = -val;
				temp[cont_line] += val;
			}
		}

		for (int g = 0; g < nump; g++){
			temp_a[nump*i + g] = temp[g];
		}

		//temp_a.push_back(temp);
		cont_line++;

		/*for (int v = 0; v < nump; v++)
			temp[v] = 0;*/
		memset(temp, 0, sizeof(double)*nump);
	}

	//temp.clear();
	//for (int v = 0; v < nump; v++)
	//	temp[v] = 0;
	memset(temp, 0, sizeof(double)*nump);

	int contRows = 0, contCols = 0;

	for (int a = 0; a < nump; a++){
		if (bcon[a] == 0){
			for (int b = 0; b < nump; b++){
				if (bcon[b] == 0){
					srcA[nump*contRows + contCols] = temp_a[a*nump + b];
					contCols++;
				}
			}
			contRows++;
			contCols = 0;
		}		
	}

	delete [] temp;
	delete [] temp_a;
	//temp.clear();
	//temp_a.clear();
}
//ALBVS end
//
//__host__ __device__ void SystemResolution::system_resolution(double dt, int MaxIterPress, int MinIterPress, double Convergence, int allsize){
//
//	//int nbdof = sizeof(B) / sizeof(B[0]);
//	//printf("a: %d, b: %d\n", sizeof(B) , sizeof(B[0]));
//	size_t nbdof = B.size();
//
//	std::vector<double> Xgmm(nbdof), Bgmm(nbdof);
//
//	gmm::row_matrix< gmm::rsvector<double> > Mgmm(nbdof, nbdof);
//
//	Bgmm = B;
//	//Mgmm = A;
//
//	for (int a = 0; a < B.size(); a++)
//	{
//		//Bgmm[a] = B[a];
//
//		for (int b = 0; b < B.size(); b++)
//		{
//			Mgmm[a][b] = A[a][b];
//			//printf("a: %d, b: %d\n", a, b);
//		}
//	}
//
//	gmm::ilut_precond< gmm::row_matrix< gmm::rsvector<double> > > P(Mgmm, 10, 1e-6);
//
//	gmm::iteration iter(1E-10);
//
//	gmm::csc_matrix<double> M2;
//
//	gmm::clean(Mgmm, 1E-16);
//
//	//gmm::copy(Mgmm, M2);
//	printf("resolvendo\n");
//	gmm::bicgstab(Mgmm, Xgmm, Bgmm, P, iter);
//
//
//	//FILE* debug_;
//	//fopen_s(&debug_, "p_puro.vems", "w");
//
//	//for (int a = 0; a < Xgmm.size(); a++)
//	//{
//	//	fprintf(debug_, "%lf\n", Xgmm[a]);
//	//}
//	//fclose(debug_);
//	//system("pause");
//
//	double *dpArray = (double*)calloc(B.size(), sizeof(double));
//
//	for (int a = 0; a < B.size(); a++)
//		dpArray[a] = Xgmm[a];
//
//	conf_dp(dpArray, allsize);
//	free(dpArray);
//
//	Xgmm.clear();
//	Bgmm.clear();
//
//	Mgmm.clear_mat();
//}
//
//
//__host__ __device__ void SystemResolution::conf_dp(double* dpArray, int size)
//{
//	int desloc = 0;
//	//int nbdof = sizeof(B) / sizeof(B[0]);
//	int nbdof = size;
//	for (int a = 0; a < size; a++)
//	{
//		if (a - desloc < B.size())
//		{
//			if (Bcon[a] == 0)
//			{
//				if (dpArray[a - desloc]>0)
//					dp[a] = dpArray[a - desloc];
//				else
//					dp[a] = 0;
//			}
//			else
//			{
//				dp[a] = 0;
//				desloc++;
//			}
//		}
//		else
//			dp[a] = 0;
//	}
//}
//
__host__ __device__ double* SystemResolution::get_dp()
{
	return this->dp;
}
//
//__host__ __device__ int* SystemResolution::GetBcon()
//{
//	return Bcon;
//}
//

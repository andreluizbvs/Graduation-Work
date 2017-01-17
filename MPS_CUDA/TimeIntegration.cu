#include "TimeIntegration.h"

#include <cusp\monitor.h>
#include <cusp\linear_operator.h>
#include <cusp\csr_matrix.h>
#include <cusp\krylov\bicgstab.h>
#include <cusp\gallery\poisson.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <cuda_runtime.h>
#include <cusparse.h>
#include <cusparse_v2.h>
#include <thrust\device_ptr.h>
#include <cuda.h>

using namespace std;

static const char *_cusparseGetErrorEnum(cusparseStatus_t error)
{
	switch (error)
	{

	case CUSPARSE_STATUS_SUCCESS:
		return "CUSPARSE_STATUS_SUCCESS";

	case CUSPARSE_STATUS_NOT_INITIALIZED:
		return "CUSPARSE_STATUS_NOT_INITIALIZED";

	case CUSPARSE_STATUS_ALLOC_FAILED:
		return "CUSPARSE_STATUS_ALLOC_FAILED";

	case CUSPARSE_STATUS_INVALID_VALUE:
		return "CUSPARSE_STATUS_INVALID_VALUE";

	case CUSPARSE_STATUS_ARCH_MISMATCH:
		return "CUSPARSE_STATUS_ARCH_MISMATCH";

	case CUSPARSE_STATUS_MAPPING_ERROR:
		return "CUSPARSE_STATUS_MAPPING_ERROR";

	case CUSPARSE_STATUS_EXECUTION_FAILED:
		return "CUSPARSE_STATUS_EXECUTION_FAILED";

	case CUSPARSE_STATUS_INTERNAL_ERROR:
		return "CUSPARSE_STATUS_INTERNAL_ERROR";

	case CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED:
		return "CUSPARSE_STATUS_MATRIX_TYPE_NOT_SUPPORTED";

	case CUSPARSE_STATUS_ZERO_PIVOT:
		return "CUSPARSE_STATUS_ZERO_PIVOT";
	}

	return "<unknown>";
}

inline void __cusparseSafeCall(cusparseStatus_t err, const char *file, const int line)
{
	if (CUSPARSE_STATUS_SUCCESS != err) {
		fprintf(stderr, "CUSPARSE error in file '%s', line %Ndims\Nobjs %s\nerror %Ndims: %s\nterminating!\Nobjs", __FILE__, __LINE__, err, \
			_cusparseGetErrorEnum(err)); \
			cudaDeviceReset(); assert(0); \
	}
}

extern "C" void cusparseSafeCall(cusparseStatus_t err) { __cusparseSafeCall(err, __FILE__, __LINE__); }

void TimeIntegration(Particle3D *particles, ReadWrite FileControl, data_in *input_data, int nump)
{
	int z, cont_bcon = 0;
	int max_nei = nump;						// Número máximo de vizinhos
	FILE *outNei = fopen("outNei.txt", "a");
	FILE *outTimes = fopen("outTimes.txt", "w");

	// Sistema
	//SystemResolution *sys;
	// --- Inicializando cuSPARSE e criando descritor da matriz esparsa A
	cusparseHandle_t handle;    cusparseSafeCall(cusparseCreate(&handle));
	cusparseMatDescr_t descrA;      cusparseSafeCall(cusparseCreateMatDescr(&descrA));
	cusparseSetMatType(descrA, CUSPARSE_MATRIX_TYPE_GENERAL);
	cusparseSetMatIndexBase(descrA, CUSPARSE_INDEX_BASE_ZERO);

	int division = nump / 1024;
	int mod = nump - (division * 1024);
	int DIM = 2;

	// Array de partículas e alocado no device
	Particle3D *particles_d = NULL;
	cudaError_t err = cudaSuccess;
	err = cudaMalloc((void **)&particles_d, nump*sizeof(Particle3D));
	
	/********************************************************/// Informação sobre a partícula mais a direita e tempo da simulação (Crista da onda x andamento da simulação)
	/*FILE *outPointsDamBreak;
	errno_t errorCode = fopen_s(&outPointsDamBreak, "outPointsMPS.txt", "w");
	double maximum = 0.0;
	int location = 0;*/
	/********************************************************/

	// Inicializando e setando os arrays de vizinhos para a iteração 0
	/*neighbor nei(nump, max_nei);
	neighbor neiICCG(nump, max_nei);

	nei.set(particles, input_data->radius2, nump);
	int *nei_1D = new int[nump*max_nei];
	for (int a = 0; a < nump; a++){
		for (int b = 0; b < max_nei; b++)	{
			nei_1D[(a * nump) + b] = nei.neighbors[a][b];
		}
	}

	neiICCG.set(particles, input_data->radius_ICCG2, nump);
	int *neiICCG_1D = new int[nump*max_nei];;
	for (int a = 0; a < nump; a++){
		for (int b = 0; b < max_nei; b++)	{
			neiICCG_1D[(a * nump) + b] = neiICCG.neighbors[a][b];
		}
	}
	*/
	int *nei_d = NULL;
	int *neiICCG_d = NULL;

	err = cudaMalloc((void **)&nei_d, nump*max_nei*sizeof(int));
	err = cudaMalloc((void **)&neiICCG_d, nump*max_nei*sizeof(int));
	cudaMemset(nei_d, 0, sizeof(int) * nump * max_nei);
	cudaMemset(neiICCG_d, 0, sizeof(int) * nump * max_nei);
	/*err = cudaMemcpy(nei_d, nei_1D, nump*max_nei*sizeof(int), cudaMemcpyHostToDevice);
	err = cudaMemcpy(neiICCG_d, neiICCG_1D, nump*max_nei*sizeof(int), cudaMemcpyHostToDevice);*/
	
	// Array de partículas transferido do host para o device
	err = cudaMemcpy(particles_d, particles, nump*sizeof(Particle3D), cudaMemcpyHostToDevice);
	//sys = new SystemResolution(particles, nump);

	// Inicializando e setando os arrays de vizinhos para a iteração 0
	set_nei_kernel << <division, 1024 >> >(0, particles_d, input_data->radius2, nump, nei_d);
	set_nei_kernel << <1, mod >> >(division * 1024, particles_d, input_data->radius2, nump, nei_d);
	set_nei_kernel << <division, 1024 >> >(0, particles_d, input_data->radius_ICCG2, nump, neiICCG_d);
	set_nei_kernel << <1, mod >> >(division * 1024, particles_d, input_data->radius_ICCG2, nump, neiICCG_d);

	// Setando os valores de n (particle number density) de cada partícula e também os valores de n0 & n0_iccg para a iteração 0
	cal_n_kernel << <division, 1024 >> >(0, neiICCG_d, particles_d, input_data->radius_ICCG, nump);
	cal_n_kernel << <1, mod >> >(division * 1024, neiICCG_d, particles_d, input_data->radius_ICCG, nump);
	cudaDeviceSynchronize();

	double n0pICCG = 0.0, *n0pICCG_d = NULL, *temp = new double[1];
	cudaMalloc((void **)&n0pICCG_d, sizeof(double));
	set_n_kernel << <1, 1 >> >(particles_d, n0pICCG_d);
	cudaDeviceSynchronize();
	cudaMemcpy(temp, n0pICCG_d, sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(n0pICCG_d);
	n0pICCG = temp[0];
	cout << n0pICCG << endl;
	system("pause");

	cal_n_kernel << <division, 1024 >> >(0, nei_d, particles_d, input_data->radius, nump);
	cal_n_kernel << <1, mod >> >(division * 1024, nei_d, particles_d, input_data->radius, nump);
	cudaDeviceSynchronize();

	double n0p = 0.0, *n0p_d = NULL;
	cudaMalloc((void **)&n0p_d, sizeof(double));
	set_n_kernel << <1, 1 >> >(particles_d, n0p_d);
	cudaDeviceSynchronize();
	cudaMemcpy(temp, n0p_d, sizeof(double), cudaMemcpyDeviceToHost);
	cudaFree(n0p_d);
	n0p = temp[0];
	cout << n0p << endl;
	system("pause");

	// Setando valor de dt
	double *dt_d = NULL;
	temp[0] = input_data->dt;
	cudaMalloc((void **)&dt_d, sizeof(double));
	cudaMemcpy(dt_d, temp, sizeof(double), cudaMemcpyHostToDevice);

	// Criando array pmin, temp_b, temp_a2, dpArray e DdvCorrector no DEVICE & bcon, contB, contTemp, srcB, temp_a, dp e srcA no HOST e no DEVICE
	double *pmin_d = NULL;
	cudaMalloc((void **)&pmin_d, nump*sizeof(double));
	Point3D *DdvCorrector_d = NULL;
	cudaMalloc((void **)&DdvCorrector_d, nump*sizeof(Point3D));
	double *temp_b_d = NULL;
	cudaMalloc((void **)&temp_b_d, nump*sizeof(double));
	int *bcon = new int[nump];
	int *bcon_d = NULL;
	cudaMalloc((void **)&bcon_d, nump*sizeof(int));
	int *bcon_cont = new int[nump];
	int *bcon_cont_d = NULL;
	cudaMalloc((void **)&bcon_cont_d, nump*sizeof(int));
	int *contB = new int[1];
	int *contB_d = NULL;
	cudaMalloc((void **)&contB_d, sizeof(int));
	cudaMemset(contB_d, 0, sizeof(int));
	int *contTemp = new int[1];
	int *contTemp_d = NULL;
	cudaMalloc((void **)&contTemp_d, sizeof(int));
	cudaMemset(contTemp_d, 0, sizeof(int));
	double *srcB = new double[nump];
	double *srcB_d = NULL;
	cudaMalloc((void **)&srcB_d, nump*sizeof(double));
	//double *temp_a = new double[nump*nump];
	double *temp_a_d = NULL;
	cudaMalloc((void **)&temp_a_d, nump*nump*sizeof(double));
	double *temp_a2_d = NULL;
	cudaMalloc((void **)&temp_a2_d, nump*nump*sizeof(double));
	double *srcA = new double[nump*nump];
	double *srcA_d = NULL;
	cudaMalloc((void **)&srcA_d, nump*nump*sizeof(double));
	double *dp = new double[nump];
	double *dp_d = NULL;
	err = cudaMalloc((void **)&dp_d, nump*sizeof(double));
	double *dpArray_d = NULL;
	err = cudaMalloc((void **)&dpArray_d, nump*sizeof(double));

	/*size_t avail;
	size_t total;
	
	double used_db, free_db, total_db;

	cudaEvent_t *t_begin, *t_end;
	int functionsNum = 20;
	t_begin = new cudaEvent_t[functionsNum];
	t_end = new cudaEvent_t[functionsNum];
	float *time_spent = new float[functionsNum];
	time_spent[19] = 0.0;

	for (int r = 0; r < functionsNum; r++) {
		cudaEventCreate(&t_begin[r]);
		cudaEventCreate(&t_end[r]);
	}*/
	
	for (int it = 1; it <= input_data->MaxIteration; it++)
	{
		//cout << it << endl;

		//cudaEventRecord(t_begin[0]);
		// Calculando o valor de dt (passo de tempo)
		//set_cal_dt <<< 1, 1 >> > (0, particles_d, dt_d, input_data->DtRatio, input_data->AveParticleDis, input_data->MaxDt, nump);
		//cudaEventRecord(t_end[0]);

		//cudaEventRecord(t_begin[1]);
		// Calculando a influência das forças externas nas partícula, seus movimentos e colisão entre elas 
		external_force_kernel << <division, 1024 >> >(0, particles_d, dt_d, input_data->Gravity);
		external_force_kernel << <1, mod >> >(division * 1024, particles_d, dt_d, input_data->Gravity);
		//cudaEventRecord(t_end[1]);
	
		//cudaEventRecord(t_begin[2]);
		mov_part_kernel << <division, 1024 >> >(0, particles_d, dt_d);
		mov_part_kernel << <1, mod >> >(division * 1024, particles_d, dt_d);
		//cudaEventRecord(t_end[2]);

		//cudaEventRecord(t_begin[3]);
		collision_kernel << <division, 1024 >> >(0, particles_d, dt_d, nei_d, input_data->radius2lim, input_data->DensityFluid0, nump);
		collision_kernel << <1, mod >> >(division * 1024, particles_d, dt_d, nei_d, input_data->radius2lim, input_data->DensityFluid0, nump);
		//cudaEventRecord(t_end[3]);
		cudaDeviceSynchronize();

		// Limpando e recalculando o array de vizinhos
		cudaMemset(nei_d, 0, sizeof(int) * nump * max_nei);
		cudaMemset(neiICCG_d, 0, sizeof(int) * nump * max_nei);

		//cudaEventRecord(t_begin[4]);
		set_nei_kernel << <division, 1024 >> >(0, particles_d, input_data->radius2, nump, nei_d);
		set_nei_kernel << <1, mod >> >(division * 1024, particles_d, input_data->radius2, nump, nei_d);
		set_nei_kernel << <division, 1024 >> >(0, particles_d, input_data->radius_ICCG2, nump, neiICCG_d);
		set_nei_kernel << <1, mod >> >(division * 1024, particles_d, input_data->radius_ICCG2, nump, neiICCG_d);
		//cudaEventRecord(t_end[4]);

		//cudaEventRecord(t_begin[5]);
		// Recalculando n (particle number density)
		cal_n_kernel << <division, 1024 >> >(0, nei_d, particles_d, input_data->radius, nump);
		cal_n_kernel << <1, mod >> >(division * 1024, nei_d, particles_d, input_data->radius, nump);
		//cudaEventRecord(t_end[5]);

		//cudaEventRecord(t_begin[6]);
		// Setando condição de contorno
		set_bcon_kernel << <division, 1024 >> >(0, particles_d, n0p, input_data->dirichlet, nump, bcon_d);
		set_bcon_kernel << <1, mod >> >(division * 1024, particles_d, n0p, input_data->dirichlet, nump, bcon_d);
		//cudaEventRecord(t_end[6]);
		cudaDeviceSynchronize();
		
		// Vetor B (source term)
		cudaMemset(temp_b_d, 0, sizeof(double) * nump);
		cudaMemset(contTemp_d, 0, sizeof(int));

		//cudaEventRecord(t_begin[7]);
		cal_B_HS_ECS_kernel << <division, 1024 >> >(0, particles_d, n0p, dt_d, nei_d, input_data->radius, nump, bcon_d, contTemp_d, temp_b_d);
		cal_B_HS_ECS_kernel << <1, mod >> >(division * 1024, particles_d, n0p, dt_d, nei_d, input_data->radius, nump, bcon_d, contTemp_d, temp_b_d);
		//cudaEventRecord(t_end[7]);
		cudaDeviceSynchronize();

		/*cal_B_kernel << <division, 1024 >> >(0, particles_d, n0p, dt_d, contTemp_d, temp_b_d);
		cal_B_kernel << <1, mod >> >(division * 1024, particles_d, n0p, dt_d, contTemp_d, temp_b_d);
		cudaDeviceSynchronize();*/


		cudaMemset(srcB_d, 0, sizeof(double) * nump);
		cudaMemset(contB_d, 0, sizeof(int));

		//cudaEventRecord(t_begin[8]);
		set_B_HS_ECS_kernel << <1, 1 >> >(0, bcon_d, contB_d, srcB_d, temp_b_d, contTemp_d);
		//cudaEventRecord(t_end[8]);
		cudaDeviceSynchronize();

		// Matriz A (Laplaciano)
		cudaMemset(temp_a_d, 0, sizeof(double) * nump*nump);
		//cudaEventRecord(t_begin[9]);
		cal_A_HL_kernel << <division, 1024 >> >(0, particles_d, neiICCG_d, n0pICCG, input_data->radius_ICCG, input_data->DensityFluid0, DIM, nump, bcon_d, temp_a_d);
		cal_A_HL_kernel << <1, mod >> >(division * 1024, particles_d, neiICCG_d, n0pICCG, input_data->radius_ICCG, input_data->DensityFluid0, DIM, nump, bcon_d, temp_a_d);
		//cudaEventRecord(t_end[9]);
		cudaDeviceSynchronize();

		/*err = cudaMemcpy(temp_a, temp_a_d, nump*nump*sizeof(double), cudaMemcpyDeviceToHost);
		for (int g = 0; g < nump*nump; g++){
			if (g > 0 && g%nump == 0) fprintf(outNei, "\n");
			fprintf(outNei,"%lf ", temp_a[g]);
		}*/

		err = cudaMemcpy(contB, contB_d, sizeof(int), cudaMemcpyDeviceToHost);
		err = cudaMemcpy(bcon, bcon_d, nump*sizeof(int), cudaMemcpyDeviceToHost);

		cont_bcon = 0;
		for (int i = 0; i < nump; i++){
			if (bcon[i] != 0){
				bcon_cont[i] = cont_bcon;
				cont_bcon++;
			}
			else{
				bcon_cont[i] = cont_bcon;
			}
		}
		
		cudaMemset(srcA_d, 0, sizeof(double) * nump*nump);
		err = cudaMemcpy(bcon_cont_d, bcon_cont, nump*sizeof(int), cudaMemcpyHostToDevice);

		//cudaEventRecord(t_begin[10]);
		set_A_HL_kernel1 << <division, 1024 >> >(0, bcon_d, temp_a2_d, temp_a_d, nump, bcon_cont_d);
		set_A_HL_kernel1 << <1, mod >> >(division * 1024, bcon_d, temp_a2_d, temp_a_d, nump, bcon_cont_d);
		//cudaEventRecord(t_end[10]);
		cudaDeviceSynchronize();

		/*err = cudaMemcpy(temp_a, temp_a_d, nump*nump*sizeof(double), cudaMemcpyDeviceToHost);
		for (int g = 0; g < nump*nump; g++){
			if (g > 0 && g%nump == 0) fprintf(outNei, "\n");
			fprintf(outNei, "%lf ", temp_a[g]);
		}*/

		if ((contB[0]/1024) < 1){ 
			//cudaEventRecord(t_begin[11]);
			set_A_HL_kernel2 << <contB[0], contB[0] >> >(0, 0, srcA_d, contB_d, temp_a2_d, nump);
			//cudaEventRecord(t_end[11]);
		}
		else{			
			//cudaEventRecord(t_begin[11]);
			//for (int y = 0; y < contB[0] / 1024; y++){
				for (z = 0; z < contB[0] / 1024; z++){
					//set_A_HL_kernel2 << < 1024, 1024 >> >(y, z, srcA_d, contB_d, temp_a2_d, nump);
					set_A_HL_kernel2 << < contB[0], 1024 >> >(0, z, srcA_d, contB_d, temp_a2_d, nump);
					cudaDeviceSynchronize();
				}
			//}
			set_A_HL_kernel2 << <contB[0], (contB[0] % 1024) >> >(0, z, srcA_d, contB_d, temp_a2_d, nump);
			//cudaEventRecord(t_end[11]);
			
			//for (int it = 0; it < contB[0] / 1024; it++){
				//set_A_HL_kernel2 << <1024, (contB[0] % 1024) >> >(it, z + 1, srcA_d, contB_d, temp_a2_d, nump);
				//set_A_HL_kernel2 << <(contB[0] % 1024), 1024 >> >(y + 1, it, srcA_d, contB_d, temp_a2_d, nump);
			//}			
			//set_A_HL_kernel2 << <(contB[0] % 1024), (contB[0] % 1024) >> >(y+1, z+1, srcA_d, contB_d, temp_a2_d, nump);
		}
		cudaDeviceSynchronize();


		/*for (int g = 0; g < nump*nump; g++){
			if (g > 0 && g%nump == 0) fprintf(outNei, "\n");
			fprintf(outNei,"%lf ", srcA[g]);
		}*/

		// Resolução do sistema linear

		/*if (nump > 10000){
			cudaMemcpy(srcA, srcA_d, nump*nump*sizeof(double), cudaMemcpyDeviceToHost);
			cudaMemcpy(srcB, srcB_d, nump*sizeof(double), cudaMemcpyDeviceToHost);
			sys->system_resolution(input_data->MaxIterPress, input_data->MinIterPress, input_data->Convergence, nump, bcon, contB[0], srcB, srcA);
			dp = sys->get_dp();
			cudaMemcpy(dp_d, dp, nump*sizeof(double), cudaMemcpyHostToDevice);
		}
		else{
*/
			/*************************************************************************************************************************************************/
			int nnz = 0;                                // Onde é guardado o número de elementos diferentes de zero na matriz A em sua forma densa
			const int lda = contB[0];                   // Dimensão principal da matriz A
			// Descobrindo o número de elementos diferentes de zero e o n.e.d.z. por linha na matriz A (em sua forma densa)
			int *d_nnzPerVector;
			cudaMalloc(&d_nnzPerVector, contB[0] * sizeof(*d_nnzPerVector));
			/*cusparseSafeCall*/(cusparseDnnz(handle, CUSPARSE_DIRECTION_ROW, contB[0], contB[0], descrA, srcA_d, lda, d_nnzPerVector, &nnz));

			// Arrays para armazenar a matriz A no formato CSR (compressed sparse row)
			double *d_A;            (cudaMalloc(&d_A, nnz * sizeof(*d_A)));
			int *d_A_RowIndices;    (cudaMalloc(&d_A_RowIndices, (contB[0] + 1) * sizeof(*d_A_RowIndices)));
			int *d_A_ColIndices;    (cudaMalloc(&d_A_ColIndices, nnz * sizeof(*d_A_ColIndices)));

			// Ponteiros para o device indicando o início de cada um dos array
			thrust::device_ptr<double> dev_ptr_d_A = thrust::device_pointer_cast(d_A);
			thrust::device_ptr<int> dev_ptr_d_A_RowIndices = thrust::device_pointer_cast(d_A_RowIndices);
			thrust::device_ptr<int> dev_ptr_d_A_ColIndices = thrust::device_pointer_cast(d_A_ColIndices);
			thrust::device_ptr<double> dev_ptr_srcB_d = thrust::device_pointer_cast(srcB_d);
			//thrust::device_ptr<int> dev_ptr_dpArray = thrust::device_pointer_cast(dpArray_d);

			// Convertendo de fato a matrix densa em CSR
			/*cusparseSafeCall*/(cusparseDdense2csr(handle, contB[0], contB[0], descrA, srcA_d, lda, d_nnzPerVector, d_A, d_A_RowIndices, d_A_ColIndices));

			cusp::csr_matrix<int, double, cusp::device_memory> csrA;
			csrA.resize(contB[0], contB[0], nnz);
			csrA.row_offsets.assign(dev_ptr_d_A_RowIndices, dev_ptr_d_A_RowIndices + contB[0] + 1);
			csrA.column_indices.assign(dev_ptr_d_A_ColIndices, dev_ptr_d_A_ColIndices + nnz);
			csrA.values.assign(dev_ptr_d_A, dev_ptr_d_A + nnz);

			cusp::array1d<double, cusp::device_memory> x(csrA.num_rows, 0);
			cusp::array1d<double, cusp::device_memory> array1dB(dev_ptr_srcB_d, dev_ptr_srcB_d + contB[0]);

			// Configurando critério de parada:
			//  Limite de iterações = 100
			//  Tolerância relativa = 1e-16
			//  Tolerância absoluta = 0
			//  Verbose            = true
			cusp::monitor<double> monitor(array1dB, 100, 1e-16, 0, false);

			// Configurando precondicionador (identidade)
			cusp::identity_operator<double, cusp::device_memory> M(csrA.num_rows, csrA.num_rows);

			// Resolvendo o sistema linear A x = b
			//cudaEventRecord(t_begin[12]);
			cusp::krylov::bicgstab(csrA, x, array1dB, monitor, M);
			//cudaEventRecord(t_end[12]);

			// Armazenando o resultado (pressões das partículas)
			//cusp::array1d<double, cusp::host_memory> array1d_dp(x);
			//for (int i = 0; i < contB[0]; i++) dp[i] = array1d_dp[i];
			dpArray_d = thrust::raw_pointer_cast(&x[0]);

			// Desalocando os arrays e matriz do device
			{
				cusp::csr_matrix<int, double, cusp::device_memory> tmp;
				csrA.swap(tmp);
				cusp::array1d<double, cusp::device_memory> temp;
				x.swap(temp);
				array1dB.swap(temp);
			}
			cudaFree(d_nnzPerVector);
			cudaFree(d_A);
			cudaFree(d_A_ColIndices);
			cudaFree(d_A_RowIndices);

			/*************************************************************************************************************************************************/

			// Setando pressões das partículas 
			//err = cudaMemcpy(dpArray_d, dp, nump*sizeof(double), cudaMemcpyHostToDevice);

			//cudaEventRecord(t_begin[13]);
			conf_dp_kernel << <1, 1 >> >(dp_d, dpArray_d, nump, bcon_d, contB_d);
			//cudaEventRecord(t_end[13]);
		//}

		//cudaEventRecord(t_begin[14]);
		set_pr_kernel << <division, 1024>> >(0, particles_d, dp_d);
		set_pr_kernel << <1, mod >> >(division * 1024, particles_d, dp_d);
		//cudaEventRecord(t_end[14]);

		//cudaEventRecord(t_begin[15]);
		// Cálculo do gradiente da velocidade -> atualização das velocidades
		set_pmin_kernel << <division, 1024 >> >(0,  particles_d, neiICCG_d, nump, pmin_d);
		set_pmin_kernel << <1, mod >> >(division * 1024, particles_d, neiICCG_d, nump, pmin_d);
		//cudaEventRecord(t_end[15]);

		//cudaEventRecord(t_begin[16]);
		set_ddv_kernel << <division, 1024 >> >(0, particles_d, neiICCG_d, input_data->DensityFluid0, dt_d, n0p, input_data->radius, DIM, nump, pmin_d, DdvCorrector_d);
		set_ddv_kernel << <1, mod >> >(division * 1024, particles_d, neiICCG_d, input_data->DensityFluid0, dt_d, n0p, input_data->radius, DIM, nump, pmin_d, DdvCorrector_d);
		//cudaEventRecord(t_end[16]);

		////cudaEventRecord(t_begin[17]);
		updateVelocity_kernel << <division, 1024 >> >(0, particles_d, dt_d, nump, DdvCorrector_d);
		updateVelocity_kernel << <1, mod >> >(division * 1024, particles_d, dt_d, nump, DdvCorrector_d);
		//cudaEventRecord(t_end[17]);
		cudaDeviceSynchronize();

		/*time_spent[18] = 0.0;
		for (int r = 0; r <= 17; r++) {
			cudaEventElapsedTime(&time_spent[r], t_begin[r], t_end[r]);
			fprintf(outTimes, "Função #%d: %lf\n", r, time_spent[r]);
			time_spent[18] = time_spent[18] + time_spent[r];
		}
		fprintf(outTimes, "Tempo total decorrido devido às funções CUDA em milissegundos: %lf\n", time_spent[18]);
		time_spent[19] = time_spent[19] + time_spent[18];*/

		// Movendo informações das partículas do device para o host para impressão nos VTUs
		err = cudaMemcpy(particles, particles_d, nump*sizeof(Particle3D), cudaMemcpyDeviceToHost);

		// Escrevendo nos VTUs a configuação do sistema no time step atual
		FileControl.WriteOut(particles, it, nump);

		// Informação sobre a partícula mais a direita e tempo da simulação (Crista da onda x andamento da simulação)
		/*************************************************************************************///ALBVS begin
		//for (int k = 1; k < nump; k++){
		//	if ((particles[k].get_po().x > maximum) && (abs(particles[k].get_v().x) + abs(particles[k].get_v().y) != 0.0)){
		//		maximum = particles[k].get_po().x;
		//		location = k + 1;
		//	}
		//}
		//fprintf(outPointsDamBreak, "%lf %lf\n"/*, (maximum - 0.02)*/, ((it-1)*0.00025)*sqrt(9.8 / 0.6), (maximum - 0.02) / 0.6); // L = 0.6 m, g = 9.8 m/s², t = it*10^-3 s
		/************************************************************************************///ALBVS end
		//sys->clear(nump);
	}
	//fprintf(outTimes, "Tempo total da simulação devido às funções CUDA em milissegundos: %lf\n", time_spent[19]);

	/*cudaMemGetInfo(&avail, &total);
	free_db = (double)avail;

	total_db = (double)total;

	used_db = total_db - free_db;

	printf("GPU memory usage: used = %lf MB, free = %lf MB, total = %lf MB\n",

	used_db / 1024.0 / 1024.0, free_db / 1024.0 / 1024.0, total_db / 1024.0 / 1024.0);

	system("pause");*/

	cudaFree(pmin_d);
	cudaFree(bcon_d);
	cudaFree(bcon_cont_d);
	cudaFree(dp_d);
	cudaFree(dpArray_d);
	cudaFree(particles_d);
	cudaFree(dt_d);
	cudaFree(bcon_d);
	cudaFree(nei_d);
	cudaFree(neiICCG_d);
	cudaFree(srcB_d);
	cudaFree(temp_b_d);
	cudaFree(contB_d);
	cudaFree(contTemp_d);
	cudaFree(srcA_d);
	cudaFree(temp_a_d);
	cudaFree(temp_a2_d);
	delete[] dp;
	/*delete[] nei_1D;
	delete[] neiICCG_1D;*/
	delete[] bcon;
	delete[] bcon_cont;
	delete[] srcB;
	//delete[] temp_a;
	delete[] contB;
	delete[] contTemp;
	delete[] srcA;
	fclose(outNei);
	fclose(outTimes);
	//fclose(outPointsDamBreak);
}
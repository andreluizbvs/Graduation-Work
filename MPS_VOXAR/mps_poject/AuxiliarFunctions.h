
void set_vector_zero(int nump, double x[3000])
{
	int i;

	for (i = 0; i<nump; i++)x[i] = 0.0;
}



void mul_matrix_vector(int nump, int*bcon, int **list, double *y, double **mat, double *x)
{
	int i, j, l;

	for (i = 0; i<nump; i++)
	{
		if (bcon[i] != 0) continue;
		y[i] = mat[i][0] * x[i];
		for (l = 1; l <= list[i][0]; l++)
		{
			j = list[i][l];
			if (bcon[j] == -1)continue;
			y[i] = y[i] + mat[i][l] * x[j];
		}
	}
}


void sub_vectors(int nump, double *y, double *x1, double *x2, double a)
{
	int i;

	for (i = 0; i<nump; i++)y[i] = x1[i] - a*x2[i];
}


void copy_vectors(int nump, double *x1, double *x2)
{
	int i;

	for (i = 0; i<nump; i++) x2[i] = x1[i];
}

void add_vectors(int nump, double *y, double *x1, double *x2, double a)
{
	int i;

	for (i = 0; i<nump; i++)y[i] = x1[i] + a*x2[i];
}


void mul_vectors(int nump, int*bcon, double*ans, double *x1, double *x2)
{
	int i;

	*ans = 0.0;
	for (i = 0; i<nump; i++)
	{
		if (bcon[i] == 0) *ans += x1[i] * x2[i];
	}
}

#include <stdio.h>
#include <math.h>
#include <string.h>
#include <assert.h>
#include <stdlib.h>
#include <cblas.h>

#include "diis.h"



int diis_step(double **err, double **v, int steps, int dim, double *v_res) 
{
	double *B;
	double *rhs;
	int info;
	int i;
	int j;
	int B_dim;

	int nrhs = 1;
	int *piv;

	B_dim = steps + 1;

	B = (double *)malloc (B_dim * B_dim * sizeof(double));
	rhs = (double *)malloc (B_dim * sizeof(double));
	piv = (int *)malloc (B_dim * sizeof(int));
	memset (B, 0, B_dim * B_dim * sizeof(double));
	memset (rhs, 0, B_dim * sizeof(double));
	memset (piv, 0, B_dim * sizeof(int));
	
	rhs[B_dim - 1] = 1.0;
	
	for (i = 0; i < B_dim - 1; i++) {
		for (j = 0; j < B_dim - 1; j++) {
			B[i + (B_dim * j)] = cblas_ddot (dim, err[i], 1, err[j], 1);
		}
	}

	for (i = 0; i < B_dim - 1; i++) {
		B[i + B_dim * (B_dim - 1)] = 1.0;
		B[B_dim - 1 + B_dim * i] = 1.0;
	}
	

	dgetrf_ (&B_dim, &B_dim, B, piv, &info);
	
	if (info != 0) return -1;

	dgetrs_ ("N", &nrhs, B, &B_dim, piv, rhs, &B_dim, &info);  

	if (info != 0) return -1;

	/*Now perform the extrapolation*/

	for (i = 0; i < steps; i++) {
	        for (j = 0; j < dim; j++) {
			v_res[j] += rhs[i] * v[i][j];
		}
	}


	free (piv);
	free (B);
	free (rhs);
	return 1;
}
	

	

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
	double test = 0.0;

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
	

	dgetrf_ (&B_dim, &B_dim, B, &B_dim, piv, &info);

	
	if (info != 0) {

		free (piv);
		free (B);
		free (rhs);
		return -1;
	}
	dgetrs_ ("N", &B_dim, &nrhs, B, &B_dim, piv, rhs, &B_dim, &info);  

	if (info != 0) {
		
		free (piv);
		free (B);
		free (rhs);
		return -1;
	}
	/*Now perform the extrapolation*/

	for (i = 0; i < steps; i++) {
		test += rhs[i];
	}
	fprintf(stderr, "\n test: %lf", test);


	for (i = 0; i < steps; i++) {
	        for (j = 0; j < dim; j++) {
			v_res[j] += rhs[i] * v[i][j];
		}
	}


	free (piv);
	free (B);
	free (rhs);
	return 1;


	/* free (piv); */
	/* free (B); */
	/* free (rhs); */
	/* return 1; */
}


void diis_extrap (diis_t *inp, double *result)
{
	int steps;
	int info;
	double *tmp; 
	
	tmp = (double *)malloc (inp->errs->vec_dim * sizeof(double));
	memset (tmp, 0, inp->errs->vec_dim * sizeof(double));
	
	if (inp->errs->num < inp->diis_lim) {
		steps = inp->errs->num;
	} else {
		steps = inp->diis_lim;
	}
	

	do {
		info = diis_step (inp->errs->vecs, inp->vecs->vecs, steps, inp->vecs->vec_dim, tmp);

		if (info == -1) {
			(inp->fails)++;
			(inp->errs->num)--;
			(inp->vecs->num)--;
			steps--;
		}
		if (info == 1) {
			memcpy (result, tmp, inp->vecs->vec_dim * sizeof(double));
			free (tmp);
			return;
			
		}
		
		} while (steps > 0);
		
	if (steps == 0) {
		memcpy (result, inp->vecs->vecs[0], inp->vecs->vec_dim * sizeof(double));
		clear_vec_stack (inp->errs);
		clear_vec_stack (inp->vecs);
		(inp->restarts)++;
		free (tmp);
		return;
	}
	
	free (tmp);
	return;


}

diis_t *init_diis(int diis_lim, int vec_len)
{
	diis_t *inp;
	
	inp = (diis_t *)malloc (sizeof(diis_t));
	inp->fails = 0;
	inp->restarts = 0;
	inp->diis_lim = diis_lim;
	inp->errs = create_vec_stack (vec_len, diis_lim);
	inp->vecs = create_vec_stack (vec_len, diis_lim);
	return inp;
	

}

void destroy_diis(diis_t *inp)
{
	destroy_vec_stack (inp->errs);
	destroy_vec_stack (inp->vecs);
	free (inp);

}


#if 0	
void diis_extrap (double **err, double **v, int dim, int diis_lim, int iter, double *v_res) 
{

	int steps;
	int info;
	double *tmp;

	tmp = (double *)malloc (dim * sizeof(double));
	memset (tmp, 0, dim * sizeof(double));

	if (iter < diis_lim) {
		steps = iter;
	} else {
		steps = diis_lim;
	}
	
	info = diis_step (err, v, steps, dim, tmp);

	if (info == -1) { 
		fprintf (stderr, "\n DIIS FAIL \n");
		free (tmp);
		return;
	} else if (info = 1) {
		memcpy (v_res, tmp, dim * sizeof(double));
		free (tmp);
		return;
	}


	return;
}
#endif	

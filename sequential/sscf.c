#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
/* #include <mkl.h> */
/* #include <mkl_lapack.h> */
#include <cblas.h>
#include <clapack.h>
/* #include <lapacke.h> */

#include "common.h"
#include "utils.h"

#include "erd_integral.h"

#include "one_electron.h"
#define TOLER 1.0e-5
#define MAX_IT 10


int sscf(double *H, double *S_sinv, int n, int n_ele, int maxit, double *D_old, double *D_new, double *F);
void compute_D(int n, int n_ele, double *F, double *D);
void print_evals(double *A, int n);

int main(int argc, char **argv) 
{
	
	double *S; /*Overlap matrix*/
	double *H; /*Core Hamiltonian*/
	double *D_new; /*Density matrix current*/
	double *D_old; /*Density matrix previous*/
	double *F; /*Fock matrix*/
	erd_t *erd_inp;
	int maxit = MAX_IT; /*Max. no of iterations*/
	int conv; /*convergence flag*/
	basis_set_t *basis; /*Basis set Structure*/
	double *scratch;
	double *S_sinv;

	/*Initialize Basis Set*/
	basis = create_basis_set();
	load_basis_set(basis, argv[1]);
	preprocess_basis_set(basis);
	
	fprintf(stderr, "\n DEBUG: Initialized basis set successfully \n");
	fflush(stderr);
	
	/*Allocate memory for matrices*/
	F = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	H = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	D_old = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	D_new = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	S = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	S_sinv = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	scratch = (double *)malloc(basis->nfunctions*basis->nfunctions*sizeof(double));
	

	fprintf(stderr, "\n DEBUG: Initialized matrices successfully \n");
	fflush(stderr);

	/*TODO: Initialize integrals package*/
	erd_inp = init_erd(basis);
	
	fprintf(stderr, "\n DEBUG: Initialized ERD successfully \n");
	fflush(stderr);

	/*Compute Core Hamiltonian and overlap matrix*/
	compute_S(S, basis, 0, basis->nshells - 1, 0, basis->nshells - 1);
	compute_H(H, basis, 0, basis->nshells - 1, 0, basis->nshells - 1);

	fprintf(stderr, "\n DEBUG: Computed One electron stuff successfully \n");
	fflush(stderr);
	/*Compute square root inverse of S*/
	sqrtinv_matrix(basis->nfunctions, S, S_sinv, scratch);

	fprintf(stderr, "\n DEBUG: Computed Square root inverse of S \n");
	fflush(stderr);
	
	/*TODO: SCF iterate until converged*/
	
	conv = sscf(H, S_sinv, basis->nfunctions, basis->nelectrons, maxit, D_old, D_new, F);
	
	if (!conv) {
		fprintf(stderr, "\n DEBUG: Convergence not achieved in %d iterations \n", maxit);
		fflush(stderr);
	} else {
		fprintf(stderr, "\n DEBUG: Convergence achieved! \n");
		fflush(stderr);
	}
		
	/*TODO: Print energy, and eigenvalues*/

	/*TODO: clean exit*/
	destroy_basis_set(basis);
	destroy_erd(erd_inp);
	free(F);
	free(H);
	free(S);
	free(S_sinv);
	free(scratch);
	free(D_old);
	free(D_new);
	
	fprintf(stderr, "\n DEBUG: All allocated memory freed, Exiting. \n");
	fflush(stderr);

	return 0;

}

int sscf(double *H, double *S_sinv, int n, int n_ele, int maxit, double *D_old, double *D_new, double *F) 
{

	/*TODO: Write this*/
	
	printmatCM("Core Hamiltonian", H, n, n);
	
	compute_D(n, n_ele, H, D_new);
	printmatCM("Test matrix", D_new, n, n);
	print_evals(D_new, n);
	return 0;
}	

void compute_D(int n, int n_ele, double *F, double *D) 
{
	
	int lwork;
	double *work;
	int info;
	// eigen values
	double *w;
	//eigen vectors
	double *ev;
	int m;

	m = n_ele/2;

	w = (double *)malloc (n * sizeof(double));
	ev = (double *)malloc (n * n * sizeof(double));
	memcpy (ev, F, n * n * sizeof(double));

	lwork = (3 * n - 1 > 1 ? 3 * n - 1 : 1);
	work = (double *)malloc (lwork * sizeof(double));

	dsyev_ ("V", "L", &n, ev, &n, w, work, &lwork, &info);
	assert (info == 0);
	
	cblas_dgemm (CblasColMajor, CblasTrans, CblasNoTrans, n, n, m,
		1.0, ev, n, ev, n, 0.0, D, n);
	

	free (ev);
	free (w);
	free (work);

	return;
}
	
void print_evals(double *A, int n) 
{
	
	
	int lwork;
	double *work;
	int info;
	// eigen values
	double *w;
	//eigen vectors
	double *ev;
	 
	int i;

	w = (double *)malloc (n * sizeof(double));
	ev = (double *)malloc (n * n * sizeof(double));
	memcpy (ev, A, n * n * sizeof(double));

	lwork = (3 * n - 1 > 1 ? 3 * n - 1 : 1);
	work = (double *)malloc (lwork * sizeof(double));

	dsyev_ ("V", "L", &n, ev, &n, w, work, &lwork, &info);
	assert (info == 0);
	
	fprintf(stderr, "\n evals:");
	for(i = 0; i<n; i++) {
		fprintf(stderr, " %lf ,", w[i]);
	}
	fprintf(stderr, "\n");
	
	free (ev);
	free (w);
	free (work);
	
	return;

}

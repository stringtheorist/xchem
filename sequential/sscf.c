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
#include "basis_set.h"
#include "utils.h"
#include "erd_integral.h"
#include "one_electron.h"
#include "sscf.h"

#define TOLER 1.0e-5
#define MAX_IT 5

int main (int argc, char **argv) 
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
	basis = create_basis_set ();
	load_basis_set (basis, argv[1]);
	preprocess_basis_set (basis);
	
	fprintf (stderr, "\n DEBUG: Initialized basis set successfully \n");
	fflush (stderr);
	
	/*Allocate memory for matrices*/
	F = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	H = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	D_old = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	D_new = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	S = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	S_sinv = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	scratch = (double *)malloc (basis->nfunctions * basis->nfunctions * sizeof(double));
	
	fprintf (stderr, "\n DEBUG: Initialized matrices successfully \n");
	fflush (stderr);

	/*TODO: Initialize integrals package*/
	erd_inp = init_erd (basis);
	
	fprintf (stderr, "\n DEBUG: Initialized ERD successfully \n");
	fflush (stderr);

	/*Compute Core Hamiltonian and overlap matrix*/
	compute_S (S, basis, 0, basis->nshells - 1, 0, basis->nshells - 1);
	compute_H (H, basis, 0, basis->nshells - 1, 0, basis->nshells - 1);

	fprintf (stderr, "\n DEBUG: Computed One electron stuff successfully \n");
	fflush (stderr);
	
	printmatCM ("S", S, basis->nfunctions, basis->nfunctions);
	printmatCM ("H", H, basis->nfunctions, basis->nfunctions);

	/*Compute square root inverse of S*/
	sqrtinv_matrix (basis->nfunctions, S, S_sinv, scratch);

	fprintf (stderr, "\n DEBUG: Computed Square root inverse of S \n");
	fflush (stderr);
	
	/*TODO: SCF iterate until converged*/
	conv = sscf (basis, erd_inp, H, S_sinv, basis->nfunctions, basis->nelectrons, maxit, D_old, D_new, F);
	
	if (!conv) {
		fprintf (stderr, "\n DEBUG: Convergence not achieved in %d iterations \n", maxit);
		fflush (stderr);
	} else {
		fprintf (stderr, "\n DEBUG: Convergence achieved! \n");
		fflush (stderr);
	}
		
	/*TODO: Print energy, and eigenvalues*/

	/*TODO: clean exit*/
	destroy_basis_set (basis);
	destroy_erd (erd_inp);
	free (F);
	free (H);
	free (S);
	free (S_sinv);
	free (scratch);
	free (D_old);
	free (D_new);
	
	fprintf (stderr, "\n DEBUG: All allocated memory freed, Exiting. \n");
	fflush (stderr);

	return 0;

}

int sscf (basis_set_t *basis, erd_t *erd_inp, double *H, double *S_sinv, int n, int n_ele, int maxit, double *D_old, 
	double *D_new, double *F) 
{
	double *int_buffer;
	double *tmp;
	double *F_t;
	int conv = 0;
	int iter = 0;
	int max_funcs;
	int max_buffer_dim;

	max_funcs =  2 * basis->max_momentum + 1;
	max_buffer_dim = max_funcs * max_funcs * max_funcs * max_funcs;
	
	int_buffer = (double *)malloc (max_buffer_dim * sizeof(double));
	tmp = (double *)malloc (n * n * sizeof(double));
	F_t = (double *)malloc (n * n * sizeof(double));
	memset (int_buffer, 0, max_buffer_dim * sizeof(double));
	memset (tmp, 0, n * n * sizeof(double));
	memset (F_t, 0, n * n * sizeof(double));
	memset (D_old, 0, n * n * sizeof(double));
	memset (D_new, 0, n * n * sizeof(double));
	memcpy (F, H, n * n * sizeof(double));


	do {
		memcpy (F, H, n * n * sizeof(double));
		/*Build F*/
		build_fock (basis, erd_inp, int_buffer, D_new, F);
		
		memcpy (D_old, D_new, n * n * sizeof(double));
		/*Transform F*/
		cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
			1.0, S_sinv, n, F, n, 0.0, tmp, n);

		cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
			1.0, tmp, n, S_sinv, n, 0.0, F_t, n);
	
		/*Compute D*/
		compute_D (n, n_ele, F_t, D_new);

		/*Transform D*/
		cblas_dgemm (CblasColMajor, CblasNoTrans, CblasNoTrans, n, n, n,
			1.0, S_sinv, n, D_new, n, 0.0, tmp, n);

		cblas_dgemm (CblasColMajor, CblasNoTrans, CblasTrans, n, n, n,
			1.0, tmp, n, S_sinv, n, 0.0, D_new, n);


		/*Extrapolate*/
		iter++;
		/*Check energy convergence*/
	} while ((iter < maxit));

	/* printmatCM ("Core Hamiltonian", H, n, n); */
	
	
	/* printmatCM ("Test matrix", D_new, n, n); */
	/* print_evals (D_new, n); */
	
	fprintf (stderr, "\n Final Energy: %lf \n", calc_hf_ene (D_new, F, H, n));
	printmatCM ("Final D", D_new, n, n);
	printmatCM ("Final F", F, n, n);
	


	free (tmp);
	free (F_t);
	free (int_buffer);
	return 0;
}	

double calc_hf_ene (double *D, double *F, double *H, int n)
{
	double ene;

	ene = cblas_ddot (n * n, D, 1, F, 1) + cblas_ddot (n * n, D, 1, H, 1);
	
	return ene;
} 


void compute_D (int n, int n_ele, double *F, double *D) 
{
	
	int lwork;
	double *work;
	int info;
	// eigen values
	double *w;
	//eigen vectors
	double *ev;
	double *tmp;
	int m;

	m = n_ele/2;

	w = (double *)malloc (n * sizeof(double));
	ev = (double *)malloc (n * n * sizeof(double));
	tmp = (double *)malloc (n * n * sizeof(double));
	memcpy (ev, F, n * n * sizeof(double));

	lwork = (3 * n - 1 > 1 ? 3 * n - 1 : 1);
	work = (double *)malloc (lwork * sizeof(double));

	dsyev_ ("V", "L", &n, ev, &n, w, work, &lwork, &info);
	assert (info == 0);
	memcpy (tmp, ev, n * n * sizeof(double));
	
	cblas_dgemm (CblasColMajor, CblasNoTrans, CblasTrans, n, n, m,
		1.0, tmp, n, ev, n, 0.0, D, n);
	
	free (tmp);
	free (ev);
	free (w);
	free (work);

	return;
}
	
void print_evals (double *A, int n) 
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
	for(i = 0; i < n; i++) {
		fprintf (stderr, " %lf ,", w[i]);
	}
	fprintf (stderr, "\n");
	
	free (ev);
	free (w);
	free (work);
	
	return;

}


void update_fock_j_q (double *integrals, int A, int B, int C, int D, basis_set_t *basis, double *dens, double *F) 
{
	int il;
	int jl;
	int kl;
	int ll;
	int ih;
	int jh;
	int kh;
	int lh;
	int i, j, k, l;
	int ix, jx, kx, lx;
	int d1, d2, d3, d4;
	int n;
	
	n = basis->nfunctions;
	
	il = basis->f_start_id[A];
	jl = basis->f_start_id[B];
	kl = basis->f_start_id[C];
	ll = basis->f_start_id[D];
	
	ih = basis->f_end_id[A];
	jh = basis->f_end_id[B];
	kh = basis->f_end_id[C];
	lh = basis->f_end_id[D];

	d1 = ih - il + 1;
	d2 = jh - jl + 1;
	d3 = kh - kl + 1;
	d4 = lh - ll + 1;

	for (i = il; i <= ih; i++) {
		for (j = jl; j <= jh; j++) {
			for (k = kl; k <= kh; k++) {
				for (l = ll; l <= lh; l++) {
					ix = i - il;
					jx = j - jl;
					kx = k - kl;
					lx = l - ll;
					F[i + n * j] += 2.0 * dens[k + n * l] 
						* integrals[ix + d1 * jx + d1 * d2 * kx + d1 * d2 * d3 * lx];
				}
			}
		}
	}
	
	
	return;
} 

void update_fock_k_q (double *integrals, int A, int B, int C, int D, basis_set_t *basis, double *dens, double *F) 
{
	int il;
	int jl;
	int kl;
	int ll;
	int ih;
	int jh;
	int kh;
	int lh;
	int d1, d2, d3, d4;
	int ix, jx, kx, lx;
	int i, j, k, l;
	int n; 

	n = basis->nfunctions;
	
	il = basis->f_start_id[A];
	jl = basis->f_start_id[B];
	kl = basis->f_start_id[C];
	ll = basis->f_start_id[D];
	
	ih = basis->f_end_id[A];
	jh = basis->f_end_id[B];
	kh = basis->f_end_id[C];
	lh = basis->f_end_id[D];
	
	d1 = ih - il + 1;
	d2 = jh - jl + 1;
	d3 = kh - kl + 1;
	d4 = lh - ll + 1;

	
	for (i = il; i <= ih; i++) {
		for (j = jl; j <= jh; j++) {
			for (k = kl; k <= kh; k++) {
				for (l = ll; l <= lh; l++) {
					ix = i - il;
					jx = j - jl;
					kx = k - kl;
					lx = l - ll;
					F[j +  n * k] -= 1.0 * dens[i + n * l] 
						* integrals[ix + d1 * jx + d1 * d2 * kx + d1 * d2 * d3 * lx];
				}
			}
		}
	}


	return;
}

void build_fock (basis_set_t *basis, erd_t *erd_inp, double *int_buffer, double *D_new, double *F) 
{
	
	int A, B, C, D;

	for (A = 0; A < basis->nshells; A++) {
		for (B = 0; B < basis->nshells; B++) {
			for (C = 0; C < basis->nshells; C++) {
				for (D = 0; D < basis->nshells; D++) {
					compute_shell_quartet (int_buffer, A, B, C, D, basis, erd_inp);
					/*contract the shell quartet into F*/
					update_fock_j_q (int_buffer, A, B, C, D, basis, D_new, F);
					update_fock_k_q (int_buffer, A, B, C, D, basis, D_new, F);
					
				}
			}
		}
	}



	return;
}

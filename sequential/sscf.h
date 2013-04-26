/* Xchem: An interface to the ERD integrals package. */

/* Copyright (C) Aftab Patel, Xing Liu  */

/* See ../COPYRIGHT and ../LISENCE */


#ifndef __SSCF_H__
#define __SSCF_H__

#include "basis_set.h"
#include "erd_integral.h"

double compute_trace (double *A, int n);
int sscf (basis_set_t *basis, erd_t *erd_inp, double *H, double *S, 
	double *S_sinv, int n, int n_ele, int maxit, int diis_lim,
	double *D_old, 
	double *D_new, double *F); 
void compute_D (int n, int n_ele, double *F, double *D);
void print_evals (double *A, int n);
void update_fock_j_q (double *integrals, int A, int B, int C, int D, basis_set_t *basis, double *dens, double *F);
void update_fock_k_q (double *integrals, int A, int B, int C, int D, basis_set_t *basis, double *dens, double *F);
void build_fock (basis_set_t *basis, erd_t *erd_inp, double *int_buffer, double *D_new, double *F);
double calc_hf_ene (double *D, double *F, double *H, int n);

#endif

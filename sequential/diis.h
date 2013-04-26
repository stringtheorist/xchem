/* Xchem: An interface to the ERD integrals package. */

/* Copyright (C) Aftab Patel, Xing Liu  */

/* See ../COPYRIGHT and ../LISENCE */


#ifndef __DIIS_H__
#define __DIIS_H__

#include "vec_stack.h"

typedef struct diis_t_ {
	
	vec_stack *errs;
	vec_stack *vecs;
	int fails;
	int restarts;
	int diis_lim;
	
} diis_t;

void destroy_diis(diis_t *inp);

diis_t *init_diis(int diis_lim, int vec_len);

int diis_step (double **err, double **v, int steps, int dim, double *v_res);

void diis_extrap (diis_t *inp, double *result);

// void diis_extrap (double **err, double **v, int dim, int diis_lim, int iter, double *v_res);


#endif

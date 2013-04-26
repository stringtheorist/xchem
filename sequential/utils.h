/* Xchem: An interface to the ERD integrals package. */

/* Copyright (C) Aftab Patel, Xing Liu  */

/* See ../COPYRIGHT and ../LISENCE */


#ifndef __UTILS_H__
#define __UTILS_H__

#if 0
double **init_vecs_queue (int num_vecs, int vec_len);

void destroy_vecs_queue (double **vecs_queue, int num_vecs);

void clear_vecs_queue (int num_vecs, int vec_len, double **vecs_queue);
#endif


void sqrt_matrix (int n, double *A, double *A2, double *workm);

void sqrtinv_matrix (int n, double *A, double *A2, double *workm);

void trace_matrix (int n, double *A, double *tr);

void transpose (double *Ato, double *Afrom, int nrows, int ncols);

void printmatRM (char *name, double *A, int nrows, int ncols);

void printmatCM (char *name, double *A, int nrows, int ncols);

void initomp (int nthreads, int verbose);


#endif /* __UTILS_H__ */

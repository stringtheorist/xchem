#ifndef __UTILS_H__
#define __UTILS_H__


void sqrt_matrix (int n, double *A, double *A2, double *workm);

void sqrtinv_matrix (int n, double *A, double *A2, double *workm);

void trace_matrix (int n, double *A, double *tr);

void transpose (double *Ato, double *Afrom, int nrows, int ncols);

void printmatRM (char *name, double *A, int nrows, int ncols);

void printmatCM (char *name, double *A, int nrows, int ncols);

void initomp (int nthreads, int verbose);


#endif /* __UTILS_H__ */

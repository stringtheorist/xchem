#ifndef __SSCF_H__
#define __SSCF_H__


int sscf (double *H, double *S_sinv, int n, int n_ele, int maxit, double *D_old, double *D_new, double *F);
void compute_D (int n, int n_ele, double *F, double *D);
void print_evals (double *A, int n);


#endif

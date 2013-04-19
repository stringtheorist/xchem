#ifndef __ERD_INTEGRAL_H__
#define __ERD_INTEGRAL_H__


#include "basis_set.h"


#define ERD_SCREEN 0
#define ERD_SPHERIC 1


// Main fortran function of the ERD package.

extern void erd__gener_eri_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2, int *ncgto3,
                                   int *ncgto4, int *npgto1, int *npgto2,
                                   int *npgto3, int *npgto4, int *shell1,
                                   int *shell2, int *shell3, int *shell4,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   double *x3, double *y3, double *z3,
                                   double *x4, double *y4, double *z4,
                                   double *alpha, double *cc, int *ccbeg,
                                   int *ccend, int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst,
                                   double *zcore);

extern void erd__memory_eri_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2, int *ncgto3,
                                    int *ncgto4, int *npgto1, int *npgto2,
                                    int *npgto3, int *npgto4, int *shell1,
                                    int *shell2, int *shell3, int *shell4,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    double *x3, double *y3, double *z3,
                                    double *x4, double *y4, double *z4,
                                    double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt, int *zmin,
                                    int *zopt);


typedef struct _erd_t
{
    int shell1;
    int shell2;
    int shell3;
    int shell4;
    int blocksize;
    int ncoef;
    int nalpha;
    int ncsum;
    int ncgto1;
    int ncgto2;
    int ncgto3;
    int ncgto4;
    int npgto1;
    int npgto2;
    int npgto3;
    int npgto4;
    int imax;
    int zmax;
    int spheric;
    int screen;
    int cc_beg[4];
    int cc_end[4];
    double *cc;
    double *alpha;
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    double x3;
    double y3;
    double z3;
    double x4;
    double y4;
    double z4;

    int nbatch;
    int nfirst;
    int *icore;
    double *zcore;

    int fp_memory_min;
    int fp_memory_opt;
    int int_memory_min;
    int int_memory_opt;
} erd_t;


void compute_shell_quartet (double *integrals, int A, int B, int C, int D,
                            basis_set_t * basis, erd_t * erd);

erd_t *init_erd (basis_set_t * basis);

void destroy_erd (erd_t * erd);


#endif /* __ERD_INTEGRAL_H__ */

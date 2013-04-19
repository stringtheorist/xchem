#ifndef __OED_INTEGRAL_H__
#define __OED_INTEGRAL_H__


#include "basis_set.h"


#define OED_SPHERIC 1
#define OED_SCREEN 0


extern void oed__gener_kin_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2,
                                   int *npgto1, int *npgto2,
                                   int *shell1, int *shell2,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   double *alpha, double *cc,
                                   int *ccbeg, int *ccend, int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst, double *zcore);


extern void oed__gener_nai_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2,
                                   int *npgto1, int *npgto2,
                                   int *shell1, int *shell2,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   int *natoms,
                                   double *xn, double *yn, double *zn,
                                   double *ncharge, double *alpha, double *cc,
                                   int *ccbeg, int *ccend,
                                   int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst, double *zcore);


extern void oed__gener_ovl_batch_ (int *imax, int *zmax,
                                   int *nalpha, int *ncoeff, int *ncsum,
                                   int *ncgto1, int *ncgto2,
                                   int *npgto1, int *npgto2,
                                   int *shell1, int *shell2,
                                   double *x1, double *y1, double *z1,
                                   double *x2, double *y2, double *z2,
                                   double *alpha, double *cc, int *ccbeg, int *ccend,
                                   int *spheric, int *screen,
                                   int *icore, int *nbatch, int *nfirst, double *zcore);


extern void oed__memory_ovl_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2,
                                    int *npgto1, int *npgto2,
                                    int *shell1, int *shell2,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt,
                                    int *zmin, int *zopt);


extern void oed__memory_kin_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2,
                                    int *npgto1, int *npgto2,
                                    int *shell1, int *shell2,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt, int *zmin, int *zopt);

extern void oed__memory_nai_batch_ (int *nalpha, int *ncoeff,
                                    int *ncgto1, int *ncgto2,
                                    int *npgto1, int *npgto2,
                                    int *shell1, int *shell2,
                                    double *x1, double *y1, double *z1,
                                    double *x2, double *y2, double *z2,
                                    int *natoms, double *alpha, double *cc, int *spheric,
                                    int *imin, int *iopt, int *zmin, int *zopt);


#define KIN 0
#define OVL 1
#define POT 2
#define HHH 3


typedef struct _oed_t
{
    int type;
    int nalpha;
    int ncoeff;
    int ncgto1;
    int ncgto2;
    int npgto1;
    int npgto2;
    int shell1;
    int shell2;
    int natoms;
    int ncsum;
    int spheric;
    int screen;
    
    double x1;
    double y1;
    double z1;
    double x2;
    double y2;
    double z2;
    double *xn;
    double *yn;
    double *zn;
    double *ncharge;
    double *cc;
    double *alpha;
    int cc_beg[2];
    int cc_end[2];

    int imax;
    int zmax;
    double *zcore;
    int *icore;
    int nfirst;
    int nbatch;
    int blocksize;

    int fp_memory_min;
    int fp_memory_opt;
    int int_memory_min;
    int int_memory_opt;
} oed_t;


oed_t *init_oed (basis_set_t *basis, int type);

void destroy_oed (oed_t *oed);
    
void compute_shell_pair_kin (double *integrals, int A, int B,
                             oed_t *oed, basis_set_t *basis);

void compute_shell_pair_ovl (double *integrals, int A, int B,
                             oed_t *oed, basis_set_t *basis);

void compute_shell_pair_pot (double *integrals, int A, int B,
                             oed_t *oed, basis_set_t *basis);

void compute_shell_pair_H (double *integrals, int A, int B,
                           oed_t *oed, basis_set_t * basis);


#endif /* __OED_INTEGRAL_H__ */

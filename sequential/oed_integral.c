#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "oed_integral.h"


static void config_oed (oed_t * oed, int A, int B, basis_set_t * basis)
{
    int coef_offset_A;
    int coef_offset_B;
    int exp_offset_A;
    int exp_offset_B;
    int cc_offset_A;
    int cc_offset_B;
    int alpha_offset_A;
    int alpha_offset_B;
    int dim_1;
    int dim_2;

    alpha_offset_A = 0;
    cc_offset_A = 0;
    alpha_offset_B = basis->nexp[A];
    cc_offset_B = basis->ncoef[A];
    oed->ncoeff = alpha_offset_B + basis->ncoef[B];
    oed->nalpha = cc_offset_B + basis->nexp[B];
    
    coef_offset_A = basis->coef_offset[A];
    exp_offset_A = basis->exp_offset[A];
    coef_offset_B = basis->coef_offset[B];
    exp_offset_B = basis->exp_offset[B];

    memcpy (&(oed->alpha[alpha_offset_A]), &(basis->exp[exp_offset_A]), sizeof(double) * basis->nexp[A]);
    memcpy (&(oed->alpha[alpha_offset_B]), &(basis->exp[exp_offset_B]), sizeof(double) * basis->nexp[B]);
    memcpy (&(oed->cc[cc_offset_A]), &(basis->coef[coef_offset_A]), sizeof(double) * basis->ncoef[A]);
    memcpy (&(oed->cc[cc_offset_B]), &(basis->coef[coef_offset_B]), sizeof(double) * basis->ncoef[B]);

    oed->npgto1 = basis->nprimitives[A];
    oed->npgto2 = basis->nprimitives[B];
    oed->cc_end[0] = oed->npgto1;
    oed->cc_end[1] = oed->npgto2;    
    oed->x1 = basis->x[A];
    oed->y1 = basis->y[A];
    oed->z1 = basis->z[A];
    oed->x2 = basis->x[B];
    oed->y2 = basis->y[B];
    oed->z2 = basis->z[B];
    oed->shell1 = basis->momentum[A];
    oed->shell2 = basis->momentum[B];

    dim_1 = 2 * oed->shell1 + 1;
    dim_2 = 2 * oed->shell2 + 1;
    oed->blocksize = dim_1 * dim_2;
}


static void oed_max_scratch (basis_set_t *basis, oed_t *oed)
{
    int max_momentum;
   
    config_oed (oed,
                basis->max_prim_id, basis->max_coef_id,
                basis); 
  
    max_momentum = basis->momentum[basis->max_momentum_id];
    oed->shell1 = max_momentum;
    oed->shell2 = max_momentum;
    oed->x1 = 1.0;
    oed->x2 = 2.0;
    oed->y1 = 1.0;
    oed->y2 = 2.0;
    oed->z1 = 1.0;
    oed->z2 = 2.0;

    if (oed->type == KIN)
    {
        oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                                &(oed->ncgto1), &(oed->ncgto2),
                                &(oed->npgto1), &(oed->npgto2),
                                &(oed->shell1), &(oed->shell2),
                                &(oed->x1), &(oed->y1), &(oed->z1),
                                &(oed->x2), &(oed->y2), &(oed->z2),
                                oed->alpha, oed->cc, &(oed->spheric),
                                &(oed->int_memory_min), &(oed->int_memory_opt),
                                &(oed->fp_memory_min), &(oed->fp_memory_opt));
    }
    else if (oed->type == OVL)
    {
        oed__memory_ovl_batch_ (&(oed->nalpha), &(oed->ncoeff),
                                &(oed->ncgto1), &(oed->ncgto2),
                                &(oed->npgto1), &(oed->npgto2),
                                &(oed->shell1), &(oed->shell2),
                                &(oed->x1), &(oed->y1), &(oed->z1),
                                &(oed->x2), &(oed->y2), &(oed->z2),
                                oed->alpha, oed->cc, &(oed->spheric),
                                &(oed->int_memory_min), &(oed->int_memory_opt),
                                &(oed->fp_memory_min), &(oed->fp_memory_opt));    
    }
    else if (oed->type == POT)
    {
        oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                                &(oed->ncgto1), &(oed->ncgto2),
                                &(oed->npgto1), &(oed->npgto2),
                                &(oed->shell1), &(oed->shell2),
                                &(oed->x1), &(oed->y1), &(oed->z1),
                                &(oed->x2), &(oed->y2), &(oed->z2),
                                &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),
                                &(oed->int_memory_min), &(oed->int_memory_opt),
                                &(oed->fp_memory_min), &(oed->fp_memory_opt));    
    }
    else if (oed->type == HHH)
    {
        int int_memory_min;
        int int_memory_opt;
        int fp_memory_min;
        int fp_memory_opt;
        oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                                &(oed->ncgto1), &(oed->ncgto2),
                                &(oed->npgto1), &(oed->npgto2),
                                &(oed->shell1), &(oed->shell2),
                                &(oed->x1), &(oed->y1), &(oed->z1),
                                &(oed->x2), &(oed->y2), &(oed->z2),
                                oed->alpha, oed->cc, &(oed->spheric),
                                &int_memory_min, &int_memory_opt,
                                &fp_memory_min, &fp_memory_opt);
                
        oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                                &(oed->ncgto1), &(oed->ncgto2),
                                &(oed->npgto1), &(oed->npgto2),
                                &(oed->shell1), &(oed->shell2),
                                &(oed->x1), &(oed->y1), &(oed->z1),
                                &(oed->x2), &(oed->y2), &(oed->z2),
                                &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),
                                &(oed->int_memory_min), &(oed->int_memory_opt),
                                &(oed->fp_memory_min), &(oed->fp_memory_opt));

        oed->int_memory_min = oed->int_memory_min > int_memory_min ?
            oed->int_memory_min : int_memory_min;
        oed->int_memory_opt = oed->int_memory_opt > int_memory_opt ?
            oed->int_memory_opt : int_memory_opt;
        oed->fp_memory_min = oed->fp_memory_min > fp_memory_min ?
            oed->fp_memory_min : fp_memory_min;
        oed->fp_memory_opt = oed->fp_memory_opt > fp_memory_opt ?
            oed->fp_memory_opt : fp_memory_opt;
    }
}


oed_t *init_oed (basis_set_t * basis, int type)
{
    oed_t *oed;

    oed = (oed_t *)malloc (sizeof(oed_t));
    assert (oed != NULL);   
    oed->type = type;
    oed->cc = (double *)malloc (2 * basis->max_coef * sizeof(double));
    assert (oed->cc != NULL);
    oed->alpha = (double *)malloc (2 * basis->max_exp * sizeof(double));
    assert (oed->alpha != NULL);

    oed->ncsum = 2;
    oed->ncgto1 = 1;
    oed->ncgto2 = 1;
    oed->cc_beg[0] = 1;
    oed->cc_beg[1] = 1;
    oed->natoms = basis->natoms;
    oed->xn = basis->xn;
    oed->yn = basis->yn;
    oed->zn = basis->zn;
    oed->ncharge = basis->ncharge;
    oed->spheric = OED_SPHERIC;
    oed->screen = OED_SCREEN;

    oed_max_scratch (basis, oed);
    oed->zcore = (double *)malloc (oed->fp_memory_opt * sizeof(double));
    oed->icore = (int *)malloc (oed->int_memory_opt * sizeof(int));
    assert (oed->zcore != NULL);
    assert (oed->icore != NULL);
    oed->zmax = oed->fp_memory_opt;
    oed->imax = oed->int_memory_opt;
    
    return oed;
}


void destroy_oed (oed_t * oed)
{
    free (oed->zcore);
    free (oed->icore);
    free (oed->alpha);
    free (oed->cc);
    free (oed);
}


void compute_shell_pair_kin (double *integrals, int A, int B,
                             oed_t *oed, basis_set_t * basis)
{
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    
    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_kin_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           oed->alpha, oed->cc,
                           oed->cc_beg, oed->cc_end, &(oed->spheric), &(oed->screen),
                           oed->icore, &(oed->nbatch), &(oed->nfirst), oed->zcore);

    if (oed->nbatch != 0)
    {
        memcpy (integrals, &oed->zcore[oed->nfirst - 1], sizeof(double) * oed->blocksize);
    }
    else
    {
        memset (integrals, 0, sizeof(double) * oed->blocksize);    
    }
}


void compute_shell_pair_ovl (double *integrals, int A, int B,
                             oed_t *oed, basis_set_t * basis)
{
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_ovl_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    
    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_ovl_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           oed->alpha, oed->cc, oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, &(oed->nbatch), &(oed->nfirst), oed->zcore);

    if (oed->nbatch != 0)
    {
        memcpy (integrals, &oed->zcore[oed->nfirst - 1], sizeof(double) * oed->blocksize);
    }
    else
    {
        memset (integrals, 0, sizeof(double) * oed->blocksize);    
    }
}


void compute_shell_pair_pot (double *integrals, int A, int B,
                             oed_t *oed, basis_set_t * basis)
{
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    
    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_nai_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           &(oed->natoms),
                           oed->xn, oed->yn, oed->zn,
                           oed->ncharge, oed->alpha, oed->cc,
                           oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, &(oed->nbatch), &(oed->nfirst), oed->zcore);

    if (oed->nbatch != 0)
    {
        memcpy (integrals, &oed->zcore[oed->nfirst - 1], sizeof(double) * oed->blocksize);
    }
    else
    {
        memset (integrals, 0, sizeof(double) * oed->blocksize);    
    }
}


void compute_shell_pair_H (double *integrals, int A, int B,
                           oed_t *oed, basis_set_t * basis)
{
    int i;
    config_oed (oed, A, B, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    oed__memory_kin_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);

    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);
    oed__memory_nai_batch_ (&(oed->nalpha), &(oed->ncoeff),
                            &(oed->ncgto1), &(oed->ncgto2),
                            &(oed->npgto1), &(oed->npgto2),
                            &(oed->shell1), &(oed->shell2),
                            &(oed->x1), &(oed->y1), &(oed->z1),
                            &(oed->x2), &(oed->y2), &(oed->z2),
                            &(oed->natoms), oed->alpha, oed->cc, &(oed->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);

    assert (fp_memory_opt <= oed->fp_memory_opt);
    assert (int_memory_opt <= oed->int_memory_opt);   
#endif

    oed__gener_kin_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           oed->alpha, oed->cc, oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, &(oed->nbatch), &(oed->nfirst), oed->zcore);
    if (oed->nbatch != 0)
    {
        memcpy (integrals, &oed->zcore[oed->nfirst - 1], sizeof(double) * oed->blocksize);
    }
    else
    {
        memset (integrals, 0, sizeof(double) * oed->blocksize);    
    }
    
 //   config_oed (oed, A, B, basis);
    oed__gener_nai_batch_ (&(oed->imax), &(oed->zmax),
                           &(oed->nalpha), &(oed->ncoeff), &(oed->ncsum),
                           &(oed->ncgto1), &(oed->ncgto2),
                           &(oed->npgto1), &(oed->npgto2),
                           &(oed->shell1), &(oed->shell2),
                           &(oed->x1), &(oed->y1), &(oed->z1),
                           &(oed->x2), &(oed->y2), &(oed->z2),
                           &(oed->natoms),
                           oed->xn, oed->yn, oed->zn,
                           oed->ncharge, oed->alpha, oed->cc,
                           oed->cc_beg, oed->cc_end,
                           &(oed->spheric), &(oed->screen),
                           oed->icore, &(oed->nbatch), &(oed->nfirst), oed->zcore);
    if (oed->nbatch != 0)
    {
        for (i = 0; i < oed->blocksize; i++)
        {
            integrals[i] += oed->zcore[oed->nfirst - 1 + i];
        }
    }
}

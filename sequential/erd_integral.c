/* Xchem: An interface to the ERD integrals package. */

/* Copyright (C) Aftab Patel, Xing Liu  */

/* See ../COPYRIGHT and ../LISENCE */



#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>

#include "common.h"
#include "erd_integral.h"


static void config_erd (erd_t * erd, int A, int B, int C, int D, basis_set_t * basis)
{
    int coef_offset_A;
    int coef_offset_B;
    int coef_offset_C;
    int coef_offset_D;
    int exp_offset_A;
    int exp_offset_B;
    int exp_offset_C;
    int exp_offset_D;
    int cc_offset_A;
    int cc_offset_B;
    int cc_offset_C;
    int cc_offset_D;
    int alpha_offset_A;
    int alpha_offset_B;
    int alpha_offset_C;
    int alpha_offset_D;
    int dim_1;
    int dim_2;
    int dim_3;
    int dim_4;
        
    alpha_offset_A = 0;
    cc_offset_A = 0;
    alpha_offset_B = basis->nexp[A];
    cc_offset_B = basis->ncoef[A];
    alpha_offset_C = alpha_offset_B + basis->nexp[B];
    cc_offset_C = cc_offset_B + basis->ncoef[B];
    alpha_offset_D = alpha_offset_C + basis->nexp[C];
    cc_offset_D = cc_offset_C + basis->ncoef[C];
    erd->nalpha = alpha_offset_D + basis->nexp[D]; 
    erd->ncoef =  cc_offset_D + basis->ncoef[D];
    
    coef_offset_A = basis->coef_offset[A];
    exp_offset_A = basis->exp_offset[A];
    coef_offset_B = basis->coef_offset[B];
    exp_offset_B = basis->exp_offset[B];
    coef_offset_C = basis->coef_offset[C];
    exp_offset_C = basis->exp_offset[C];
    coef_offset_D = basis->coef_offset[D];
    exp_offset_D = basis->exp_offset[D];
    
    memcpy (&(erd->alpha[alpha_offset_A]), &(basis->exp[exp_offset_A]), sizeof(double) * basis->nexp[A]);
    memcpy (&(erd->alpha[alpha_offset_B]), &(basis->exp[exp_offset_B]), sizeof(double) * basis->nexp[B]);
    memcpy (&(erd->alpha[alpha_offset_C]), &(basis->exp[exp_offset_C]), sizeof(double) * basis->nexp[C]);
    memcpy (&(erd->alpha[alpha_offset_D]), &(basis->exp[exp_offset_D]), sizeof(double) * basis->nexp[D]); 
    memcpy (&(erd->cc[cc_offset_A]), &(basis->coef[coef_offset_A]), sizeof(double) * basis->ncoef[A]);
    memcpy (&(erd->cc[cc_offset_B]), &(basis->coef[coef_offset_B]), sizeof(double) * basis->ncoef[B]);
    memcpy (&(erd->cc[cc_offset_C]), &(basis->coef[coef_offset_C]), sizeof(double) * basis->ncoef[C]);
    memcpy (&(erd->cc[cc_offset_D]), &(basis->coef[coef_offset_D]), sizeof(double) * basis->ncoef[D]);

    erd->npgto1 = basis->nprimitives[A];
    erd->npgto2 = basis->nprimitives[B];
    erd->npgto3 = basis->nprimitives[C];
    erd->npgto4 = basis->nprimitives[D];
    erd->cc_end[0] = erd->npgto1;
    erd->cc_end[1] = erd->npgto2;
    erd->cc_end[2] = erd->npgto3;
    erd->cc_end[3] = erd->npgto4;
    erd->x1 = basis->x[A];
    erd->y1 = basis->y[A];
    erd->z1 = basis->z[A];
    erd->x2 = basis->x[B];
    erd->y2 = basis->y[B];
    erd->z2 = basis->z[B];
    erd->x3 = basis->x[C];
    erd->y3 = basis->y[C];
    erd->z3 = basis->z[C];
    erd->x4 = basis->x[D];
    erd->y4 = basis->y[D];
    erd->z4 = basis->z[D];
    erd->shell1 = basis->momentum[A];
    erd->shell2 = basis->momentum[B];
    erd->shell3 = basis->momentum[C];
    erd->shell4 = basis->momentum[D]; 
    dim_1 = 2 * erd->shell1 + 1;
    dim_2 = 2 * erd->shell2 + 1;
    dim_3 = 2 * erd->shell3 + 1;
    dim_4 = 2 * erd->shell4 + 1;
    erd->blocksize = dim_1 * dim_2 * dim_3 * dim_4;
}


static void erd_max_scratch (basis_set_t *basis, erd_t *erd)
{
    int max_momentum;
    
    config_erd (erd,
                basis->max_prim_id, basis->max_prim_id,
                basis->max_prim_id, basis->max_prim_id,
                basis); 
  
    max_momentum = basis->momentum[basis->max_momentum_id];
    erd->shell1 = max_momentum;
    erd->shell2 = max_momentum;
    erd->shell3 = max_momentum;
    erd->shell4 = max_momentum;
    erd->x1 = 1.0;
    erd->x2 = 2.0;
    erd->x3 = 3.0;
    erd->x4 = 4.0;
    erd->y1 = 1.0;
    erd->y2 = 2.0;
    erd->y3 = 3.0;
    erd->y4 = 4.0;
    erd->z1 = 1.0;
    erd->z2 = 2.0;
    erd->z3 = 3.0;
    erd->z4 = 4.0;
  
    erd__memory_eri_batch_ (&(erd->nalpha), &(erd->ncoef),
                            &(erd->ncgto1), &(erd->ncgto2),
                            &(erd->ncgto3), &(erd->ncgto4),
                            &(erd->npgto1), &(erd->npgto2),
                            &(erd->npgto3), &(erd->npgto4),
                            &(erd->shell1), &(erd->shell2),
                            &(erd->shell3), &(erd->shell4),
                            &(erd->x1), &(erd->y1), &(erd->z1),
                            &(erd->x2), &(erd->y2), &(erd->z2),
                            &(erd->x3), &(erd->y3), &(erd->z3),
                            &(erd->x4), &(erd->y4), &(erd->z4),
                            erd->alpha, erd->cc, &(erd->spheric),
                            &(erd->int_memory_min), &(erd->int_memory_opt),
                            &(erd->fp_memory_min), &(erd->fp_memory_opt));
}


erd_t *init_erd (basis_set_t * basis)
{
    erd_t *erd;

    erd = (erd_t *)malloc (sizeof(erd_t));
    assert (erd != NULL);

    erd->ncgto1 = 1;
    erd->ncgto2 = 1;
    erd->ncgto3 = 1;
    erd->ncgto4 = 1;    
    erd->cc_beg[0] = 1;
    erd->cc_beg[1] = 1;
    erd->cc_beg[2] = 1;
    erd->cc_beg[3] = 1;
    erd->ncsum = 4;
    erd->spheric = ERD_SPHERIC;
    erd->screen = ERD_SCREEN;

    erd->cc = (double *)malloc (4 * basis->max_coef * sizeof(double));
    erd->alpha = (double *)malloc (4 * basis->max_exp * sizeof(double));
    assert (erd->cc != NULL);
    assert (erd->alpha != NULL);

    erd_max_scratch (basis, erd);
    erd->zcore = (double *)malloc (erd->fp_memory_opt * sizeof(double));
    erd->icore = (int *)malloc (erd->int_memory_opt * sizeof(int));   
    assert (erd->zcore != NULL);
    assert (erd->icore != NULL);
    erd->zmax = erd->fp_memory_opt;
    erd->imax = erd->int_memory_opt;
    
    return erd;
}


void destroy_erd (erd_t * erd)
{
    free (erd->zcore);
    free (erd->icore);
    free (erd->alpha);
    free (erd->cc);
    free (erd);
}


void compute_shell_quartet (double *integrals, int A, int B, int C, int D,
                            basis_set_t * basis, erd_t * erd)
{
    config_erd (erd, A, B, C, D, basis);

#if ( _DEBUG_LEVEL_ == 3 )
    int int_memory_min;
    int int_memory_opt;
    int fp_memory_min;
    int fp_memory_opt;
    erd__memory_eri_batch_ (&(erd->nalpha), &(erd->ncoef),
                            &(erd->ncgto1), &(erd->ncgto2),
                            &(erd->ncgto3), &(erd->ncgto4),
                            &(erd->npgto1), &(erd->npgto2),
                            &(erd->npgto3), &(erd->npgto4),
                            &(erd->shell1), &(erd->shell2),
                            &(erd->shell3), &(erd->shell4),
                            &(erd->x1), &(erd->y1), &(erd->z1),
                            &(erd->x2), &(erd->y2), &(erd->z2),
                            &(erd->x3), &(erd->y3), &(erd->z3),
                            &(erd->x4), &(erd->y4), &(erd->z4),
                            erd->alpha, erd->cc, &(erd->spheric),
                            &int_memory_min, &int_memory_opt,
                            &fp_memory_min, &fp_memory_opt);
    assert (fp_memory_opt <= erd->fp_memory_opt);
    assert (int_memory_opt <= erd->int_memory_opt);   
#endif

    erd__gener_eri_batch_ (&(erd->imax), &(erd->zmax),
                           &(erd->nalpha), &(erd->ncoef),
                           &(erd->ncsum), &(erd->ncgto1),
                           &(erd->ncgto2), &(erd->ncgto3),
                           &(erd->ncgto4), &(erd->npgto1),
                           &(erd->npgto2), &(erd->npgto3),
                           &(erd->npgto4), &(erd->shell1),
                           &(erd->shell2), &(erd->shell3),
                           &(erd->shell4),
                           &(erd->x1), &(erd->y1), &(erd->z1),
                           &(erd->x2), &(erd->y2), &(erd->z2),
                           &(erd->x3), &(erd->y3), &(erd->z3),
                           &(erd->x4), &(erd->y4), &(erd->z4),
                           erd->alpha, erd->cc,
                           erd->cc_beg, erd->cc_end,
                           &(erd->spheric), &(erd->screen),
                           erd->icore, &(erd->nbatch),
                           &(erd->nfirst), erd->zcore);
    
    if (erd->nbatch != 0)
    {
        memcpy (integrals, &erd->zcore[erd->nfirst - 1], sizeof(double) * erd->blocksize);
    }
    else
    {
        memset (integrals, 0, sizeof(double) * erd->blocksize);    
    }
}

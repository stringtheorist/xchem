#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <sys/time.h>

#include "basis_set.h"
#include "common.h"


static int sum_int_array (int *array, int length)
{
    int sum;
    int i;
    
    sum = 0;
    for (i = 0; i < length; i++)
    {
        sum += array[i];
    }
    
    return sum;
}


static int read_double_array (double *in, FILE * fp, int in_length)
{
    int length;

    fread (&length, sizeof(int), 1, fp);
    if (length != in_length)
    {
        return -1;
    }
    fread (in, sizeof(double), length, fp);

    return length;
}


static int read_int_array (int *in, FILE * fp, int in_length)
{
    int length;

    fread (&length, sizeof(int), 1, fp);
    if (length != in_length)
    {
        return -1;
    }
    fread (in, sizeof(int), length, fp);

    return length;
}


static int probe_read (FILE * fp)
{
    int length;
    
    fread (&length, sizeof(int), 1, fp);
    
    return length;
}


static void normalize_coefficients (basis_set_t * basis)
{
    double sum;
    double temp;
    double temp2;
    double temp3;
    double xnorm;
    double a1;
    double a2;
    int o_exp;
    int o_coef;

    int i;
    int j;
    int k;
    int l;

    for (i = 0; i < basis->nshells; i++)
    {
        sum = 0.0;
        for (j = 0; j < basis->ncoef[i]; j++)
        {
            for (k = 0; k <= j; k++)
            {
                o_exp = sum_int_array (basis->nexp, i);
                o_coef = sum_int_array (basis->ncoef, i);
                a1 = basis->exp[o_exp + j];
                a2 = basis->exp[o_exp + k];
                temp = (basis->coef[o_coef + j]) * (basis->coef[o_coef + k]);
                temp2 = ((double) basis->momentum[i] + 1.5);
                temp3 = (2.0 * sqrt (a1 * a2) / (a1 + a2));
                temp3 = pow (temp3, temp2);
                temp = temp * temp3;
                sum = sum + temp;
                if (j != k)
                {
                    sum = sum + temp;
                }
            }
        }
        xnorm = 1.0 / sqrt (sum);
        for (l = 0; l < basis->ncoef[i]; l++)
        {
            o_coef = sum_int_array (basis->ncoef, i);
            basis->coef[o_coef + l] *= xnorm;
        }
    }
}


basis_set_t *create_basis_set (void)
{
    basis_set_t *basis;
    basis = (basis_set_t *)malloc (sizeof(basis_set_t));
    assert (basis != NULL);

    return basis;
}


// load basis set, only called by proc 0
void load_basis_set (basis_set_t *basis, char *dir)
{
    FILE *fp;
    char filename[1024]; 
    double basissize; 
    int i;
    struct timeval tv1;
    struct timeval tv2;
    double timepass;

    DPRINTF (1, "Loading basis set ...\n");
    gettimeofday (&tv1, NULL);
    
    sprintf (filename, "%s/alpha.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    basis->nshells = probe_read (fp);
    fclose (fp);

    sprintf (filename, "%s/exps.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    basis->nexp_tot = probe_read (fp);
    fclose (fp);

    sprintf (filename, "%s/cc.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    basis->ncoef_tot = probe_read (fp);
    fclose (fp);

    sprintf (filename, "%s/nchg.xcb", dir);
    fp = fopen (filename, "rb");
    basis->natoms = probe_read (fp);
    fclose (fp);

    /* initialize metametadata and Allocate memory */
    basis->ncharge = (double *)malloc (basis->natoms * sizeof(double));
    assert (basis->ncharge != NULL);
    basis->xn = (double *)malloc (basis->natoms * sizeof(double));
    assert (basis->xn != NULL);
    basis->yn = (double *)malloc (basis->natoms * sizeof(double));
    assert (basis->yn != NULL);
    basis->zn = (double *)malloc (basis->natoms * sizeof(double));
    assert (basis->zn != NULL);
    basis->ncoef = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->ncoef != NULL);
    basis->nexp = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->nexp != NULL);
    basis->momentum = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->momentum != NULL);
    basis->nprimitives = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->nprimitives != NULL);
    basis->f_start_id = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->f_start_id != NULL);
    basis->f_end_id = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->f_end_id != NULL);
    basis->x = (double *)malloc (basis->nshells * sizeof(double));
    assert (basis->x != NULL);
    basis->y = (double *)malloc (basis->nshells * sizeof(double));
    assert (basis->y != NULL);
    basis->z = (double *)malloc (basis->nshells * sizeof(double));
    assert (basis->z != NULL);
    basis->exp = (double *)malloc (basis->nexp_tot * sizeof(double));
    assert (basis->exp != NULL);    
    basis->coef = (double *)malloc (basis->ncoef_tot * sizeof(double));
    assert (basis->coef != NULL);

    basissize = basis->nshells * (3.0 * sizeof(double) + 8.0 * sizeof(int)) +
                (4.0 * basis->natoms + basis->nexp_tot + basis->ncoef_tot) * sizeof(double);
    basissize = basissize / 1024.0 / 1024.0;
    DPRINTF (1, "  using %.3lf MB\n", basissize);

    sprintf (filename, "%s/cc.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_double_array (basis->coef, fp, basis->ncoef_tot);
    fclose (fp);

    sprintf (filename, "%s/exps.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_double_array (basis->exp, fp, basis->nexp_tot);
    fclose (fp);

    sprintf (filename, "%s/centx.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_double_array (basis->x, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/centy.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_double_array (basis->y, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/centz.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_double_array (basis->z, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/alpha.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_int_array (basis->nexp, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/cc.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_int_array (basis->ncoef, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/am.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_int_array (basis->momentum, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/si.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_int_array (basis->f_start_id, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/ei.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_int_array (basis->f_end_id, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/np.meta.xcb", dir);
    fp = fopen (filename, "rb");
    assert (fp != NULL);
    read_int_array (basis->nprimitives, fp, basis->nshells);
    fclose (fp);

    sprintf (filename, "%s/nchg.xcb", dir);
    fp = fopen (filename, "rb");
    read_double_array (basis->ncharge, fp, basis->natoms);
    fclose (fp);

    sprintf (filename, "%s/nucx.xcb", dir);
    fp = fopen (filename, "rb");
    read_double_array (basis->xn, fp, basis->natoms);
    fclose (fp);

    sprintf (filename, "%s/nucy.xcb", dir);
    fp = fopen (filename, "rb");
    read_double_array (basis->yn, fp, basis->natoms);
    fclose (fp);

    sprintf (filename, "%s/nucz.xcb", dir);
    fp = fopen (filename, "rb");
    read_double_array (basis->zn, fp, basis->natoms);
    fclose (fp);
    basis->nfunctions = 0;
    for (i = 0; i < basis->nshells; i++)
    {
        basis->nfunctions += 2 * (basis->momentum[i]) + 1;
    }
    normalize_coefficients (basis);

    gettimeofday (&tv2, NULL);
    timepass = (tv2.tv_sec - tv1.tv_sec) +
        (tv2.tv_usec - tv1.tv_usec) / 1000.0 / 1000.0;
    DPRINTF (1, "  takes %.3lf secs\n", timepass);
    
    DPRINTF (1, "  Done\n");
}


void preprocess_basis_set (basis_set_t *basis)
{
    int countexp;
    int countcoef;   
    double sumcharge;
    int i;

    basis->coef_offset = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->coef_offset != NULL);
    basis->exp_offset = (int *)malloc (basis->nshells * sizeof(int));
    assert (basis->exp_offset != NULL);
    
    basis->max_momentum = 0;
    basis->max_momentum_id = 0;
    basis->max_coef_id = 0;
    basis->max_exp_id = 0;
    basis->max_coef = 0;
    basis->max_exp = 0;
    basis->max_prim_id = 0;
    basis->max_nprim = 0;
    countcoef = 0;
    countexp = 0;

    for (i = 0; i < basis->nshells; i++)
    {
        basis->coef_offset[i] = countcoef; 
        basis->exp_offset[i] = countexp;
        if (basis->momentum[i] > basis->max_momentum)
        {
            basis->max_momentum_id = i;
            basis->max_momentum = basis->momentum[i];
        }
        if (basis->ncoef[i] > basis->max_coef)
        {
            basis->max_coef_id = i;
            basis->max_coef = basis->ncoef[i];
        }
        if (basis->nexp[i] > basis->max_exp)
        {
            basis->max_exp_id = i;
            basis->max_exp = basis->nexp[i];
        }
        if (basis->nprimitives[i] > basis->max_nprim)
        {
            basis->max_prim_id = i;
            basis->max_nprim = basis->nprimitives[i];
        }
        countcoef = basis->ncoef[i] + countcoef;
        countexp = basis->nexp[i] + countexp;       
    }
    sumcharge = 0.0;
    for (i = 0; i < basis->natoms; i++)
    {
        sumcharge += basis->ncharge[i];
    }
    basis->nelectrons = (int)sumcharge;
    
#if 0
    printf ("%d:%d, %d:%d, %d:%d, %d:%d\n",
            basis->max_momentum_id, basis->max_momentum,
            basis->max_coef_id, basis->max_coef,
            basis->max_exp_id, basis->max_exp,
            basis->max_prim_id, basis->max_nprim);
#endif
}


void destroy_basis_set (basis_set_t * basis)
{
    free (basis->momentum);
    free (basis->ncoef);
    free (basis->nexp);
    free (basis->nprimitives);
    free (basis->f_start_id);
    free (basis->f_end_id);
    free (basis->exp);
    free (basis->coef);
    free (basis->x);
    free (basis->y);
    free (basis->z);
    free (basis->xn);
    free (basis->yn);
    free (basis->zn);
    free (basis->ncharge);
    free (basis->coef_offset);
    free (basis->exp_offset);

    free (basis);

    return;
}

#ifndef __BASIS_SET_H__
#define __BASIS_SET_H__


typedef struct _basis_set_t
{
    int nshells;  
    int natoms;
    int ncoef_tot;
    int nexp_tot;
    int nfunctions;
    
    int *ncoef;
    int *nexp;
    double *exp;
    double *coef;  
    int *momentum;
    int *nprimitives;
    int *f_start_id;
    int *f_end_id;    
    double *x;
    double *y;
    double *z;
    double *xn;
    double *yn;
    double *zn;
    double *ncharge;

    int *coef_offset;
    int *exp_offset;
    int nelectrons;
    int max_coef;
    int max_exp;
    int max_coef_id;
    int max_exp_id;
    int max_momentum;
    int max_momentum_id;
    int max_nprim;
    int max_prim_id;
} basis_set_t;


basis_set_t *create_basis_set (void);

void load_basis_set (basis_set_t *basis, char *dir);

void destroy_basis_set (basis_set_t * basis);

void preprocess_basis_set (basis_set_t *basis);




#endif /* __BASIS_SET_H__ */

#ifndef __DIIS_H__
#define __DIIS_H__

int diis_step(double **err, double **v, int steps, int dim, double *v_res);
#endif

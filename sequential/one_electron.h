#ifndef __ONE_ELECTRON_H__
#define __ONE_ELECTRON_H__


#include "basis_set.h"


void compute_S (double *S, basis_set_t * basis,
                int startshellrow, int endshellrow,
                int startshellcol, int endshellcol);

void compute_H (double *H, basis_set_t * basis,
                int startshellrow, int endshellrow,
                int startshellcol, int endshellcol);


#endif /* __ONE_ELECTRON_H__ */

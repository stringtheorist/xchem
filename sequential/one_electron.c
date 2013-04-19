#include <stdio.h>
#include <stdlib.h>
#include <assert.h>

#include "oed_integral.h"
#include "one_electron.h"


static void matrix_block_write (double *matrix, int startrow, int startcol, int ldm,
                                double *block, int nrows, int ncols)
{
	int i;
	int j;
	int k;
	int l;

	for (k = 0; k < nrows; k++)
	{
		for (l = 0; l < ncols; l++)
		{
			i = startrow + k;
			j = startcol + l;
			matrix[i * ldm + j] = block[k + nrows * l];
		}
	}
}


void compute_S (double *S, basis_set_t * basis,
                int startshellrow, int endshellrow,
                int startshellcol, int endshellcol)
{
	oed_t *oed;
	int A;
	int B;
	int dim_1;
	double *integrals;
	int row_id_1;
	int row_id_2;
	int col_id_1;
	int col_id_2;
	int start_row_id;
	int start_col_id;
	int end_col_id;
	int ldS;
	int nrows;
	int ncols;
	int startrow;
	int startcol;
    
	oed = init_oed (basis, OVL);
    
	dim_1 = 2 * basis->max_momentum + 1;
	integrals = (double *)malloc (dim_1 * dim_1 * sizeof(double));
	assert (integrals != NULL);

	start_row_id = basis->f_start_id[startshellrow];
	start_col_id = basis->f_start_id[startshellcol];
	end_col_id = basis->f_end_id[endshellcol];
	ldS = end_col_id - start_col_id + 1;
	for (A = startshellrow; A <= endshellrow; A++)
	{
		row_id_1 = basis->f_start_id[A];
		row_id_2 = basis->f_end_id[A];
		startrow = row_id_1 - start_row_id;
		nrows = row_id_2 - row_id_1 + 1;
		for (B = startshellcol; B <= endshellcol; B++)
		{
		
			col_id_1 = basis->f_start_id[B];
			col_id_2 = basis->f_end_id[B];
			startcol = col_id_1 - start_col_id;
			ncols = col_id_2 - col_id_1 + 1;
			compute_shell_pair_ovl (integrals, A, B, oed, basis);
			matrix_block_write (S, startrow, startcol, ldS,
					integrals, nrows, ncols);
		}
	}

	free (integrals);
	destroy_oed (oed);
}


void compute_H (double *H, basis_set_t * basis,
                int startshellrow, int endshellrow,
                int startshellcol, int endshellcol)
{
	oed_t *oed;
	int A;
	int B;
	int dim_1;
	double *integrals;
	int row_id_1;
	int row_id_2;
	int col_id_1;
	int col_id_2;
	int start_row_id;
	int start_col_id;
	int end_col_id;
	int ldH;
	int nrows;
	int ncols;
	int startrow;
	int startcol;
    
	oed = init_oed (basis, HHH);
    
	dim_1 = 2 * basis->max_momentum + 1;
	integrals = (double *)malloc (dim_1 * dim_1 * sizeof(double));
	assert (integrals != NULL);

	start_row_id = basis->f_start_id[startshellrow];
	start_col_id = basis->f_start_id[startshellcol];
	end_col_id = basis->f_end_id[endshellcol];
	ldH = end_col_id - start_col_id + 1;
	for (A = startshellrow; A <= endshellrow; A++)
	{
		row_id_1 = basis->f_start_id[A];
		row_id_2 = basis->f_end_id[A];
		startrow = row_id_1 - start_row_id;
		nrows = row_id_2 - row_id_1 + 1;
		for (B = startshellcol; B <= endshellcol; B++)
		{
			col_id_1 = basis->f_start_id[B];
			col_id_2 = basis->f_end_id[B];
			startcol = col_id_1 - start_col_id;
			ncols = col_id_2 - col_id_1 + 1;
			compute_shell_pair_H (integrals, A, B, oed, basis);
			matrix_block_write (H, startrow, startcol, ldH,
					integrals, nrows, ncols);
		}
	}

	free (integrals);
	destroy_oed (oed);
}

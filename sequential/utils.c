#include <string.h>
#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include <omp.h>
/* #include <mkl.h> */
/* #include <mkl_lapack.h> */
#include <cblas.h>


#include "common.h"
#include "utils.h"

void trace_matrix (int n, double *A, double *tr)
{
	int i;

	*tr = 0.0;
	for (i = 0; i < n; i++)
	{
		*tr += A[i * n + i];
	}
}


// A2 = sqrtm (A)
void sqrt_matrix (int n, double *A, double *A2, double *tmpA)
{
	int lwork;
	double *work;
	int info;
	// eigen values
	double *w;
	//eigen vectors
	double *ev;
	double lambda;
	int i;
	int j;

	DPRINTF (3, "sqrt matrix %dx%d\n", n, n);

	w = (double *)malloc (n * sizeof(double));
	ev = (double *)malloc (n * n * sizeof(double));
	tmpA = (double *)malloc (n * n * sizeof(double));
	memcpy (ev, A, n * n * sizeof(double));

	lwork = (3 * n - 1 > 1 ? 3 * n - 1 : 1);
	work = (double *)malloc (lwork * sizeof(double));

	DPRINTF (3, "  step 1 ...\n");
	dsyev_ ("V", "L", &n, ev, &n, w, work, &lwork, &info);
	assert (info == 0);

	DPRINTF (3, "  step 2 ...\n");
	for (j = 0; j < n; j++)
	{
		lambda = sqrt (w[j]);
		for (i = 0; i < n; i++)
		{
			tmpA[j * n + i] = ev[j * n + i] * lambda;
		}
	}

	DPRINTF (3, "  step 3 ...\n");

	cblas_dgemm (CblasColMajor, CblasNoTrans, CblasTrans, n, n, n,
		1.0, tmpA, n, ev, n, 0.0, A2, n);

	DPRINTF (3, "  Done\n");
}


// A2 = inv (sqrtm (A))
// dim(work) = dim (A)
void sqrtinv_matrix (int n, double *A, double *A2, double *workm)
{
	int lwork;
	double *work;
	int info;
	// eigen values
	double *w;
	//eigen vectors
	double *ev;
	double lambda;
	int i;
	int j;

	DPRINTF (3, "sqrt matrix %dx%d\n", n, n);

	w = (double *)malloc (n * sizeof(double));
	ev = (double *)malloc (n * n * sizeof(double));
	memcpy (ev, A, n * n * sizeof(double));

	lwork = (3 * n - 1 > 1 ? 3 * n - 1 : 1);
	work = (double *)malloc (lwork * sizeof(double));

	DPRINTF (3, "  step 1 ...\n");
        
	dsyev_ ("V", "L", &n, ev, &n, w, work, &lwork, &info);
	assert (info == 0);

	DPRINTF (3, "  step 2 ...\n");

	for (j = 0; j < n; j++)
	{
		lambda = 1.0 / sqrt (w[j]);
		for (i = 0; i < n; i++)
		{
			workm[j * n + i] = ev[j * n + i] * lambda;
		}
	}

	DPRINTF (3, "  step 3 ...\n");

	cblas_dgemm (CblasColMajor, CblasNoTrans, CblasTrans, n, n, n,
		1.0, workm, n, ev, n, 0.0, A2, n);

	DPRINTF (3, "  Done\n");

	free (ev);
	free (w);
	free (work);
}



void printmatRM (char *name, double *A, int nrows, int ncols)
{
	DPRINTF (1, "name[RM]: %s\n", name);
	int r;
	int c;
	for (r = 0; r < nrows; r++)
	{
		for (c = 0; c < ncols; c++)
		{
			DPRINTF (1, "%10.6e ", A[r * ncols + c]);
		}
		DPRINTF (1, "\n");
	}
}


void printmatCM (char *name, double *A, int nrows, int ncols)
{
	printf ("name[CM]: %s\n", name);
	int r;
	int c;
	for (c = 0; c < ncols; c++)
	{
		for (r = 0; r < nrows; r++)
		{
			DPRINTF (1, "%10.6e ", A[r * ncols + c]);
		}
		printf ("\n");
	}
}


void transpose (double *Ato, double *Afrom, int nrows, int ncols)
{
	int r;
	int c;
	for (r = 0; r < nrows; r++)
	{
		for (c = 0; c < ncols; c++)
		{
			Ato[r * ncols + c] = Afrom[r + nrows * c];
		}
	}
}

double **init_vecs_queue (int num_vecs, int vec_len) 
{

	double **vecs_queue;
	int i;
	vecs_queue = (double **)malloc (num_vecs * sizeof(double *));
	
	for (i = 0; i < num_vecs; i++) {
		vecs_queue[i] = (double *)malloc (vec_len * sizeof (double));
		memset (vecs_queue[i], 0, vec_len * sizeof(double));
	}

	return vecs_queue;
	
}

void destroy_vecs_queue (double **vecs_queue, int num_vecs) 
{

	int i;
	for (i = 0; i < num_vecs; i++) {
		free (vecs_queue[i]);
	}

	free (vecs_queue);

	return;
}

void clear_vecs_queue (int num_vecs, int vec_len, double **vecs_queue) 
{
	int i;

	for (i = 0; i < num_vecs; i++) {
		memset (vecs_queue[i], 0, vec_len * sizeof(double));
	}

	return;
}

void enqueue_vec (double *vec, int vec_len, double **vecs_queue)
{

}


#if 0
void initomp (int nthreads, int verbose)
{
	char schedule[1024];

	if (verbose == 1)
	{
		sprintf (schedule, "KMP_AFFINITY=granularity=fine,compact,verbose");
	}
	else
	{
		sprintf (schedule, "KMP_AFFINITY=granularity=fine,compact");
	}
	kmp_set_defaults (schedule);
	mkl_set_num_threads (nthreads);
	omp_set_num_threads (nthreads);
}
#endif


/* Xchem: An interface to the ERD integrals package. */

/* Copyright (C) Aftab Patel, Xing Liu  */

/* See ../COPYRIGHT and ../LISENCE */

#include <stdlib.h>
#include <string.h>
#include "vec_stack.h"


vec_stack *create_vec_stack (int vec_dim, int stack_len) 
{
	vec_stack *stk;
	int i;

	stk = (vec_stack *)malloc (sizeof(vec_stack));
	stk->len = stack_len;
	stk->num = 0;
	stk->vec_dim = vec_dim;
	
	stk->vecs = (double **)malloc (stk->len * sizeof(double *));

	for (i = 0; i < stk->len; i++) {
		stk->vecs[i] = (double *)malloc (stk->vec_dim * sizeof(double));
		memset (stk->vecs[i], 0, stk->vec_dim * sizeof(double));
	}
	
	return stk;
	
}

void destroy_vec_stack (vec_stack *stk)
{

	int i;
	for (i = 0; i < stk->len; i++) {
		free (stk->vecs[i]);
	}

	free (stk->vecs);
	free (stk);
	return;
}

void push_vec (double *vec, vec_stack *stk)
{
	int i;
	for (i = 0; i < stk->len - 1; i++) {
		memcpy (stk->vecs[i + 1], stk->vecs[i], stk->vec_dim * sizeof(double));
	}
	
	memcpy (stk->vecs[0], vec, stk->vec_dim * sizeof(double));
	if (stk->num < stk->len) {
		(stk->num)++;
	} 

	return;
}
	
void pop_vec (vec_stack *stk, double *vec) 
{
	int i;
	memcpy (vec, stk->vecs[0], stk->vec_dim * sizeof(double));
	for (i = 0; i < stk->len - 1; i++) {
		memcpy (stk->vecs[i], stk->vecs[i + 1], stk->vec_dim * sizeof(double));
	}
	
	if (stk->num > 0) {
		(stk->num)--;
	}
	return;

}
	
void clear_vec_stack (vec_stack *stk)
{

	int i;

	for (i = 0; i < stk->len; i++) {
		memset (stk->vecs[i], 0, stk->vec_dim * sizeof (double));
	}
	stk->num = 0;
	
	return;
}
	

	

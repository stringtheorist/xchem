#ifndef __STACK_H__
#define __STACK_H__


typedef struct vec_stack_ {
	
	double **vecs;
	int len;
	int num;
	int vec_dim;

} vec_stack;

vec_stack *create_vec_stack (int vec_dim, int stack_len);
void destroy_vec_stack (vec_stack *stk);
void push_vec (double *vec, vec_stack *stk);
void pop_vec (vec_stack *stk, double *vec);
void clear_vec_stack (vec_stack *stk);



#endif

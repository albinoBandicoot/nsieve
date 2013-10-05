#ifndef MATRIX_H
#define MATRIX_H

#include "common.h"

typedef struct {	// One of these gets associated with each row in the matrix. r2 is non-null if this row was construted by combining partials.
	rel_t *r1;	// This will have to be changed to a list should we implement the double large prime variant.
	rel_t *r2;
	uint64_t *row;	// the row of the matrix, packed into 64-bit ints.
} matrel_t;	

void solve_matrix (nsieve_t);	// the entire matrix solving and square root step. Will print out the factors (with very high probability).

#endif
#endif

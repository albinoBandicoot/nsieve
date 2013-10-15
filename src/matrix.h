#ifndef MATRIX_H
#define MATRIX_H

#include "common.h"
#include "poly.h"

void solve_matrix (nsieve_t *);	// the entire matrix solving and square root step. Will print out the factors (with very high probability).

void multiply_in_lhs (mpz_t, rel_t *, nsieve_t *);	// subroutine in the factor determination, for multiplying in a single relation.
int  construct_rhs   (uint16_t *, mpz_t, nsieve_t *);
void add_factors_to_table (uint16_t *, rel_t *);

#endif

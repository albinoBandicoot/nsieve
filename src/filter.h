#ifndef FILTER_H
#define FILTER_H

#include "common.h"
#include "sieve.h"

/* This section of the program will do the matrix building and filtering */

void build_matrix (nsieve_t *);
void combine_partials (nsieve_t *);
void filter (nsieve_t *);

#endif

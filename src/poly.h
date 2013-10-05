#ifndef POLY_H
#define POLY_H

#include <gmp.h>
#include "common.h"

#define KMAX 6
#define BMAX 1 << (KMAX - 1)	// 2**(kmax - 1)

typedef struct {
	mpz_t a;
	mpz_t b;
	mpz_t c;
	mpz_t istart;		// the start of the sieving interval. The end is istart + M. This way we can keep rel_t short by only storing an int offset from istart.
} poly_t;

typedef struct {
	mpz_t a;
	mpz_t bvals[BMAX];
	uint32_t gvals[KMAX];
	uint32_t *ainverses;	// the values of a^-1 (mod p) for each p in the factor base. This stays the same for the whole group.
} poly_group_t;

void generate_polygroup (poly_group_t *);		// this will pick some G values, compute the b values, and also precompute the inverses.
void generate_poly (poly_t *, poly_group_t *, nsieve_t *, int);	// generate the polynomial with the the i'th value of 'b' in the list in the poly_group_t. This will also compute the starting values (it needs the nsieve_t to get the square roots stored there).

void poly (poly_t *, uint32_t offset);	// evaluate the polynomial at poly->istart + offset.

#endif

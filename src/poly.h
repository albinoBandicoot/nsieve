#ifndef POLY_H
#define POLY_H

#include "common.h"

#define PG_REL_STORAGE 512

// the structures are defined in common.h

void generate_polygroup (poly_gpool_t *, poly_group_t *, nsieve_t *);		// this will pick some G values, compute the b values, and also precompute the inverses.
void generate_poly (poly_t *, poly_group_t *, nsieve_t *, int);	// generate the polynomial with the the i'th value of 'b' in the list in the poly_group_t. This will also compute the starting values (it needs the nsieve_t to get the square roots stored there).

void gpool_init (poly_gpool_t *gpool, nsieve_t *);
void advance_gpool (poly_gpool_t *, poly_group_t *);

void polygroup_init (poly_group_t *pg, nsieve_t *);
void polygroup_free (poly_group_t *pg, nsieve_t *);
void poly_init (poly_t *);
void poly_free (poly_t *);

void poly_print (poly_t *);

void poly (mpz_t res, poly_t *, int32_t offset);	// evaluate the polynomial at poly->istart + offset.

#endif

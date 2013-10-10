#ifndef SIEVE_H
#define SIEVE_H

#include "common.h"
#include "poly.h"


typedef struct {
	uint8_t sieve[BLOCKSIZE];
} block_data_t;

uint8_t fast_log (uint32_t);

void add_polygroup_relations (poly_group_t *, nsieve_t *);
void sieve_poly (block_data_t *, poly_group_t *, poly_t *, nsieve_t *);	// sieves a single polynomial completely, adding its results to relns.
void sieve_block (block_data_t *, poly_group_t *, poly_t *, nsieve_t *, int offset);	// offset is the starting offset (block# * BLOCKSIZE). 
void extract_relations (block_data_t *, poly_group_t *, poly_t *, nsieve_t *, int offset);
void construct_relation (mpz_t qx, int32_t x, poly_t *p, nsieve_t *ns);	// builds a relation and adds it to the matrix (or the hashtable if it's a partial)

int fb_factor_rel (rel_t *, uint64_t *, nsieve_t *);
int fb_factor     (mpz_t, uint64_t *, nsieve_t *);
#endif

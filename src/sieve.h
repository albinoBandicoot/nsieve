#ifndef SIEVE_H
#define SIEVE_H

#include "common.h"
#include "poly.h"

#define BLOCKSIZE 8192

typedef struct {
	uint16_t sieve[BLOCKSIZE];
} block_data_t;

uint8_t fast_log (uint32_t);

void sieve_poly (block_data_t *, poly_group_t *, poly_t *, nsieve_t *);	// sieves a single polynomial completely, adding its results to relns.
void sieve_block (block_data_t *, poly_group_t *, poly_t *, nsieve_t *, int offset);	// offset is the starting offset (block# * BLOCKSIZE). 
void extract_relations (block_data_t *, poly_group_t *, poly_t *, nsieve_t *, int offset);
#endif

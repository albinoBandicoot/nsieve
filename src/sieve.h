#ifndef SIEVE_H
#define SIEVE_H

#include "common.h"

#define BLOCKSIZE 8192

typedef struct {
	int16_t sieve[BLOCKSIZE];
} block_data_t;

void sieve_poly (block_data_t *, poly_t *, nsieve_t *);	// sieves a single polynomial completely, adding its results to relns.
void sieve_block (block_data_t *, poly_t *, nsieve_t *, int offset);	// offset is the starting offset (block# * BLOCKSIZE). 

#endif

#include "sieve.h"

void sieve_poly (block_data_t *data, poly_t *p, nsieve_t *ns){
	for (int i=0; i < ns->M/BLOCKSIZE; i++){
		sieve_block (data, p, ns, i * BLOCKSIZE);
	}
}

/* This is the real heart of the Quadratic Sieve. */
void sieve_block (block_data_t *data, poly_t *p, nsieve_t *ns, int offset){
}

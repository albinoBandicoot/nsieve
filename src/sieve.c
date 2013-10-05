#include "sieve.h"

/* Sieve an entire polynomial */
void sieve_poly (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns){
	for (int i=0; i < ns->M/BLOCKSIZE; i++){
		sieve_block (data, pg, p, ns, -ns->M + i * BLOCKSIZE);
	}
}

uint8_t fast_log (uint32_t x){	// very rough aprox. to log_2 (x)
	uint8_t res = 0;
	while (x > 0){
		x = x >> 1;
		res++;
	}
	return res - 1;
}

#define CHECK		// define this to check to make sure p | Q(x) when it should.
/* This is the real heart of the Quadratic Sieve. */
void sieve_block (block_data_t *data, poly_group_t *pg, poly_t *q, nsieve_t *ns, int block_start){
	memset (data->sieve, 0, sizeof(uint16_t) * BLOCKSIZE);	// clear the sieve
	mpz_t temp;
	mpz_init (temp);
	for (int i = 0; i < ns->fb_len; i++){
		// first compute the offsets. 
		uint32_t p = ns->fb[i];
		/* We compute start (the value for which poly(start) ~= 0 (mod p)) as 
		 *	[A^-1 * (sqrt(N) - B)] % p,  where the sqrt is the modular one we precomputed.
		 * However, we really want this relativce to the current sieve block, so 
		 * we need to adjust. If Q(x) = 0 (mod p), then compute z = (x - block_start) % p.
		 * Q(z+block_start) = 0 (mod p), and this corresponds to position z in the sieve data.
		*/
		// we should really get rid of this multiprecision code here and precompute the b % p, so then
		// everything fits in 4 bytes.
		mpz_set_ui (temp, ns->roots[i]);
		mpz_sub (temp, temp, q->b);
		mpz_mul_ui (temp, temp, pg->ainverses[i]);
		mpz_sub_ui (temp, temp, block_start);
		uint32_t z = mpz_mod_ui (temp, temp, p);

		while (z < BLOCKSIZE){
			data->sieve[z] += ns->fb_logs[i];
			z += p;
#ifdef CHECK
			poly(temp, q, z + block_start);
			mpz_mod_ui(temp, temp, p);
			if (mpz_cmp_ui(temp, 0) != 0){
				printf("Fatal error! Check failed: Q(x) mod p != 0. z = %d, block_start = %d\n", z, block_start);
				exit(1);
			}
#endif
		}

		// now do the other root (p - ns->roots[i])
		mpz_set_ui (temp, p - ns->roots[i]);
		mpz_sub (temp, temp, q->b);
		mpz_mul_ui (temp, temp, pg->ainverses[i]);
		mpz_sub_ui (temp, temp, block_start);
		z = mpz_mod_ui (temp, temp, p);

		while (z < BLOCKSIZE){
			data->sieve[z] += ns->fb_logs[i];
			z += p;
		}
	}
	extract_relations (data, pg, q, ns, block_start);
}

void extract_relations (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns, int block_start){
}


#include "sieve.h"

/* Sieve an entire polynomial */
void sieve_poly (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns){
	for (int i=0; i < 2*p->M/BLOCKSIZE; i++){
		sieve_block (data, pg, p, ns, -p->M + i * BLOCKSIZE);
		ns->sieve_locs += BLOCKSIZE;
	}
//	free (p->bmodp);	// free these temporary precomputed values.
//	free (p->offsets1);
//	free (p->offsets2);
}


/* This gets called after all of the polynomials in a block have been sieved to collect the relations together, and do the multiplying
 * through by the victim. */

void add_polygroup_relations (poly_group_t *pg, nsieve_t *ns){
	// now we get to pick our victim. It must be a full relation (though I guess theoretically if there were no fulls but two
	// partials from this poly group shared a cofactor it could be used, but that's more thinking and debugging than it's worth)
	// Theoretically, it would probably be best to pick the sparsest relation as the victim, but we'll just pick the first one.
	for (int i=0; i < pg->nrels; i++){
//		printf("Looking at i = %d; cofactor = %d\n", i, pg->relns[i]->cofactor);
		if (pg -> relns[i]->cofactor == 1){
			pg->victim = pg->relns[i];
			break;
		}
	}

	if (pg->victim != NULL){	// we found a full relation
		// factor victim over the fb:
		mpz_t temp;
		mpz_init (temp);
		fb_factor_rel(pg->victim, pg->victim_factors, ns);
		mpz_clear (temp);
		// now we will loop through the relations. For each full relation, we construct a matrel_t for it, and fill in the 
		// row (multiplying (xor-ing) by the victim). The partials will have this done to them whenever they are added to
		// the matrix at the end of sieving. 
		for (int i=0; i < pg->nrels; i++){
			if (ns->nfull >= ns->rels_needed) return;	// we're done sieving.
			if (pg->relns[i] == pg->victim) continue;	// we don't want to add the victim to the list.
			if (pg->relns[i]->cofactor == 1){	// full relation
				matrel_t *m = &ns->relns[ns->nfull];
				m -> r1 = pg->relns[i];
				m -> r2 = NULL;
				fb_factor_rel (m->r1, m->row, ns);
				xor_row (m->row, pg->victim_factors, ns->row_len);
				ns->nfull ++;
			}
		}
	} else {
		// we did not find one. This is not good, but not an error either - we were just unlucky. However, we should probably be
		// doing either larger sieve intervals or a larger k or something. 
		printf("There are no full relations for this polygroup! We must throw away the partials.\n");
	}
}

uint8_t fast_log (uint32_t x){	// very rough aprox. to log_2 (x)
	x = (x*3)/2;		// we do this so that some of the logs of the primes are too high and some are too low. If we left x as is, 
				// they would be consistently too low, and the errors would accumulate, proportionally to the number of factors found.
				// The idea is that the errors will now be essentially random, and average to roughly 0. 
				// We are using 3/2 ~= sqrt(2).
	uint8_t res = 0;
	while (x > 0){
		x = x >> 1;
		res++;
	}
	return res - 1;
}

void mpz_sub_si (mpz_t res, mpz_t a, int32_t x){
	if (x < 0){
		mpz_add_ui(res, a, -x);
	} else {
		mpz_sub_ui(res, a, x);
	}
}

inline uint32_t get_offset (uint32_t p, int32_t i, int block_start, int rn, poly_t *q, poly_group_t *pg, nsieve_t *ns){
	int64_t z = q->bmodp[i];
	if (rn == 0){
		z += ns->roots[i];
	} else {
		z += p - ns->roots[i];
	}
	z *= pg->ainverses[i];
	z -= block_start;
	if (z < 0){
		z = p + (z % p);
		if (z == p) z = 0;
	} else {
		z %= p;
	}
	return (uint32_t) z;
}

//#define CHECK		// define this to check to make sure p | Q(x) when it should. It slows things down a LOT, and should only be used for debugging.
/* This is the real heart of the Quadratic Sieve. */
void sieve_block (block_data_t *data, poly_group_t *pg, poly_t *q, nsieve_t *ns, int block_start){
	memset (data->sieve, 0, sizeof(uint8_t) * BLOCKSIZE);	// clear the sieve
	mpz_t temp;
	mpz_init (temp);
	for (int i = 25; i < ns->fb_len; i++){
		// first compute the offsets. 
		if (ns->fb[i] <= pg->gvals[ns->k-1] && ns->fb[i] >= pg->gvals[0]){
			while (ns->fb[i] <= pg->gvals[ns->k-1]){
				i++;
			}
		}
		uint32_t p = ns->fb[i];
		/* We compute start (the value for which poly(start) ~= 0 (mod p)) as 
		 *	[A^-1 * (sqrt(N) - B)] % p,  where the sqrt is the modular one we precomputed.
		 * However, we really want this relativce to the current sieve block, so 
		 * we need to adjust. If Q(x) = 0 (mod p), then compute z = (x - block_start) % p.
		 * Q(z+block_start) = 0 (mod p), and this corresponds to position z in the sieve data.
		*/

		// first we check to see if p is one of the g_i; if so, the sieve has the potential to break on those values.
		// In fact, p will either divide none of the values, or all of them.
		/*
		int bad = (p >= pg->gvals[0] && p <
		for (int g=0; g < ns->k; g++){
			if (p == pg->gvals[g]){
				bad = 1;
				break;
			}
		}
		if (bad == 1) continue;
		*/
		// we should really get rid of this multiprecision code here and precompute the b % p, so then
		// everything fits in 4 bytes.
//		uint32_t z = mod(q->offsets1[i] - (q->M + block_start),  p);
		uint32_t z = get_offset (p, i, block_start, 0, q, pg, ns);
		/*
		int64_t z = ns->roots[i] + q->bmodp[i];	// this is a + here because we actually precomputed -b % p.
							// this has to be 64 bits, as the intermediate values could get above 2^32 easily.
		z *= pg->ainverses[i];
		z -= block_start;		// there is the possibility that z is now negative. We must then be careful with the mod.
		if (z < 0){
			z = p + (z % p);	// because C's mod function doesn't actually do a mathematicla mod for negative numbers.
		} else {
			z %= p;
		}
		*/
		/*
		mpz_set_ui (temp, ns->roots[i]);
		mpz_sub (temp, temp, q->b);
		mpz_mul_ui (temp, temp, pg->ainverses[i]);
		mpz_sub_si (temp, temp, block_start);
		uint32_t z = mpz_mod_ui (temp, temp, p);
		*/

		while (z < BLOCKSIZE){
			data->sieve[z] += ns->fb_logs[i];
#ifdef CHECK
			poly(temp, q, z + block_start);
			mpz_mod_ui(temp, temp, p);
			if (mpz_cmp_ui(temp, 0) != 0){
				printf("Fatal error! Check failed: Q(x) mod %d != 0. z = %d, block_start = %d\n", p, z, block_start);
				exit(1);
			}
#endif
			z += p;
		}

		// now do the other root (p - ns->roots[i])
		if (p == 2) continue;	// there's only one root of 2, so continuing would erroneously duplicate the sieve on p=2.
//		z = mod(q->offsets2[i] - (q->M + block_start), p); 
		z = get_offset (p, i, block_start, 1, q, pg, ns);
/*
		z = p - ns->roots[i] + q->bmodp[i];
		z *= pg->ainverses[i];
		z -= block_start;
		if (z < 0){
			z = p + (z % p);
		} else {
			z %= p;
		}
		*/
		/*
		mpz_set_ui (temp, p - ns->roots[i]);
		mpz_sub (temp, temp, q->b);
		mpz_mul_ui (temp, temp, pg->ainverses[i]);
		mpz_sub_si (temp, temp, block_start);
		z = mpz_mod_ui (temp, temp, p);
		*/
		while (z < BLOCKSIZE){
			data->sieve[z] += ns->fb_logs[i];
#ifdef CHECK
			poly(temp, q, z + block_start);
			mpz_mod_ui(temp, temp, p);
			if (mpz_cmp_ui(temp, 0) != 0){
				printf("Fatal error! Check failed: Q(x) mod %d != 0. z = %d, block_start = %d\n", p, z, block_start);
				exit(1);
			}
#endif
			z += p;
		}
	}
	extract_relations (data, pg, q, ns, block_start);
}

void extract_relations (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns, int block_start){
	/* We want to estimate the size of the unfactored part of each number in the sieve, which is about
	 * r = log(Q(x+block_start)) - sieve[x]. If this is close to 0, we fully factored the relation. If we
	 * let F = log(fb[fb_len-1]), and if F < r < F^2, then the remaining factor is prime, and we have a partial.
	 */
//	int cutoff = (int) (fast_log(ns->fb_bound) * 1.4);	// for now, this is very ad hoc	
	int cutoff = (int) (fast_log(ns->lp_bound) * ns->T);
	mpz_t temp;
	mpz_init(temp);
	poly(temp, p, block_start + BLOCKSIZE/2);
	mpz_abs(temp, temp);
	int logQ = mpz_sizeinbase (temp, 2);
//	printf("poly(middle) = ");
//	mpz_out_str(stdout, 10, temp);
//	printf("  has %d digits in base 2\n", logQ);
#define CHUNK_SCANNING	
#ifdef CHUNK_SCANNING
	uint64_t *chunk = (uint64_t *) data->sieve;
	uint32_t nchunks = BLOCKSIZE/8;
	uint8_t maskchar = 1;
	while (maskchar < logQ - cutoff){
		maskchar *= 2;
	}
	maskchar /= 2;
	maskchar = ~ (maskchar - 1);
	uint64_t mask = maskchar;
	for (int i=0; i<8; i++){
		mask = (mask << 8) | maskchar;
	}
//	printf("logq = %d, cutoff = %d, maskchar = %d, mask = %llx\n", logQ, cutoff, (int) maskchar, mask);
	// we have constructed our mask.
	for (int i=0; i < nchunks; i++){
		if (mask & chunk[i]){
			for (int j=0; j<8; j++){
				if (logQ - data->sieve[i*8+j] < cutoff){
					poly (temp, p, block_start + i*8 + j);
					construct_relation (temp, block_start + i*8 + j, p, ns);
				}
			}
		}
	}
#else
	for (int i = 0; i < BLOCKSIZE; i++){
//		printf("x = %d (i = %d): log is %d; logQ = %d, thresh = %d\n", (i+block_start), i, data->sieve[i], logQ, cutoff);
		if (logQ - data->sieve[i] < cutoff){
			// we probably have at least a partial relation.
			poly (temp, p, block_start + i);
			construct_relation (temp, i+block_start, p, ns);
		}
	}
#endif
}

/* NOTE WELL: This routine assumes x factors over the factor base. Cofactors are divided out ahead of time */
void fb_factor_rel (rel_t *rel, uint64_t *row, nsieve_t *ns){
	mpz_t x;
	mpz_init (x);
	poly (x, rel->poly, rel->x);
	mpz_divexact_ui(x, x, rel->cofactor);

	clear_row (row, ns);
	if (mpz_cmp_ui(x, 0) < 0){
		mpz_neg (x,x);
		flip_bit (row, 0);
	}
	uint64_t q;
	int i;
	while (mpz_divisible_ui_p(x, 2)){
		mpz_divexact_ui(x, x, 2);
		flip_bit (row, 1);
	}
	for (i=1; i < ns->fb_len; i++){
		if (ns->fb[i] >= rel->poly->group->gvals[0] && ns->fb[i] <= rel->poly->group->gvals[ns->k-1]) goto mainloop;	// we don't know if the values in this range will divide or not, so we go to the standard loop without the initial divide.
		if (get_offset (ns->fb[i], i, rel->x, 0, rel->poly, rel->poly->group, ns) == 0 || get_offset (ns->fb[i], i, rel->x,  1, rel->poly, rel->poly->group, ns) == 0){
			mpz_divexact_ui (x, x, ns->fb[i]);
			flip_bit (row, i+1);
			if (mpz_fits_64 (x)) goto fp_tdiv;
mainloop:
			while (mpz_divisible_ui_p(x, ns->fb[i])){
				mpz_divexact_ui(x, x, ns->fb[i]);
				flip_bit (row, i+1);
				if (mpz_fits_64 (x)){
					goto fp_tdiv;
				}
			}
		}
	}
fp_tdiv:
	q = mpz_get_64 (x);
	mpz_clear(x);
	while (i < ns->fb_len){
		while (q % ns->fb[i] == 0){
			q /= ns->fb[i];
			flip_bit (row, i+1);
			if (q < ns->fb[i] * ns->fb[i]){
				if (q == 1) return;
				if (q < ns->fb_bound){	// q is a prime in the factor base.
					// here we do binary search to find the index of q.
					int high = ns->fb_len - 1;
					int low = i;
					int guess = (high + low) / 2;
					if (q < ns->fb[low] || q > ns->fb[high]){
						printf("Q OUT OF RANGE AAAAH!\n");
					}
					while (high - low > 1){
						if (ns->fb[guess] < q){
							low = guess;
							guess = (high + low) / 2;
						} else if (ns->fb[guess] > q){
							high = guess;
							guess = (high + low) / 2;
						} else {
							flip_bit (row, guess + 1);
							return;
						}
					}
					for (guess = low; guess <= high; guess++){
						if (q == ns->fb[guess]){
							flip_bit (row, guess + 1);
							return;
						}
					}
					printf("UH OH BAD BAD.\n");
					return;
				}
			}
		}
		i++;
	}
}


/* Allocates the rel_t object, and adds it to the list in the polygroup if it factored or was a partial.
 * The trial division here is solely to determine if it actually factors; it does not fill in the matrix rows.
*/

void construct_relation (mpz_t qx, int32_t x, poly_t *p, nsieve_t *ns){
	ns->tdiv_ct ++;
	rel_t *rel = (rel_t *)(malloc(sizeof(rel_t)));
	rel->poly = p;
	rel->x = x;
	mpz_abs(qx, qx);
	uint64_t q = 0;
	int i;
	while (mpz_divisible_ui_p(qx, 2)){
		mpz_divexact_ui(qx, qx, 2);
	}
	for (i=1; i < ns->fb_len; i++){
		if (get_offset (ns->fb[i], i, x, 0, p, p->group, ns) == 0 || get_offset (ns->fb[i], i, x, 1, p, p->group, ns) == 0){
//		if ((p->offsets1[i] -(p->M + x)) % ns->fb[i] == 0 || (p->offsets2[i] -(p->M + x)) % ns->fb[i] == 0){
			mpz_divexact_ui(qx, qx, ns->fb[i]);
			if (mpz_fits_64 (qx)) goto fixedprec_tdiv;
			while (mpz_divisible_ui_p (qx, ns->fb[i])){
				mpz_divexact_ui(qx, qx, ns->fb[i]);
				if (mpz_fits_64 (qx)){	// we don't need to test for anything else here - all the relevant cases happen only when qx is well under 2^64.
					goto fixedprec_tdiv;
				}
			}
		}
	}
fixedprec_tdiv:
	q = mpz_get_64 (qx);
	if (q < ns->fb[i] * ns->fb[i]){
		if (q < ns->fb_bound) goto add_full_rel;
		if (q < ns->lp_bound) goto add_partial_rel;
	}
	while (i < ns->fb_len){
		while (q % ns->fb[i] == 0){
			q /= ns->fb[i];
			if (q < ns->fb[i] * ns->fb[i]){
				if (q < ns->fb_bound) goto add_full_rel;
				if (q < ns->lp_bound) goto add_partial_rel;
			}
		}
		i++;
	}
	// if we're here, we weren't able to do anything with this relation.
	free (rel);
	return;

add_full_rel:
	rel->cofactor = 1;
	p->group->relns[ p->group->nrels ] = rel;
	p->group->nrels ++;
	return;
add_partial_rel:
	rel->cofactor = q;
	ht_add (&ns->partials, rel);
	return;

}

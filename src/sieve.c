#include "sieve.h"

/* Sieve an entire polynomial */
void sieve_poly (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns){
	for (int i=0; i < 2*p->M/BLOCKSIZE; i++){
		sieve_block (data, pg, p, ns, -p->M + i * BLOCKSIZE);
	}
}


/* This gets called after all of the polynomials in a block have been sieved to collect the relations together, and do the multiplying
 * through by the victim. */

void add_polygroup_relations (poly_group_t *pg, nsieve_t *ns){
	// now we get to pick our victim. It must be a full relation (though I guess theoretically if there were no fulls but two
	// partials from this poly group shared a cofactor it could be used, but that's more thinking and debugging than it's worth)
	// Theoretically, it would probably be best to pick the sparsest relation as the victim, but we'll just pick the first one.
	for (int i=0; i < pg->nrels; i++){
		printf("Looking at i = %d; cofactor = %d\n", i, pg->relns[i]->cofactor);
		if (pg -> relns[i]->cofactor == 1){
			pg->victim = pg->relns[i];
			break;
		}
	}

	if (pg->victim != NULL){	// we found a full relation
		// factor victim over the fb:
		mpz_t temp;
		mpz_init (temp);
		poly (temp, pg->victim->poly, pg->victim->x);
		fb_factor (temp, pg->victim_factors, ns);
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

//#define CHECK		// define this to check to make sure p | Q(x) when it should. It slows things down a LOT, and should only be used for debugging.
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

		// first we check to see if p is one of the g_i; if so, the sieve has the potential to break on those values.
		// In fact, p will either divide none of the values, or all of them.
		int bad = 0;
		for (int g=0; g < ns->k; g++){
			if (p == pg->gvals[g]){
				bad = 1;
				break;
			}
		}
		if (bad == 1) continue;
		// we should really get rid of this multiprecision code here and precompute the b % p, so then
		// everything fits in 4 bytes.
		mpz_set_ui (temp, ns->roots[i]);
		mpz_sub (temp, temp, q->b);
		mpz_mul_ui (temp, temp, pg->ainverses[i]);
		mpz_sub_si (temp, temp, block_start);
		uint32_t z = mpz_mod_ui (temp, temp, p);

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
		mpz_set_ui (temp, p - ns->roots[i]);
		mpz_sub (temp, temp, q->b);
		mpz_mul_ui (temp, temp, pg->ainverses[i]);
		mpz_sub_si (temp, temp, block_start);
		z = mpz_mod_ui (temp, temp, p);

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
	
	for (int i = 0; i < BLOCKSIZE; i++){
//		printf("x = %d (i = %d): log is %d; logQ = %d, thresh = %d\n", (i+block_start), i, data->sieve[i], logQ, cutoff);
		if (logQ - data->sieve[i] < cutoff){
			// we probably have a full relation.
			poly (temp, p, block_start + i);
	//		mpz_mul(temp, temp, p->a);	// I think I need to do this.
	//		printf("Passed cutoff test: x = %d; logQ - sieve[i] = %d; Q(x) = ", i + block_start, logQ - data->sieve[i]);
	//		mpz_out_str(stdout, 10, temp);
	//		printf("\n");
			if (ns->nfull >= ns->rels_needed) return;
			construct_relation (temp, i+block_start, p, ns);
		}
	}
}

int fb_factor_rel (rel_t *rel, uint64_t *row, nsieve_t *ns){
	mpz_t x;
	mpz_init (x);
	poly (x, rel->poly, rel->x);
//	mpz_divexact_ui(x, x, rel->cofactor);
	int rval = fb_factor (x, row, ns);
	mpz_clear (x);
	return rval;
}

int fb_factor (mpz_t x, uint64_t *row, nsieve_t *ns){	// returns 0 if factors completely. At the end, x contains the unfactored portion.
	clear_row (row, ns);
	if (mpz_cmp_ui(x, 0) < 0){
		mpz_neg (x,x);
		flip_bit (row, 0);
	}
	for (int i=0; i < ns->fb_len; i++){
		while (mpz_divisible_ui_p(x, ns->fb[i])){
			mpz_divexact_ui(x, x, ns->fb[i]);
			flip_bit (row, i+1);
			if (mpz_cmp_ui(x, 1) == 0){
				return 0;
			}
		}
	}
	return -1;
}


/* Allocates the rel_t object, and adds it to the list in the polygroup if it factored or was a partial.
 * The trial division here is solely to determine if it actually factors; it does not fill in the matrix rows.
*/
void construct_relation (mpz_t qx, int32_t x, poly_t *p, nsieve_t *ns){
	rel_t *rel = (rel_t *)(malloc(sizeof(rel_t)));
	rel->poly = p;
	rel->x = x;
	mpz_abs(qx, qx);
	for (int i=0; i < ns->fb_len; i++){
		while (mpz_divisible_ui_p (qx, ns->fb[i])){
			mpz_divexact_ui(qx, qx, ns->fb[i]);
			if (mpz_cmp_ui(qx, 1) == 0){
				break;
			}
		}
	}
	if (mpz_cmp_ui(qx, 1) == 0){
		if (p->group->nrels < PG_REL_STORAGE){
//			printf("Adding full relation; x = %d to position %d\n", x, p->group->nrels);
			rel->cofactor = 1;
			p->group->relns[ p->group->nrels ] = rel;
			p->group->nrels ++;
			return;
		}
	} else {
		if (mpz_cmp_ui(qx, ns->lp_bound) < 0){
			if (mpz_probab_prime_p(qx, 10)){	// if lp_bound < fb_bound^2, this primality test could be omitted.
				/*
				rel->cofactor = mpz_get_ui(qx);
				if (p->group->nrels < PG_REL_STORAGE){
					p->group->relns[ p->group->nrels ++ ] = rel;
					return;
				}
				*/
		//		ht_add (&ns->partials, rel);
			}
		}
	}
	// if we're here, we weren't able to do anything with this relation.
	free (rel);
}

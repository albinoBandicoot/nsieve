#include "sieve.h"

/* Sieve an entire polynomial */
void sieve_poly (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns){
	int start = (p->M * BLOCKSIZE / 2);
	start = -start;
	int i = 0;
	for (i=0; i < p->M; i++){
		sieve_block (data, pg, p, ns, start + i * BLOCKSIZE);
		ns->sieve_locs += BLOCKSIZE;
	}
	free (p->bmodp);	// free these temporary precomputed values.
}


/* This gets called after all of the polynomials in a block have been sieved to collect the relations 
 * together and do the multiplying through by the victim. */

void add_polygroup_relations (poly_group_t *pg, nsieve_t *ns){
	/* now we get to pick our victim. It must be a full relation (though I guess theoretically if 
	 * there were no fulls but two partials from this poly group shared a cofactor it could be used, 
	 * but that's more thinking and debugging than it's worth, especially since if we're finding 
	 * groups without full relations, our parameters are way off.
	 *
	 * Theoretically, it would probably be best to pick the sparsest relation as the victim, but 
	 * we'll just pick the first one.
	*/
	for (int i=0; i < pg->nrels; i++){
		if (pg -> relns[i]->cofactor == 1){
			if (fl_check (pg->relns[i], ns)){
				pg->victim = pg->relns[i];
				break;
			}
		}
	}

	/* Since we will me modifying state in the nsieve_t, which is the same for all threads, 
	 * we need to acquire the lock first to make sure we don't have multiple threads writing
	 * all over each other. */
	pthread_mutex_lock (&ns->lock);

	if (pg->victim != NULL){	// we found a full relation
		for (int i=0; i < pg->nrels; i++){
			if (ns->nfull >= ns->rels_needed) {	// we're done sieving.
				ns->info_npg ++;
				ns->info_npoly += ns->bvals;
				pthread_mutex_unlock (&ns->lock);
				return;	
			}
			if (pg->relns[i] == pg->victim) continue;	// we don't want to add the victim to the list.
			if (!fl_check (pg->relns[i], ns)) continue;	// don't add bad relations
			if (pg->relns[i]->cofactor == 1){		// full relation
				matrel_t *m = &ns->relns[ns->nfull];
				m -> r1 = pg->relns[i];
				m -> r2 = NULL;
				fl_concat  (m->r1, pg->victim);
				ns->nfull ++;
			} else {
				ht_add (&ns->partials, pg->relns[i]);
			}
		}
	} else {
		// we did not find one. This is not good, but not an error either - we were just unlucky. 
		// However, we should probably be doing either larger sieve intervals or a larger k or something. 
		printf("There are no full relations for this polygroup! We must throw away the partials.\n");
	}
	ns->info_npg ++;
	ns->info_npoly += ns->bvals;
	pthread_mutex_unlock (&ns->lock);
}

/* This is a very cheap and somewhat imprecise approximation to log_2 (x). It is used to compute the
 * logs of each factor base prime. It is these values that will be added up in the sieve bins. */
uint8_t fast_log (uint32_t x){
	/* First, multiply x by about sqrt(2), so that our truncation will not produce values that
	 * have systematic error in them. As effectively random values add up (from the various primes
	 * that are being multiplied together), the errors will average to 0.
	*/
	x = (x*3)/2;
	uint8_t res = 0;
	while (x > 0){
		x = x >> 1;
		res++;
	}
	return res - 1;
}

/* Find the offset relative to block_start of the first value x_0 of the polynomial q that has the
 * property that p | q(x_0). This computation is done with the aid of the precomputed values ns->roots,
 * pg->ainverses, and q->bmodp. Acceleration of this routine would be most welcome, since it accounts
 * for perhaps a quarter or a third of the sieve time (note its use also in the trial division code)
 *
 * We compute start (the value for which poly(start) ~= 0 (mod p)) as 
 *	[A^-1 * (sqrt(N) - B)] % p,  where the sqrt is the modular one we precomputed.
 * However, we really want this relativce to the current sieve block, so we need to adjust. If 
 * Q(x) = 0 (mod p), then compute z = (x - block_start) % p. Q(z+block_start) = 0 (mod p), and 
 * this corresponds to position z in the sieve data.
*/

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

/* This is the real heart of the Quadratic Sieve. */
void sieve_block (block_data_t *data, poly_group_t *pg, poly_t *q, nsieve_t *ns, int block_start){
	
	memset (data->sieve, 0, sizeof(uint8_t) * BLOCKSIZE);	// initialize the sieve

	/* Notice that we don't start sieving with the first prime in the FB; instead with start
	 * with the 25th (a somewhat arbitrary choice). This is because sieving takes time on the 
	 * order of 1/p, since only 2/p sieve locations will be divisible by p. Hence, the loop will
	 * run on average BLOCKSIZE/p times for each block. Thus the smallest few primes will take the
	 * majority of the time. However, they also contribute the least to the sieve values, so it
	 * makes sense to skip them entirely.
	 *
	 * Do we lose relations this way? Absolutely, but fewer than you might expect. It turns out 
	 * that the greatly increased sieving speed more than makes up for it. 
	*/
	for (int i = 25; i < ns->fb_len; i++){
		/* This checks to see whether we would be trying to sieve with the primes g that make up A.
		 * If we are, skip until we're beyond that range. */
		if (ns->fb[i] <= pg->gvals[ns->k-1] && ns->fb[i] >= pg->gvals[0]){
			while (ns->fb[i] <= pg->gvals[ns->k-1]){
				i++;
			}
		}

		uint32_t p = ns->fb[i];
		uint32_t z = get_offset (p, i, block_start, 0, q, pg, ns);
		// p | q(z+block_start)

		while (z < BLOCKSIZE){
			data->sieve[z] += ns->fb_logs[i];	// we keep in each sieve bin a running
			z += p;					// total of the logs of each prime that
					// divided it; if this is close to 0 at the end, we probably
					// have a useable relation.
		}

		// now do the other root (p - ns->roots[i]), since there will be 2 modular square roots
		// for each prime. Note that 2 is an exception, but we aren't sieving with 2 anyway.
		z = get_offset (p, i, block_start, 1, q, pg, ns);
		while (z < BLOCKSIZE){
			data->sieve[z] += ns->fb_logs[i];
			z += p;
		}
	}
	extract_relations (data, pg, q, ns, block_start);
}

/* This routine will scan over the sieve for values that look promising, and then attempt to factor
 * them over the factor base. If it finds one that factors (or is a partial), it will add it to the
 * list that we're accumulating for this poly_group.
*/
void extract_relations (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns, int block_start){
	/* The sieve now contains estimates of the size (in bits) of the unfactored portion of the
	 * polynomial values. We scan the sieve for values less than this cutoff, and trial divide
	 * the ones that pass the test. */

	mpz_t temp;
	mpz_init (temp);
	poly(temp, p, block_start + BLOCKSIZE/2);
	mpz_abs(temp, temp);
	uint8_t logQ = (uint8_t) mpz_sizeinbase (temp, 2);

	int cutoff = (int) (fast_log(ns->lp_bound) * ns->T);

	/* To accelerate the sieve scanning for promising values, instead of comparing each 8-bit entry
	 * one at a time, we instead look at 64 bits at a time. We are looking for sieve values such that
	 * sieve[x] > logQ - cutoff. We can use as a first approximation to this the test that 
	 * sieve[x] > S, where S is the smallest power of 2 larger than (logQ - cutoff). Create an 8-bit
	 * mask - say S is 32 - the mask is 11100000. Then if you AND the mask with any sieve value, and 
	 * the result is nonzero, it cannot possibly be below the cutoff. Make 8 copies of the mask, put
	 * it in a 64-bit int, cast the pointer to the sieve block, and do this AND on 64-bit chunks (8
	 * sieve locations at a time). Since it will be relatively rare for a value to pass the cutoff, 
	 * most of the time this test will immediately reject all 8 values. If it doesn't, test each one
	 * against the cutoff individually.
	*/
	
	uint64_t *chunk = (uint64_t *) data->sieve;
	uint32_t nchunks = BLOCKSIZE/8;
	uint8_t maskchar = 1;		// this is the 8-bit version of the mask
	while (maskchar < logQ - cutoff){
		maskchar *= 2;
	}
	maskchar /= 2;
	maskchar = ~ (maskchar - 1);
	uint64_t mask = maskchar;
	for (int i=0; i<8; i++){	// make 8 copies of it
		mask = (mask << 8) | maskchar;
	}

	// now loop over the sieve, a chunk at a time
	for (int i=0; i < nchunks; i++){
		if (mask & chunk[i]){	// then some value *might* have passed the test
			for (int j=0; j<8; j++){	
				if (logQ - data->sieve[i*8+j] < cutoff){	// check them all
					poly (temp, p, block_start + i*8 + j);
					construct_relation (temp, block_start + i*8 + j, p, ns);
				}
			}
		}
	}
	mpz_clear (temp);
}

/* Allocates the rel_t object, and adds it to the list in the polygroup if it factored or was a partial.
 * It determines the factors by trial division, and it also adds the factors to the linked list in the
 * rel_t. If it didn't factor and wasn't a partial, the relation is freed.
*/
void construct_relation (mpz_t qx, int32_t x, poly_t *p, nsieve_t *ns){
	ns->tdiv_ct ++;
	rel_t *rel = (rel_t *)(malloc(sizeof(rel_t)));
	if (rel == NULL){
		printf ("Malloc failed\n");
		exit(1);
	}
	rel->poly = p;
	rel->x = x;
	rel->cofactor = 1;
	rel->factors = NULL;
	if (mpz_cmp_ui (qx, 0) < 0){
		fl_add (rel, 0);
	}
	mpz_abs(qx, qx);
	uint64_t q = 0;
	int i;
	while (mpz_divisible_ui_p(qx, 2)){	// handle 2 separately
		mpz_divexact_ui(qx, qx, 2);
		fl_add (rel, 1);
	}
	for (i=1; i < ns->fb_len; i++){
		// instead of doing a multi-precision divisiblilty test, we can use the get_offset method to
		// detect if 'x' is in the arithmetic progression of sieve values divisible by ns->fb[i].
		if (get_offset (ns->fb[i], i, x, 0, p, p->group, ns) == 0 || get_offset (ns->fb[i], i, x, 1, p, p->group, ns) == 0){
			mpz_divexact_ui(qx, qx, ns->fb[i]);
			fl_add (rel, i+1);	// add to the factor list
			// If the result of the division fits in 64 bits, ditch the arbitrary precision.
			if (mpz_fits_64 (qx)) goto fixedprec_tdiv;
			while (mpz_divisible_ui_p (qx, ns->fb[i])){	// the sieve doesn't tell us
				mpz_divexact_ui(qx, qx, ns->fb[i]);	// how many times the factor divided
				fl_add (rel, i+1);
				if (mpz_fits_64 (qx)){
					goto fixedprec_tdiv;
				}
			}
		}
	}
fixedprec_tdiv:
	q = mpz_get_64 (qx);
	if (q < ns->fb[i] * ns->fb[i]){	// q must be prime
		if (q < ns->fb_bound){	// if it's less than the factor base bound, it had better be in the FB.
			fl_add (rel, fb_lookup (q, ns));	// look it up and add it to the list.
			goto add_rel;
		}
		if (q < ns->lp_bound) {	// in this case we have a partial relation.
			rel->cofactor = q;
			goto add_rel;
		}
	}
	while (i < ns->fb_len){		// continue the trial division
		while (q % ns->fb[i] == 0){	// it is no longer efficient to compute offsets here (that 
			q /= ns->fb[i];		// calculation involved mods!)
			fl_add (rel, i+1);
			if (q < ns->fb[i] * ns->fb[i]){
				if (q < ns->fb_bound){
					fl_add (rel, fb_lookup (q, ns));
					goto add_rel;
				}
				if (q < ns->lp_bound) {
					rel->cofactor = q;
					goto add_rel;
				}
			}
		}
		i++;
	}
	// if we're here, we weren't able to do anything with this relation.
//	rel_free (rel);
	return;

add_rel:	// add the relation to the list in the poly_group_t we're working with.
	if (p->group->nrels < PG_REL_STORAGE){
		p->group->relns[ p->group->nrels ] = rel;
		p->group->nrels ++;
		return;
	} else {	// if we ran out of space, just let it go. 
//		rel_free(rel);
		return;
	}
}

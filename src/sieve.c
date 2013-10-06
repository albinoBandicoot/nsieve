#include "sieve.h"

/* Sieve an entire polynomial */
void sieve_poly (block_data_t *data, poly_group_t *pg, poly_t *p, nsieve_t *ns){
	for (int i=0; i < 2*p->M/BLOCKSIZE; i++){
		sieve_block (data, pg, p, ns, -p->M + i * BLOCKSIZE);
		if (ns->nfull >= ns->rels_needed) return;
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

void mpz_sub_si (mpz_t res, mpz_t a, int32_t x){
	if (x < 0){
		mpz_add_ui(res, a, -x);
	} else {
		mpz_sub_ui(res, a, x);
	}
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
	int cutoff = (int) (fast_log(ns->lp_bound) * 1.4);
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

/* Allocates the rel_t and matrel_t objects, and adds the matrel_t to the list of full relations.
 * Performs the trial division, producing the row of the matrix. It may be revealed by doing this that
 * the relation in fact does not factor over the factor base. If it is a partial, it is instead added
 * to the hashtable.
*/
void construct_relation (mpz_t qx, int32_t x, poly_t *p, nsieve_t *ns){
	rel_t *rel = (rel_t *)(malloc(sizeof(rel_t)));
	rel->poly = p;
	rel->x = x;
	clear_row (&ns->relns[ns->nfull], ns);
	if (mpz_cmp_ui(qx, 0) < 0){
		mpz_abs(qx, qx);
		flip_bit (&ns->relns[ns->nfull], 0);
	}
#ifdef FULL_TDIV
	long t = 2;
	while (t < ns->fb_bound){
		while (mpz_divisible_ui_p(qx, t)){
			mpz_divexact_ui(qx, qx, t);
		}
		if (t == 2){
			t ++;
		} else {
			t += 2;
		}
		if (mpz_cmp_ui(qx, 1) == 0){
			break;
		}
	}
#endif
#ifndef FULL_TDIV
	for (int i=0; i < ns->fb_len; i++){
		while (mpz_divisible_ui_p (qx, ns->fb[i])){
			mpz_divexact_ui(qx, qx, ns->fb[i]);
			flip_bit (&ns->relns[ns->nfull], i+1);
		}
		if (mpz_cmp_ui(qx, 1) == 0){
			break;
		}
	}
#endif 
	if (mpz_cmp_ui(qx, 1) == 0){
		printf("Found full relation; x = %d\n", x);
		ns->relns[ns->nfull].r1 = rel;
		ns->relns[ns->nfull].r2 = NULL;
		rel->cofactor = 1;
		ns->nfull ++;
	} else {
		if (mpz_cmp_ui(qx, ns->lp_bound) < 0){
			if (mpz_probab_prime_p(qx, 10)){	// if lp_bound < fb_bound^2, this primality test could be omitted.
				ns->npartial++;
	//			printf("Found partial relation; x = %d, cofactor = %ld\n", x, mpz_get_ui(qx));

	//			rel->cofactor = mpz_get_ui(qx);
	//			ht_add (&ns->partials, rel);
	//			return;
			}
		}
	}
	// if we're here, we weren't able to do anything with this relation.
	free (rel);
}

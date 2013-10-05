#include "poly.h"

void generate_polygroup (poly_group_t *pg, nsieve_t *ns){
	/* This is a tricky one. First we must choose A, by picking k primes g_i, and multiplying them together.
	 * Then we need to find all of the values of B which satisfy  B^2 = N (mod a). There will be 2^(k-1) of them. 
	 * Then we will compute the values of A^-1 (mod p) for each p in the factor base. This is really a precomputation
	 * to speed up the computation of Q(x) = 0 (mod p). 
	*/
	
	// We want A to be about sqrt (2N) / M. Thus the primes we're looking at should be around [sqrt(2n)/M)] ^ (1/k).
	// However, there is no particular requirement that they be of comparable sizes, so this is just a rough guide.
	
	/* Now that we've chosen A, we can compute the values of B. The goal is to produce all values of B which satisfy 
	 * B^2 = N (mod A). Since A is composite, this is a little tricky, but we know it's prime factorization (that's how we
	 * constructed it). First, we find for each g_i the values r_i_1 and r_i_2 such that r_i^2 = N (mod g_i). Then, we set
	 * up a system of modular congruences 
	 *
	 * 	b ~= r_0_? (mod g_0)
	 * 	b ~= r_1_? (mod g_1)
	 * 	...
	 * 	b ~= r_k_? (mod g_k).
	 *
	 * which can be solved via the Chinese Remainder Theorem. The Wolfram Mathworld page on the CRT provides a helpful formula:
	 *
	 * 	b ~= r_0_? * j_0 * (A / g_0) + ... r_k_? * j_k * (A / g_k)  (mod A). We just take B = that expression.
	 *
	 * where the j_i are determined as 
	 * 	
	 * 	j_i * (A / g_i) ~= 1 (mod g_i), or in other words, j_i is the inverse of (A / g_i) (mod g_i).
	 *
	 *
	 * We repeat this procedure, 2^k times, taking every possible combination of which r_i_? we select 
	 * (? is independently replaced by either 1 or 2 in each congruence). 
	 * We disregard values of B > A/2, since these are the negations of other, smaller B. Hence, we end up with 2^(k-1) values for B.
	*/

	uint32_t r[k][2];
	uint32_t j[k];

	// first find r_i^2 ~= N (mod g_i)
	for (int i=0; i<k; i++){
		r[i][0] = find_root (ns->N, pg->gvals[i]);
		r[i][1] = pg->gvals[i] - r[i][0];
	}

	// now do the main loop
	mpz_t t1, t2;		// temps
	mpz_inits (t1, t2, NULL);
	for (int z = 0; z < (1 << k); z++){	// the ith bit of z will control which root r_i_? is used for the current B computation.
		mpz_set_ui (pg->bvals[z], 0);	// clear the B-value
		for (int i=0; i < k; i++){	// for each G-value
			mpz_set_ui (t2, pg->gvals[i]);			// t2 = g_i
			mpz_divexact_ui (t1, pg->a, pg->gvals[i]);	// t1 = A / g_i
			mpz_invert (t2, t1, t2);	// t2 = inverse of (A / g_i) (mod g_i) = j_i
			mpz_mul (t2, t2, t1);		// t2 = j_i * (A / g_i)
			// pick the root according to the i'th least siginificant bit of z.
			int root = (z & (1 << i)) == 0 ? 0 : 1;
			mpz_mul_ui (t2, t2, r[i][root]);
			mpz_add (pg->bvals[z], pg->bvals[z], t2);
		}
		mpz_mod (pg->bvals[z], pg->bvals[z], pg->a);	// take the sum mod A.
	}

	// Now that we've chosen A and determined the values of B, we compute A^-1 (mod p) for each prime in the factor base.
	mpz_t p, temp;
	mpz_init (p, temp, NULL);
	for (int i=0; i<ns->fb_len; i++){
		mpz_set_ui (p, ns->fb[i]);
		int invertable = mpz_invert (temp, pg->a, p);
		if (invertable){	// good
			pg -> ainverses[i] = mpz_get_ui (temp);
		} else {
			// hmmm. This is interesting. This should happen only when fb[i] == g_j. I don't know yet what exactly should happen here.
		}
	}
	mpz_clears (p, temp, NULL);
}

void generate_poly (poly_t *p, poly_group_t *pg, nsieve_t *ns, int i){
	p->a = pg->a;
	p->b = p->bvals[i];
	// compute C = (b^2 - n) / a
	mpz_mul (p->c, b, b);		 // C = b^2
	mpz_sub (p->c, p->c, ns->N);	 // C = b^2 - N
	mpz_divexact (p->c, p->c, p->a); // C = (b^2 - N) / a
}

void poly (mpz_t res, poly_t *p, uint32_t offset){
	// evaluate Ax^2 + 2Bx + C  where x = istart + offset
	// = (((A * X) + 2B) * X) + C
	mpz_t x;
	mpz_inits  (x, NULL);
	mpz_add_ui (x, p->istart, offset);	// x = istart + offset
	mpz_mul    (res, p->a, x);
	mpz_add    (res, res, p->b);
	mpz_add    (res, res, p->b);
	mpz_mul    (res, res, x);
	mpz_add    (res, res, p->c);
	mpz_clear  (x);
}

#include "poly.h"

/* Initialize a polynomial group structure */
void polygroup_init (poly_group_t *pg, nsieve_t *ns){
	mpz_init (pg->a);
	pg->ainverses = (uint32_t *) malloc(ns->fb_len * sizeof(uint32_t));
	pg->bvals = (mpz_t *) malloc( ns->bvals * sizeof (mpz_t));
	pg->relns = (rel_t **) calloc (PG_REL_STORAGE, sizeof (rel_t *));
	pg->nrels = 0;
	pg->victim = NULL;
	for (int i=0; i < ns->bvals; i++){
		mpz_init (pg->bvals[i]);
	}
}

void polygroup_free (poly_group_t *pg, nsieve_t *ns){
	mpz_clear (pg->a);
	free (pg->ainverses);
	free (pg->relns);
	for (int i=0; i < ns->bvals; i++){
		mpz_clear (pg->bvals[i]);
	}
	free (pg->bvals);
}

void poly_init (poly_t *p){
	mpz_inits (p->a, p->b, p->c, NULL);
}

void poly_free (poly_t *p){
	mpz_clears (p->a, p->b, p->c, NULL);
}

void poly_print (poly_t *p){
	printf("Q(x) = ");
	mpz_out_str(stdout, 10, p->a);
	printf("*x^2 + 2*");
	mpz_out_str(stdout, 10, p->b);
	printf("*x + ");
	mpz_out_str(stdout, 10, p->c);
	printf("\n");
}

// rough estimate of the prime-counting function. Implemented as x/ln(x). Used in selection of 'k'
uint32_t pi (uint32_t x){	
	return (uint32_t) (x / log(x));
}

/* Analogue of the mpz_nextprime function, but going backwards. */
void mpz_prevprime (mpz_t res){
	mpz_sub_ui(res, res, 1);
	while (!mpz_probab_prime_p(res, 10)){	// this will be fast for the even values, since 
						// mpz_probab_prime_p does a little trial division first.
		mpz_sub_ui(res, res, 1);
	}
}

/* This routine will pick a value of k and an appropriate set of g values. 
 * We want to pick the largest value of k that will give us enough possible polynomial groups. If k is 
 * too large, the g values become very small, and to produce enough distinct A, we end up with most A well 
 * outside the optimal range for our N and M, because the differences between the g_i are significant fractions
 * of their values. Consider an optimal A value of 10^20 and k = 10; each g is about 100; to get enough values
 * of g, we'd need to go very far away from 100.  We take the following approach:
 *
 * Pick a minimum number of polynomial groups, minP, say 25000. Pick a rather large k. Also pick a value 
 * 0 < c < 1. c will constrain the range of allowable g values. Let Aopt be the optimal value of A (that 
 * is, sqrt(2N)/M). Then we only take g between (c * Aopt)^(1/k) and (1/c * Aopt)^(1/k). 
 *
 * Determine the smallest q such that q C k > minP. (Note that these q can be precomputed (they don't depend 
 * on N) and stored as compile-time constants). The range must contain at least q valid values of g. Instead 
 * of actually computing the number of valid g, we can estimate it, using an approximation to pi(x), the 
 * prime-counting function. If we have a range, [gmin..gmax], then we can estimate the # of valid g as
 * 	
 * 	   pi(gmax) - pi(gmin)
 * #g ~=   ------------------            pi(x) ~= x / ln(x), for instance.
 * 	            2
 * 
 * where the 2 is present because only about half of the time will N be a quadratic residue mod p. 
 *
 * If #g is too small (i.e., less than q), decrement k and try again.
 *
 * Once a suitable k has been found, find q/2 values of g on either side of Aopt^(1/k). This is the gpool.
*/

/* The precomputed values of Q for minP = 25000. Thank you Mathematica for solving the equations. It is 
 * unlikely that a factorization will need more than 25000 _groups_ of polynomials, especially since for
 * larger factorizations, k will be high, so each group might contain 1024 or more polys. */
uint32_t q[] = {25000, 224, 54, 30, 22, 19, 18, 18, 18, 18, 18, 19};

/* Initialize the gpool, which involves selecting k and computing all of the allowable g values, as 
 * described above. */
void gpool_init (poly_gpool_t *gp, nsieve_t *ns){
	mpz_t aopt, temp, g;
	mpz_inits(aopt, temp, g, NULL);

	mpz_mul_ui    (aopt, ns->N, 2);
	mpz_sqrt      (aopt, aopt);
	mpz_tdiv_q_ui (aopt, aopt, ns->M * BLOCKSIZE / 2);	// yay, aopt = sqrt(2N)/M.

	int c_num = 6;
	int c_den = 10;	// c = 0.6

	int k = 12;
	while (k > 2) {	// we never really want k = 1.
		// compute gmin and gmax
		mpz_mul_ui(temp, aopt, c_num);
		mpz_tdiv_q_ui(temp, temp, c_den);
		mpz_root (temp, temp, k);
		uint32_t gmin = mpz_get_ui(temp);

		mpz_mul_ui(temp, aopt, c_den);
		mpz_tdiv_q_ui(temp, temp, c_num);
		mpz_root (temp, temp, k);
		uint32_t gmax = mpz_get_ui(temp);

		// approximate #g and compare to q[k-1]	(since q[0] has data for k=1)
		uint32_t approx_ng = (pi(gmax) - pi(gmin)) / 2;
		if (approx_ng >= q[k-1]){	// then this value of k is good.
			break;
		}
		k--;
	}
	if (k == 0){
		printf("Fatal error: no values of k were deemed viable.\n");
		exit(1);
	}
	mpz_root (temp, aopt, k);	// temp contains the 'central' value for our range.
	int ng = q[k-1];	// we must find this many g values.
	gp->gpool = (uint32_t *) malloc (ng * sizeof (uint32_t));	// allocate the g pool

	// fill in the upper half of the pool
	int pos = ng/2;
	mpz_set (g, temp);
	while (pos < ng){
		mpz_nextprime(g, g);
		if (mpz_kronecker (ns->N, g) == 1){
			gp->gpool[pos] = mpz_get_ui(g);
			pos++;
		}
	}
	// now do the lower half
	pos = ng/2 - 1;
	mpz_set (g, temp);
	while (pos >= 0){
		mpz_prevprime(g);
		if (mpz_kronecker (ns->N, g) == 1){
			gp->gpool[pos] = mpz_get_ui(g);
			pos--;
		}
	}
	// initialize the rest of the fields, including the frogs (see block comment in front of advance_gpool)
	gp->ng = ng;
	gp->frogs = (uint32_t *)(malloc(k * sizeof(uint32_t)));	// these hold indices into gpool.
	for (int i=0; i < k; i++){
		gp->frogs[i] = ng - (k - i - 1) - 1;	// frogs are stored in decreasing order
	}
	gp->k = k;
	ns->k = k;
	ns->bvals = 1 << (k-1);
}

/* We need a way to sequentially generate unique values of A. We initialize k 'frogs' to be the last 
 * (highest) k primes in the gpool. Then, whenever advance_gpool is called, the highest frog (in the k'th 
 * position of frogs) hops one slot backwards. If it can't (it's at the end of gpool), the next one back 
 * tries to hop backwards. This continues until a frog is found who can hop (Billy). Then, all the frogs 
 * higher than Billy (which were jammed up at the bottom of gpool) are set immediately following Billy. 
 * A large enough gpool will be chosen so that we won't run out of combinations. 
 *
 * This also provides a mechanism for distributing the sieving work to the different threads. We make
 * one gpool by calling gpool_init, then make a copy for each thread. Each copy will, as part of its 
 * initialization, be advanced as many times as its thread number. All of them now have different A values.
 * When a thread needs a new poly group (it's done sieving the previous one), it calls advance_gpool as 
 * many times as there are sieve threads. This way all of the threads are guaranteed to produce distinct
 * A values, without resorting to random selection and the potential for duplication of work.
*/
void advance_gpool (poly_gpool_t *gp, poly_group_t *group){	// advances the frogs, and sets the gvals field of 'group'
	int k = gp->k;
	// update group->gvals
	if (group != NULL){
		for (int i=0; i<k; i++){
			group->gvals[i] = gp->gpool[gp->frogs[i]];
		}
	}

	// now advance the frogs.
	int j = 0;
	while (gp->frogs[j] == j){
		j++;
		if (j == k){	// then we ran out of polynomials.
			printf("Fatal error: Out of polynomials!!!\n");
			exit(1);
		}
	}
	// j is now the index of the first frog that is not jammed against the end of the pool. 
	// Or edge of the pond, I suppose. Its name is Billy.
	// I guess for the frogs' sake its a good thing I'm not writing this in Python.
	gp->frogs[j] --;	// advance Billy.
	// now move all of the beached frogs back to being just ahead of billy, with no space between any of them.
	int billy = j;
	j--;
	while (j >= 0){	// most of the time, this won't even execute, since we'll just be advancing the highest frog.
		gp->frogs[j] = gp->frogs[billy] + (j - billy);
		j--;
	}
}

/* This will perform all of the work to set up the polygroup so that we may pull out polynomials from it. It
 * does some precomputation (of A^-1 (mod p)) as well. */
void generate_polygroup (poly_gpool_t *gp, poly_group_t *pg, nsieve_t *ns){
	/* This is a tricky one. First we must choose A, by picking k primes g_i, and multiplying them together.
	 * Then we need to find all of the values of B which satisfy  B^2 = N (mod a). There will be 2^(k-1) of them. 
	 * Then we will compute the values of A^-1 (mod p) for each p in the factor base. This is really a precomputation
	 * to speed up the computation of Q(x) = 0 (mod p). 
	*/
	
	for (int i=0; i < ns->nthreads; i++){	// advance the gpool nthreads times to get our next set of gvals
		advance_gpool (gp, pg);
	}
	mpz_set_ui (pg->a, pg->gvals[0]);	// multiply all of the g-values together to get A.
	for (int i=1; i < gp->k; i++){
		mpz_mul_ui(pg->a, pg->a, pg->gvals[i]);
	}
	
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
	int k = ns->k;

	uint32_t r[k][2];

	// first find r_i^2 ~= N (mod g_i)
	for (int i=0; i<k; i++){
		r[i][0] = find_root (ns->N, pg->gvals[i]);
		r[i][1] = pg->gvals[i] - r[i][0];
	}

	// now do the main loop
	mpz_t t1, t2;		// temps
	mpz_inits (t1, t2, NULL);
	int w = 0;
	for (int z = 0; z < (1 << k); z++){	// the ith bit of z will control which root r_i_? is used for the current B computation.

		if (w == 1 << (k-1)) break;	// then we've found enough, and we don't need to look at the rest. This is really mostly just
						// a safegard for us to not keep setting pg->bvals[w] when w is out of range.

		mpz_set_ui (pg->bvals[w], 0);	// clear the B-value
		for (int i=0; i < k; i++){	// for each G-value
			mpz_set_ui (t2, pg->gvals[i]);			// t2 = g_i
			mpz_divexact_ui (t1, pg->a, pg->gvals[i]);	// t1 = A / g_i
			mpz_invert (t2, t1, t2);	// t2 = inverse of (A / g_i) (mod g_i) = j_i
			mpz_mul (t2, t2, t1);		// t2 = j_i * (A / g_i)
			// pick the root according to the i'th least siginificant bit of z.
			int root = (z & (1 << i)) == 0 ? 0 : 1;
			mpz_mul_ui (t2, t2, r[i][root]);
			mpz_add (pg->bvals[w], pg->bvals[w], t2);
		}
		mpz_mod (pg->bvals[w], pg->bvals[w], pg->a);	// take the sum mod A.
		// if it's over A/2, reject it as the negation of another square root.
		mpz_mul_ui(pg->bvals[w], pg->bvals[w], 2);
		if (mpz_cmp (pg->bvals[w], pg->a) < 1){	// it's good
			mpz_divexact_ui (pg->bvals[w], pg->bvals[w], 2);
			w ++;
		}
	}

	// Now that we've chosen A and determined the values of B, we compute A^-1 (mod p) for each prime in the factor base.
	mpz_t p, temp;
	mpz_inits (p, temp, NULL);
	for (int i=0; i<ns->fb_len; i++){
		mpz_set_ui (p, ns->fb[i]);
		int invertable = mpz_invert (temp, pg->a, p);
		if (invertable){	// good
			pg -> ainverses[i] = mpz_get_ui (temp);
		} else {
			// hmmm. This is interesting. Does this ever occur? Maybe I should figure this out. 
			// It seems to work fine, though.
		}
	}
	mpz_clears (p, temp, NULL);
}

extern uint32_t get_offset (uint32_t, int, int, int, poly_t *, poly_group_t *, nsieve_t *);

/* Generate a polynomial from a group. generate_polygroup should have been called on the group before
 * this method is called. It will get polynomial #i, corresponding to the i'th computed b value */
void generate_poly (poly_t *p, poly_group_t *pg, nsieve_t *ns, int i){
	p->group = pg;

	mpz_set(p->a, pg->a);
	mpz_set (p->b, pg->bvals[i]);
	// compute C = (b^2 - n) / a
	mpz_mul (p->c, p->b, p->b);	 // C = b^2
	mpz_sub (p->c, p->c, ns->N);	 // C = b^2 - N
	mpz_divexact (p->c, p->c, p->a); // C = (b^2 - N) / a

	mpz_t temp;
	mpz_init (temp);

	p->M = ns->M;

	// compute -B % p for each prime in the factor base. This is another precomputation for get_offset.
	p -> bmodp = (uint32_t *) malloc (ns->fb_len * sizeof(uint32_t));
	for (int i=0; i < ns->fb_len; i++){
		mpz_neg (temp, p->b);
		p->bmodp[i] = mpz_mod_ui (temp, temp, ns->fb[i]);
	}
	mpz_clear(temp);
}

void poly (mpz_t res, poly_t *p, int32_t x){
	// evaluate Ax^2 + 2Bx + C 
	// = (((A * X) + 2B) * X) + C
	mpz_mul_si (res, p->a, x);	// res = ax
	mpz_add    (res, res, p->b);	// res = ax + b
	mpz_add    (res, res, p->b);	// res = ax + 2b
	mpz_mul_si (res, res, x);	// res = ax^2 + 2bx
	mpz_add    (res, res, p->c);	// res = ax^2 + 2bx + c
}

/* Self-check on the validity of a relation. Makes sure the LHS and RHS correspond. This
 * will detect errors in the factor lists, or perhaps in polynomial stuff as well. */
int rel_check (rel_t *rel, nsieve_t *ns){
	mpz_t facprod, pol, temp;
	mpz_inits (facprod, pol, temp, NULL);

	poly (pol, rel->poly, rel->x);
	poly (temp, rel->poly->group->victim->poly, rel->poly->group->victim->x);
	mpz_mul (pol, pol, temp);

	mpz_set_ui (facprod, rel->cofactor);
	fl_entry_t *entry = rel->factors;
	while (entry != NULL){
		if (entry->fac == 0){
			mpz_neg (facprod, facprod);
		} else {
			if (entry->fac < 0 || entry->fac > ns->fb_len){
				printf ("fac out of bounds error: fac = %d\n", entry->fac);
				mpz_clears (pol, temp, facprod, NULL);
				return 0;
			}
			mpz_mul_ui (facprod, facprod, ns->fb[entry->fac - 1]);
			if (!mpz_divisible_ui_p(pol, ns->fb[entry->fac - 1])){
				mpz_clears (facprod, pol, NULL);
				printf ("divisibility failed\n");
				return 0;
			}
		}
		entry = entry->next;
	}
//	mpz_mod (pol, pol, ns->N);
//	mpz_mod (facprod, facprod, ns->N);
	int res = (mpz_cmp (pol, facprod) == 0) ? 1 : 0;
	mpz_clears (facprod, pol, NULL);
	return res;
}

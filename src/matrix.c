#include "matrix.h"

/* The Gaussian Elimination code is modeled on a post on Programming Praxis. It proceeds basically as follows:
 *
 * Construct a square history matrix, with dimension equal to the # of rows in the exponent matrix. Initialize it
 * to the identity matrix. This records which exponent vectors have been combined together. Now:
 *
 * for each column in the matrix, working right to left:
 * 	find the 'pivot' vector closest to the top with rightmost 1 in the current column. If none exists, go to the next column.
 * 		for each other vector (Yolanda) with rightmost 1 in the current column:
 * 			Yolanda = Yolanda XOR pivot	(for both the exponent and history matricies)
 *
 * Any zero row-vectors in the exponent matrix represent dependencies that can be converted into a congruence of squares. 
*/

void solve_matrix (nsieve_t *ns){
	printf ("\nStarting gaussian elimination... \n");
	long start = clock();

	const int hmlen = (ns->rels_needed -1)/64 + 1;
	const int hmsize = ns->rels_needed;
	uint64_t *history[hmsize];
	uint32_t rmos[hmsize];	// for keeping track of the position of the rightmost 1 in each exponent vector.
				// caching these values greatly accelerates the matrix solving, since they are
				// used often but change rarely.
	for (int i=0; i < hmsize; i++){
		history[i] = (uint64_t *) calloc (hmlen, sizeof(uint64_t));
		flip_bit (history[i], i);
		rmos[i] = rightmost_1 (ns->relns[i].row, ns->fb_len);
	}
#define MAT_CHECK	// enables a couple of self-checks on the matrix solving. However, it requires making a 
#ifdef MAT_CHECK	// copy of the exponent side of the matrix, so should be turned off for 'release' versions.
	uint64_t *expm[ns->rels_needed];
	for (int i=0; i < hmsize; i++){
		expm[i] = (uint64_t *) malloc (ns->row_len * sizeof(uint64_t));
		for (int j=0; j<ns->row_len; j++){
			expm[i][j] = ns->relns[i].row[j];
		}
//		print_row(expm[i], 64 * ns->row_len);
	}
#endif
	const int expm_rows = ns->rels_needed;
	const int expm_cols = ns->fb_len + 1;

	// this is the main loop of the Gaussian Elimination
	// col starts at fb_len because there are fb_len + 1 columns in the matrix.
	for (int col = expm_cols-1; col >= 0; col --){	
		if ((expm_cols - col) % 50 == 0){	// print progress report
			printf("Column %d of %d\r", expm_cols - col, expm_cols);
			fflush(stdout);
		}
		int pivot = 0;
		while (pivot < expm_rows && rmos[pivot] != col){
			pivot ++;
		}
		for (int yolanda = pivot + 1; yolanda < expm_rows; yolanda++){
			if (rmos[yolanda] == col){
				xor_row (ns->relns[yolanda].row, ns->relns[pivot].row, ns->row_len);
				xor_row (history[yolanda], history[pivot], hmlen);

				/* update the cached value of the rightmost 1. We know the rightmost 1 of the 
				 * xor'ed result has to be to the left of 'col', since that was the position
				 * of the rightmost 1 of both yolanda and pivot, and the xor sets that to 0. */
				rmos[yolanda] = rightmost_1 (ns->relns[yolanda].row, col);	
			}
		}
	}
	printf("\nMatrix solved; deducing factors...\n");
	ns->timing.matsolve_time = clock() - start;
	start = clock();

/* Factor deduction and dealing with multiple polynomials, etc.
 *
 * First, note that we are trying to produce congruences of the form X^2 ~= Y^2 (mod N). gcd (X-Y, N) will (about half
 * of the time) yield a nontrivial factor of N. Each polynomial Q(x) is constructed so that 
 *
 * 	a Q(x) = (ax + b)^2 - N   ~= H^2 (mod N), for H = (ax + b) (mod n) 
 *
 * The sieve will find many values of x such that Q(x) = y is smooth over the factor base. Note that a*Q(x) may or may
 * not actually factor over the factor base (this depends on whether the g_i are in the factor base or not, and we cannot
 * assume that they will be), but we really need the LHS to be a square mod N, so we must keep this factor a. This creates
 * a puzzle, until we realize that we can select one relation, say   a * Q(x_0) = y_0   and multiply all of the other 
 * relations that share the same 'a' value by it. We then have a list of relations like this:
 *
 * 	a^2 Q(x_1) Q(x_0) = y_1 * y_0
 * 	...
 * 	a^2 Q(x_n) Q(x_0) = y_n * y_0
 *
 * Now all of the LHS's are squares mod N, and the RHS's are still smooth over the factor base. Define
 *
 * 	H_q,i = [(ax_i + b) * (ax_0 + b)] % N	 for the a,b values of the polynomial q. 
 *
 * Then we may rewrite the relations above like this:
 *
 * 	(H_Q,i)^2 ~= y_i * y_0   (mod N).
 *
 * The linear algebra step will then select a subset of these relations so that the RHS is a square (in the integers). 
 * Note that the relations may come from different polynomials. We get
 *
 * 	(H_p,i)^2 * (H_q,j)^2 * ... ~= y_p,i * y_q,j ...  (mod N)
 * 	   "           "         "  ~= v^2		  (mod N)
 *
 * This is our desired congruence of squres. 
*/

/* It is impractical to compute the values of lhs^2 and rhs^2 explicitly, then take their square roots and reduce mod N;
 * since several hundred relations may be combined into one congruence, these numbers might end up being many megabytes
 * long (ask me how I know), which can slow things down a lot. It is fairly clear how to compute the lhs without
 * computing lhs^2; since the coefficients of the polynomials are all known, it is just a matter of multiplying by
 * H_q,i instead of (H_q,i)^2; furthermore, the value can be reduced mod N at each step.
 *
 * For the RHS, each relation's piece is not a square, so a similar approach is not possible. Instead, we keep a table,
 * with one entry for each prime in the FB (and -1), and for each relation, we walk down its factor list, adding 1
 * to the corresponding entries of the table. Doing this for all the relations will yield a table with all entries 
 * being even (if all has gone well); the entries represent the exponents on the primes that when multiplied together
 * produce the value of rhs^2. Looping over the values, and multiplying rhs by fb[i-1]^(table[i]/2) (and reducing mod
 * N - this is very important that we can do this now!) will result in the correct computation of rhs.
*/
	mpz_divexact_ui (ns->N, ns->N, ns->multiplier);	// otherwise we will uncover the multiplier as a factor.
	mpz_t ncopy, temp;
	mpz_init_set (ncopy, ns->N);
	mpz_inits (temp, NULL);
	uint16_t factor_counts[ns->fb_len+1];	// this is the table we mentioned in the above comment.
	for (int row = 0; row < expm_rows; row ++){
		if (is_zero_vec (ns->relns[row].row, ns->row_len)){	// we found a dependency
			memset (factor_counts, 0, 2 * (ns->fb_len + 1));	// clear the factor_counts table
			int relct = 0;
			int partialct = 0;
#ifdef MAT_CHECK
			// This code will verify that the matrix solving worked; that is, it will xor together all of the rows
			// specified in the history matrix, and verify that it is indeed the zero vector.
			uint64_t check [ns->row_len];
			clear_row (check, ns);
			for (int i=0; i < hmsize; i++){
				if (get_bit (history[row], i) == 1){
					xor_row (check, expm[i], ns->row_len);
				}
			}
			if (is_zero_vec (check, ns->row_len)){
//				printf("Check succeeded.\n");
			} else {
				printf("Check FAILED for row = %d\n", row);
			}
#endif
			// yay! we have a dependency. Now the ugly math begins.
			mpz_t lhs, rhs;	// we will end up with lhs^2 ~= rhs^2 (mod N)
					// the left hand side is the H_p,i and the right side is the y_p,i.
			mpz_inits (lhs, rhs, NULL);
			mpz_set_ui(lhs, 1);
			mpz_set_ui(rhs, 1);
			for (int relnum = 0; relnum < hmsize; relnum ++){
				if (get_bit (history[row], relnum) == 1){	// the relation numbered 'relnum' is included in the dependency
					matrel_t *m = &ns->relns[relnum];

					if (!rel_check (m->r1, ns)){		// one can never have too much checking.
						printf ("relation failed check. [%s]\n", m->r2==NULL?"full":"partial, r1");
					} else {
					//	printf ("relation passed check.\n");
					}
					multiply_in_lhs (lhs, m->r1, ns);
					add_factors_to_table (factor_counts, m->r1);
					relct ++;
					if (m->r2 != NULL){	// partial
						partialct ++;
						if (m->r1->cofactor != m->r2->cofactor){
							printf("AAAH - cofactors disagree! (%d and %d)\n", m->r1->cofactor, m->r2->cofactor);
						}
						if (!rel_check (m->r2, ns)){
							printf ("relation failed check. [partial, r2]\n");
						}
						multiply_in_lhs (lhs,  m->r2, ns);
						add_factors_to_table (factor_counts, m->r2);

						mpz_mul_ui (rhs, rhs, m->r1->cofactor);	// the cofactors aren't stored
									// in the lists, so we have to do them separately.
					}
				}
			}
			int is_good = construct_rhs (factor_counts, rhs, ns);
//			printf ("multiplied together %d relations, %d of which were from partials\n", relct, partialct);
			if ( ! is_good ) {	// more self-checks.
				printf ("construct_rhs check failed.\n");
				continue;
			}
			mpz_mod (lhs, lhs, ns->N);
			mpz_mod (rhs, rhs, ns->N);

			mpz_sub (temp, rhs, lhs);
			mpz_gcd (temp, temp, ncopy);	// take the gcd with ncopy instead of N, to avoid reprinting already found factors.
			mpz_abs (temp, temp);	// probably unneceesary
			/* Now we check to see if the factor we found was nontrivial (1 or n) */
			if (mpz_cmp_ui (temp, 1) > 0){
				if (mpz_cmp (temp, ns->N) != 0){	// then it's a nontrivial factor!!!
					if (mpz_divisible_p(ncopy, temp)){	// then we haven't found it before.
						if (mpz_probab_prime_p (temp, 10)){	// verify its primality
							mpz_out_str (stdout, 10, temp);
							printf (" (prp)\n");
							mpz_divexact(ncopy, ncopy, temp);
							/* If the cofactor is prime, print it out too */
							if (mpz_probab_prime_p (ncopy, 10)){
								mpz_out_str(stdout, 10, ncopy);
								printf (" (prp)\n");
								mpz_set_ui(ncopy, 1);
							}
							if (mpz_cmp_ui(ncopy, 1) == 0){	// we're done!
								mpz_clears(lhs, rhs, temp, ncopy, NULL);
								ns->timing.facdeduct_time = clock() - start;
								return;
							}
						}
					}
				}
			}
			mpz_clears(lhs, rhs, NULL);
		}
	}

	if (mpz_cmp_ui(ncopy, 1) != 0){
		mpz_out_str (stdout, 10, ncopy);
		if (mpz_probab_prime_p (ncopy, 10)){
			printf (" (prp)\n");
		} else {
			/* It is a sad day. Most likely there is a bug. */
			printf (" (c)\n");
		}
	}
	ns->timing.facdeduct_time = clock() - start; 
}

/* For each entry of rel, increment position rel->fac of the table (the factor lists are storing
 * matrix row positions, not the actual primes) */
void add_factors_to_table (uint16_t *table, rel_t *rel){
	fl_entry_t *factor = rel->factors;
	while (factor != NULL){
		table[factor->fac] ++;
		factor = factor->next;
	}
}

/* Given a filled out table, construct the right hand side of our fancy congruence */
int construct_rhs (uint16_t *table, mpz_t rhs, nsieve_t *ns){	// returns nonzero on success, zero on failure.
	if (table[0] % 2 != 0){
		printf("Error: table[%d] is not even (=%d)\n", 0, table[0]);
		return 0;
	}
	if ((table[0]/2) % 2 == 1){	// then we need to negate rhs
		mpz_neg (rhs, rhs);
	}
	mpz_t temp;
	mpz_init (temp);
	for (int i=1; i < ns->fb_len + 1; i++){
		if (table[i] % 2 != 0){
			printf("Error: table[%d] is not even (=%d)\n", i, table[i]);
			mpz_clear (temp);
			return 0;
		}
		if (table[i] > 0){
			mpz_ui_pow_ui (temp, ns->fb[i-1], table[i]/2);
			mpz_mul (rhs, rhs, temp);
			mpz_mod (rhs, rhs, ns->N);
		}
	}
	mpz_clear(temp);
	return 1;
}

/* Take care of what needs to be done to the LHS for this relation. This is called once for each
 * component of the partial. */
void multiply_in_lhs (mpz_t lhs, rel_t *rel, nsieve_t *ns) {	// this should work just as well for partials as for fulls.
	mpz_t temp;
	mpz_init (temp);

	rel_t *victim = rel->poly->group->victim;

	// multiply the (Ax_victim + B_victim) for the victim that was used for this reln's poly group.
	mpz_set_si (temp, victim->x);
	mpz_mul (temp, temp, victim->poly->a);
	mpz_add (temp, temp, victim->poly->b);
	mpz_mul (lhs, lhs, temp);

	// multiply the (Ax_i + B_i) for our relation
	mpz_set_si (temp, rel->x);
	mpz_mul (temp, temp, rel->poly->a);
	mpz_add (temp, temp, rel->poly->b);
	mpz_mul (lhs, lhs, temp);

	// It's easier to deal with the A^2 factor here in the LHS, since we're actually looping over
	// the relations here already. We multiply by the modular multiplicative inverse A^-1 (mod N)
	// on the left.
	mpz_invert (temp, victim->poly->a, ns->N);
	mpz_mul (lhs, lhs, temp);

	mpz_mod (lhs, lhs, ns->N);	// reduce mod N to keep things small.

	mpz_clear (temp);
}

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
	for (int i=0; i < hmsize; i++){
		history[i] = (uint64_t *) calloc (hmlen, sizeof(uint64_t));
		flip_bit (history[i], i);
		rmos[i] = rightmost_1 (ns->relns[i].row, ns->fb_len);
	}
//#define MAT_CHECK
#ifdef MAT_CHECK
	uint64_t *expm[ns->rels_needed];
	for (int i=0; i < hmsize; i++){
		expm[i] = (uint64_t *) malloc (ns->row_len * sizeof(uint64_t));
		for (int j=0; j<ns->row_len; j++){
			expm[i][j] = ns->relns[i].row[j];
		}
		print_row(expm[i], 64 * ns->row_len);
	}
#endif
	const int expm_rows = ns->rels_needed;
	const int expm_cols = ns->fb_len + 1;
	for (int col = expm_cols-1; col >= 0; col --){	// col starts at fb_len because there are fb_len + 1 columns in the matrix.
		if ((expm_cols - col) % 50 == 0){
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
				rmos[yolanda] = rightmost_1 (ns->relns[yolanda].row, expm_cols-1);	// update the cached value of the rightmost 1.
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

	mpz_t ncopy, temp;
	mpz_init_set (ncopy, ns->N);
	mpz_inits (temp, NULL);
	for (int row = 0; row < expm_rows; row ++){
		if (is_zero_vec (ns->relns[row].row, ns->row_len)){
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
				printf("Check succeeded.\n");
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

					multiply_in_relation (lhs, rhs, m->r1, ns, 0);
					if (m->r2 != NULL){	// partial
						multiply_in_relation (lhs, rhs, m->r2, ns, 1);
					//	mpz_divexact_ui (rhs, rhs, m->r1->cofactor * m->r1->cofactor);	// divide by the cofactor
						mpz_set_ui (temp, m->r1->cofactor);
						mpz_invert (temp, temp, ns->N);
						mpz_mul (temp, temp, temp);
						mpz_mul (lhs, lhs, temp);
						mpz_mod (lhs, lhs, ns->N);
					}
					// is this all there is to it?

				}
			}
			// ok, now we have lhs^2 ~= rhs (mod N). rhs should be a square if all went well, so we shall take its sqrt to get
			// 	lhs^2 ~= v^2 (mod N)
			if (mpz_cmp_ui(rhs, 0) < 0){	// this shouldn't happen, since there should be an even number of (-1) factors.
				mpz_neg(rhs, rhs);
				printf("WARNING - found negative rhs; this should not happen; proceeding with taking |rhs|.\n");
			} else {
//				printf ("GOOD - rhs was positive.\n");
			}
			mpz_sqrt (rhs, rhs);

			mpz_mod (rhs, rhs, ns->N);	// this is only necessary to print out the congruence of squares.
		//	fprintf("Found congruence of squares: ");
		/*
			mpz_out_str(stderr, 10, lhs);
			fprintf(stderr, " ");
			mpz_out_str(stderr, 10, rhs);
			fprintf(stderr, "\n");
		*/
			mpz_sub (temp, rhs, lhs);
			mpz_gcd (temp, temp, ncopy);	// take the gcd with ncopy instead of N, to avoid reprinting already found factors.
			mpz_abs (temp, temp);	// probably unneceesary
			if (mpz_cmp_ui (temp, 1) > 0){
				if (mpz_cmp (temp, ns->N) != 0){	// then it's a nontrivial factor!!!
					if (mpz_divisible_p(ncopy, temp)){	// then we haven't found it before.
						if (mpz_probab_prime_p (temp, 10)){
							mpz_out_str (stdout, 10, temp);
							printf (" (prp)\n");
							mpz_divexact(ncopy, ncopy, temp);
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
			printf (" (c)\n");
		}
	}
	ns->timing.facdeduct_time = clock() - start; }

void multiply_in_relation (mpz_t lhs, mpz_t rhs, rel_t *rel, nsieve_t *ns, int partial){	// this should work just as well for partials as for fulls.
	mpz_t temp;
	mpz_init (temp);

	rel_t *victim = rel->poly->group->victim;

	// lhs *= victimH * relH;  if victim poly = p and rel poly = q, this is (p.a * victim.x + p.b) * (q.a * rel.x + q.b)
	mpz_set_si (temp, victim->x);
	mpz_mul (temp, temp, victim->poly->a);
	mpz_add (temp, temp, victim->poly->b);
	mpz_mul (lhs, lhs, temp);

	mpz_set_si (temp, rel->x);
	mpz_mul (temp, temp, rel->poly->a);
	mpz_add (temp, temp, rel->poly->b);
	mpz_mul (lhs, lhs, temp);

	mpz_invert (temp, victim->poly->a, ns->N);
	mpz_mul (lhs, lhs, temp);

	mpz_mod (lhs, lhs, ns->N);

	// now update RHS
	// rhs *= p(victim.x) * q(rel.x)
	poly (temp, victim->poly, victim->x);
	mpz_mul (rhs, rhs, temp);

	poly (temp, rel->poly, rel->x);
	mpz_mul (rhs, rhs, temp);
	// note that we can't reduce rhs mod N, since we have to take the sqrt of its final value in the integers.
}

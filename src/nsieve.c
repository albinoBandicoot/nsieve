#include "nsieve.h"

void era_sieve (nsieve_t *ns, char *vals){
	// assumes vals has been allocated and cleared with sufficient space for fb_bound chars.
	// 0 is prime, 1 is composite.
	for (int skip=2; skip < (int) (sqrt(ns->fb_bound) + 1); skip++){
		if (vals[skip-2] == 1){
			continue;
		}
		for (int pos=2*skip; pos < ns->fb_bound; pos += skip){
			vals[pos-2] = 1;
		}
	}
}

void extract (nsieve_t *ns, char *vals){
	// first count the number of primes found such that (n/p) = 1. (n/p) is the Legendre symbol.
	int count = 0;
	for (int i=0; i < ns->fb_bound; i++){
		if (vals[i] == 0){
			if (mpz_kronecker_ui(ns->N, i+2) == 1){
				count++;
			}
		}
	}
	ns->fb_len = count;
	ns->rels_needed = ns->fb_len + ns->extra_rels;
	ns->fb = (uint32_t *)(malloc(count * sizeof(uint32_t)));
	int w = 0;
	for (int i=0; i < ns->fb_bound; i++){
		if (vals[i] == 0 && mpz_kronecker_ui(ns->N, i+2) == 1){
			ns->fb[w] = i+2;
			w++;
		}
	}
}

/* Generates the factor base, given the bound fb_bound; will also compute the roots of N mod each p. */
void generate_fb (nsieve_t *ns){
	// first we use the sieve of Eratosthenes and the mpz_kronecker function to find a list of primes p < fb_bound such that (n/p) = 1.
	char *erasieve_vals = (char *)(calloc(ns->fb_bound, sizeof(char)));
	era_sieve (ns, erasieve_vals);
	extract (ns, erasieve_vals);
	free (erasieve_vals);

	// now we compute the modular square roots of n mod each p. 
	ns->roots = (uint32_t *)(malloc(ns->fb_len * sizeof(uint32_t)));
	for (int i=0; i<ns->fb_len; i++){
		ns->roots[i] = find_root (ns->N, ns->fb[i]);	// see common.c for the implementation of this method.
		// note that there are actually 2 square roots; however, the second may be obtained readily as p - sqrt#1, so only one is stored.
		ns->fb_logs[i] = fast_log (ns->fb[i]);
	}
}

/* Initialization and selection of the parameters for the factorization. This will fill allocate space for and initialize all of the 
 * fields in the nsieve_t. Notably, it will generate the factor base, compute the square roots, and allocate space for the matrix and partials. 
 */
void nsieve_init (nsieve_t *ns, mpz_t n){
	// all of the parameters here are complete BS for now.
	mpz_init_set (ns->N, n);
	ns->k = 3;
	ns->bvals = 1 << (ns->k - 1);
	ns->M = 8 * BLOCKSIZE;	
	ns->fb_bound = 2000;
	ns->extra_rels = 48;

	generate_fb (ns);

	ns->relns = (matrel_t *)(malloc((ns->fb_len + ns->extra_rels) * sizeof(matrel_t)));
	ht_init (ns);
}

/* Once ns has been initialized (by calling nsieve_init), this method is called to actaully perform the bulk of the factorization */
void factor (nsieve_t *ns){
	poly_gpool_t gpool;
	gpool_init (&gpool, ns);

	poly_group_t curr_polygroup;
	poly_t curr_poly;
	polygroup_init (&curr_polygroup, ns);
	poly_init (&curr_poly);

	block_data_t sievedata;
//	while (ns->nfull + ns->npartial < ns->rels_needed){
		// while we don't have enough relations, sieve another poly group.
		generate_polygroup (&gpool, &curr_polygroup, ns);
		printf("\nThe inverses of A (mod p) for each prime in the factor base: \n");
		for (int i=0; i<ns->fb_len; i++){
			printf("p = %d \tA^-1 = %d \tsqrt(n) mod p: %d\n", ns->fb[i], curr_polygroup.ainverses[i], ns->roots[i]);
		}
		for (int i = 0; i < ns -> bvals; i ++){
			generate_poly (&curr_poly, &curr_polygroup, ns, i);
			poly_print (&curr_poly);
			printf("\n");
		//	sieve_poly (&sievedata, &curr_poly, ns);
		}
		return;
//	}
	// now we have enough relations, so we build the matrix (combining the partials).
	build_matrix (ns);	// this includes combining the partials (and cycle-finding if we do double large primes)

	// Filter the matrix to reduce its size without reducing its yummy content. This will accelerate the matrix solving step, and also reduce the memory usage.
	filter (ns);
	
	// now we solve the matricies. The rest of the guts are in matrix.c, in the function solve_matrix.
	solve_matrix (ns);
}

int main (int argc, const char *argv[]){
	nsieve_t ns;
	mpz_t n;
	mpz_init (n);
	if (argc >= 2){
		mpz_set_str (n, argv[1], 10);
	} else {
		mpz_inp_str (n, stdin, 10);
	}
	nsieve_init (&ns, n);

	factor (&ns);
}

#include "nsieve.h"

/* This file contains the routines that coordinate the pieces defined in all of the other files. It
 * also does some initialization work */


/* First some, routines for generating the factor base */

/* The familiar sieve of Eratosthenes */
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

/* Extract the values of the primes from the sieve. Only extract those primes for which N is a
 * quadratic residue (mod p). This test was conveniently provided by the gmp routine 
 * mpz_kronecker_ui. */
void extract (nsieve_t *ns, char *vals){
	// first count the number of primes found such that (n/p) = 1. (n/p) is the Legendre symbol.
	int count = 0;
	for (int i=0; i < ns->fb_bound-2; i++){
		if (vals[i] == 0){
			if (i == 0 || i+2 == ns->multiplier || mpz_kronecker_ui(ns->N, i+2) == 1){	// we must admit 2, since N is always a QR mod 2, but there's something wierd about the kronecker symbol for 2  because it's not odd.
				count++;
			}
		}
	}
	// now we know how much space to allocate for our list
	ns->fb_len = count;
	ns->rels_needed = ns->fb_len + ns->extra_rels;
	ns->fb = (uint32_t *)(malloc(count * sizeof(uint32_t)));
	int w = 0;
	for (int i=0; i < ns->fb_bound-2; i++){
		if (vals[i] == 0 && (i == 0 || i+2 == ns->multiplier || mpz_kronecker_ui(ns->N, i+2) == 1)){
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

	// now we compute the modular square roots of n mod each p, and approximations to log_2(p).
	ns->roots = (uint32_t *)(malloc(ns->fb_len * sizeof(uint32_t)));
	ns->fb_logs = (uint8_t *)(malloc(ns->fb_len * sizeof(uint8_t)));
	for (int i=0; i<ns->fb_len; i++){
		if (ns->fb[i] == ns->multiplier){
			ns->roots[i] = 0;
		} else {
			ns->roots[i] = find_root (ns->N, ns->fb[i]);	// see common.c for the implementation of this method.
			if (ns->roots[i] > ns->fb[i]/2){
				ns->roots[i] = ns->fb[i] - ns->roots[i];	// normalize these to be the smaller of the two roots.
			}
			// note that there are actually 2 square roots; however, the second may be 
			// obtained readily as p - sqrt#1, so only one is stored.
		}
		ns->fb_logs[i] = fast_log (ns->fb[i]);
	}
}

/* Now some things for automatic selection of parameters. By 'automatic' this is more of a reflection
 * of my twiddling the parameters for various sizes of N and recording what worked best rather than
 * anything more mathematically motivated. 
*/
const int PARAM_FBBOUND = 1;
const int PARAM_LPBOUND = 2;
const int PARAM_M = 3;
const int PARAM_T = 4;

#define NPLEVELS  10	// number of entries in the params table
#define  NPARAMS  5
//						bits   FBB    LPB  M   T
const double params[NPLEVELS][NPARAMS] =   { 	{ 80,  1600,  50,  1 , 1.4},
						{100,  5000,  70,  1 , 1.45},
						{120,  8000,  90,  1 , 1.5},
						{140, 18000,  120, 1 , 1.5},
					    	{160, 36000,  120, 1 , 1.45},
						{180, 66000,  120, 1 , 1.45},
						{200, 120000, 150, 2 , 1.5}, 
						{220, 200000, 180, 2 , 1.55},/* beyond this point these are guesses. */ 	
						{230, 280000, 195, 2 , 1.55},
						{240, 360000, 210, 2 , 1.57}
					   };

/* Linearly interpolate parameters that were not manually overriden by the user between the adjacent
 * values in the parameter list. This should be called by select_parameters. 
*/
void set_params (nsieve_t *ns, int p1, int p2, double fac){
	// -1 indicates that the property was not manually overridden by the user via a command line argument.
	if (ns->fb_bound == -1) ns -> fb_bound = (uint32_t) (params[p1][PARAM_FBBOUND] * fac + params[p2][PARAM_FBBOUND] * (1 - fac));
	if (ns->lp_bound == -1){
		ns -> lp_bound = ns->fb_bound * (uint32_t) (params[p1][PARAM_LPBOUND] * fac + params[p2][PARAM_LPBOUND] * (1 - fac));
	} else {
		ns -> lp_bound *= ns->fb_bound;
	}
	if (ns->M == -1) ns -> M        = (uint32_t) (params[p1][PARAM_M] * fac + params[p2][PARAM_M] * (1 - fac));
	if (ns->T == -1) ns -> T        = (float)    (params[p1][PARAM_T] * fac + params[p2][PARAM_T] * (1 - fac));
	printf("Selected parameters: \n\tfb_bound = %d \n\tlp_bound = %d \n\tM = %d\n\tT - %f\n", ns->fb_bound, ns->lp_bound, ns->M, ns->T);
}

/* Perform automatic parameter selection. Only parameters not specified by the user will be chosen automatically */
void select_parameters (nsieve_t *ns){
	/* All of the choices are dependent solely on the number of bits in N */
	int bits = mpz_sizeinbase (ns->N, 2);
	printf("Choosing parameters for %d bit number... \n", bits);
	if (bits <= params[0][0]){	// smaller than the bottom of the table
		set_params(ns, 0, 0, 0);
	} else if (bits >= params[NPLEVELS-1][0]){	// above the end of the table
							// this will probably result in a very bad choice of 
							// parameters if you're much beyond the end.
		set_params(ns, NPLEVELS - 1, NPLEVELS - 1, 0);
	} else {
		int i = 0;
		while (i < NPLEVELS && params[i][0] < bits){
			i++;
		}
		set_params (ns, i, i-1, (bits - params[i-1][0]) / (params[i][0] - params[i-1][0]));
	}
}

const int nsmall_primes = 18;
uint32_t small_primes[] = {2,3,5,7,11,13,17,19,23,29,31,37,41,43,47,53,59,61};

/* In many instances it is actually better to try to factor tN for some small squarefree (usually prime)
 * t. The reasoning behind this is that if one can find a value of t for which tN is a quadratic residue
 * mod many small primes, the chances of a particular sieve value being smooth go up, since there are 
 * more small values in the factor base. This effect can be surprisingly large; the difference between
 * selecting a particularly bad multiplier and a particularly good one can amount to a factor of almost
 * 3 in the sieving time. It is therefore worthwhile to select wisely.
*/

/* This function rates a multiplier, giving it a score. The multiplier with the highest score is selected */
double get_multiplier_score (uint32_t mult, nsieve_t *ns){
	/* The modifier score is based on the 'Knuth-Schroeppel function,' defined as:
	 *
	 * f(t, N) = SUM_i g(p_i, tN) * log(p) - 0.5 log(t)	for all p_i in the factor base.
	 * 
	 * where	g(p, tN) = 2/p   if p does not divide t
	 * 			 = 1/p	 if p divides t
	 * 		g(2, tN) = 2	 if N ~= 1 (mod 8)
	 * 			 = 0	 otherwise
	 *
	 * The components of this function make a lot of intuitive sense when one considers it as an 
	 * approximation to the average amount the prime 'p' will contribute to a sieve bucket.
	 *
	 * In practice, we only sum over the smallest elements of the factor base (the ones in small_primes),
	 * since the 1/p factors quickly diminsh towards 0 as p increases. 
	*/
	mpz_t temp;
	mpz_init (temp);
	mpz_mul_ui (ns->N, ns->N, mult);
	double res = -0.5 * log(mult);
	if (mpz_mod_ui (temp, ns->N, 8) == 1){
		res += 2 * log(2.0);
	}
	mpz_clear (temp);
	for (int i=1; i < nsmall_primes; i++){
		if (small_primes[i] == mult){
			res += (1.0/small_primes[i]) * log(small_primes[i]);
		} else {
			if (mpz_kronecker_ui(ns->N, small_primes[i]) == 1){
				res += (2.0/small_primes[i]) * log(small_primes[i]);
			}
		}
	}
	mpz_divexact_ui (ns->N, ns->N, mult);
	return res;
}

/* Try multipliers up to about 60; pick the one with the highest score */
void select_multiplier (nsieve_t *ns){
	int idx = -1;
	double best = get_multiplier_score (1, ns);
//	printf("\tmultiplier %d has score %f\n", 1, best);
	for (int i=0; i < nsmall_primes; i++){
		double score = get_multiplier_score (small_primes[i], ns);
//		printf("\tmultiplier %d has score %f\n", small_primes[i], score);
		if (score > best){
			best = score;
			idx = i;
		}
	}
	if (idx == -1){
		ns->multiplier = 1;
	} else {
		ns->multiplier = small_primes[idx];
		mpz_mul_ui (ns->N, ns->N, small_primes[idx]);
	}
	printf("Selected multiplier %d.\n", ns->multiplier);
}
	
/* Initialization and selection of the parameters for the factorization. This will fill allocate space 
 * for and initialize all of the fields in the nsieve_t. Notably, it will generate the factor base, 
 * compute the square roots, and allocate space for the matrix and partials. 
*/
void nsieve_init (nsieve_t *ns, mpz_t n){
	long start = clock();
	mpz_init_set (ns->N, n);

	select_parameters (ns);
	if (ns->multiplier == -1){	// not preselected
		select_multiplier (ns);
	} else {
		mpz_mul_ui (ns->N, ns->N, ns->multiplier);
	}

	ns->nthreads = 1;
	pthread_mutex_init (&ns->lock, NULL);
	ns->info_npoly = 0;
	ns->info_npg = 0;

	ns->nfull = 0;
	ns->npartial = 0;
	ns->tdiv_ct = 0;
	ns->sieve_locs = 0;
	ns->extra_rels = 120;

	generate_fb (ns);

	/* Allocate space for the matrix relations; note that the actual rows of the matrix containing
	 * the packed bits are not allocated until the matrix building phase. */
	ns->relns = (matrel_t *)(malloc((ns->fb_len + ns->extra_rels) * sizeof(matrel_t)));
	ns->row_len = (ns->fb_len)/(8*sizeof(uint64_t)) + 1;	// we would need that to be ns->fb_len - 1, except we need to throw in the factor -1 into the FB. 
	if (ns->row_len % 2 == 1) ns->row_len ++;	// to take advantage of SSE instructions, we want to chunk by 128 bits.

	printf("There are %d primes in the factor base, so we will search for %d relations. The matrix rows will have %d 8-byte chunks in them.\n", ns->fb_len, ns->rels_needed, ns->row_len);

	ht_init (ns);
	ns->timing.init_time = clock() - start;
}

/* Run the SIQS with nthreads sieving threads. All other phases are single-threaded. Must have called 
 * nsieve_init prior to calling this, so that everything is set up. */
void multithreaded_factor (nsieve_t *ns, int nthreads){
	long start = clock ();
	
	/* Set up the threads; most of this code is in getting the various copies of the gpool set up. */
	ns->nthreads = nthreads;
	ns->threads = (pthread_t *) malloc(nthreads * sizeof (pthread_t));

	thread_data_t *td = (thread_data_t *) malloc (nthreads * sizeof (thread_data_t));

	poly_gpool_t gpool;
	gpool_init (&gpool, ns);
	printf("Using k = %d; gvals range from %d to %d.\n", ns->k, gpool.gpool[0], gpool.gpool[gpool.ng-1]);

	for (int i=0; i<nthreads; i++){
		td[i].ns = ns;
		td[i].gpool = gpool;	// this does a block copy. However, the frogs are stored as a pointer, so we have to copy frogs manually.
		td[i].gpool.frogs = (uint32_t *) malloc(ns->k * 4);
		for (int k=0; k<ns->k; k++){
			td[i].gpool.frogs[k] = gpool.frogs[k];
		}
		for (int j=0; j < i; j++){
			advance_gpool (&td[i].gpool, NULL);
		}
	}
	long sievestart = clock();

	/* Set things in motion */
	for (int i=0; i<nthreads; i++){
		pthread_create (&ns->threads[i], NULL, run_sieve_thread, &td[i]);
	}

	/* Wait for all threads to finish */
	for (int i=0; i<nthreads; i++){
		pthread_join (ns->threads[i], NULL);
	}
	ns->timing.sieve_time = clock() - sievestart;
	printf("\n");

	/* Now proceed with the rest of the factorization in this thread */
	build_matrix (ns);
	solve_matrix (ns);
	ns->timing.total_time = clock() - start;
}

/* Each sieve thread runs this method as its task. When it returns, the thread dies. This method
 * performs sieving until enough relations have been collected. The odd return and parameter types
 * are mandated by the pthreads specification. */
void *run_sieve_thread (void *args){
	thread_data_t *td = (thread_data_t *) args;
	nsieve_t *ns = td->ns;
	
	block_data_t sievedata;		// allocate a sieve block.
	while (ns->nfull + ns->npartial < ns->rels_needed){	// while we don't have enough relations
		/* Allocate, initialize, and generate a new poly group */
		poly_group_t *curr_polygroup = (poly_group_t *) malloc (sizeof (poly_group_t));
		polygroup_init (curr_polygroup, ns);
		generate_polygroup (&td->gpool, curr_polygroup, ns);
		
		/* Loop over the polynomials our group can generate, and sieve them */
		for (int i=0; i < ns->bvals; i++){
			poly_t *curr_poly = (poly_t *) malloc (sizeof (poly_t));
			poly_init (curr_poly);
			generate_poly (curr_poly, curr_polygroup, ns, i);

			sieve_poly (&sievedata, curr_polygroup, curr_poly, ns);
		}
		/* Once our group is done, we can add the relations to the main repository for them inside
		 * the nsieve_t. This method will acquire the lock on the mutex stored in the nsieve_t, so
		 * two threads don't try to do this at the same time */
		add_polygroup_relations (curr_polygroup, ns);

		free (curr_polygroup->ainverses);	// free our precomputed values, we don't need them anymore.

		/* Get the lock, count the partials in the hashtable, and print our status */
		pthread_mutex_lock (&ns->lock);	
		ns->npartial = ht_count (&ns->partials);
		printf("Have %d of %d relations (%d full + %d combined from %d partial); sieved %d polynomials from %d groups. \r", ns->nfull + ns->npartial, ns->rels_needed, ns->nfull, ns->npartial, ns->partials.nentries, ns->info_npoly, ns->info_npg);
		fflush(stdout);
		pthread_mutex_unlock (&ns->lock);
	}
	return NULL;
}

/* Behold - the main method. You knew it was here somewhere. */
int main (int argc, const char *argv[]){
	nsieve_t ns;
	mpz_t n;
	mpz_init (n);

	int pos = 1;
	int nspecd = 0;

	ns.T = -1;
	ns.fb_bound = -1;
	ns.lp_bound = -1;
	ns.M = -1;
	ns.multiplier = -1;
	int nthreads = 1;
	/* Parse command line arguments that override parameters or specify N */
	while (pos < argc){
		if (!strcmp(argv[pos], "-T")){
			ns.T = atof (argv[pos+1]);
			pos ++;
		} else if (!strcmp(argv[pos], "-fbb")){
			ns.fb_bound = atoi (argv[pos+1]);
			pos ++;
		} else if (!strcmp(argv[pos], "-lpb")){
			ns.lp_bound = atoi (argv[pos+1]);
			pos++;
		} else if (!strcmp(argv[pos], "-M")){
			ns.M = atoi (argv[pos+1]);
			pos++;
		} else if (!strcmp(argv[pos], "-np")){
			ns.lp_bound = 1;
		} else if (!strcmp(argv[pos], "-mult")){
			ns.multiplier = atoi (argv[pos+1]);
			pos++;
		} else if (!strcmp(argv[pos], "-threads")){
			nthreads = atoi (argv[pos+1]);
			pos++;
		} else {
			mpz_set_str (n, argv[pos], 10);
			nspecd = 1;
		}
		pos ++;
	}
	if (!nspecd){
		mpz_inp_str (n, stdin, 10);
	}

	printf ("Removing small factors of N by trial division and pollard rho... \n");
	tdiv (n, 32768);
	rho  (n, 65536, 0);
	
	/* If we found all of the factors by trial division or rho, or the cofactor after doing that
	 * is prime, then we're done and we don't need to start the quadratic sieve. */
	printf("Will factor N = ");
	mpz_out_str (stdout, 10, n);
	printf("\n");
	if (mpz_cmp_ui (n, 1) == 0){
		return 0;
	} else if (mpz_probab_prime_p (n, 15)){
		mpz_out_str(stdout, 10, n);
		printf("\n");
		return 0;
	}

	/* Otherwise, start the sieving */
	printf ("Starting the quadratic sieve... \n");
	long start = clock();
	nsieve_init (&ns, n);

	multithreaded_factor (&ns, nthreads);

	ns.timing.total_time = clock() - start;

	printf ("\nTiming summary: \
		 \n\tInitialization:    %ldms \
		 \n\tSieving:           %ldms \
		 \n\tMatbuild + Filter: %ldms \
		 \n\tMatrix solving:    %ldms \
		 \n\tFactor deduction:  %ldms \
		 \n\tTOTAL:             %ldms\n", ns.timing.init_time/1000, ns.timing.sieve_time/1000, ns.timing.filter_time/1000, ns.timing.matsolve_time/1000, ns.timing.facdeduct_time/1000, ns.timing.total_time/1000);

}

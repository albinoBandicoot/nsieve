#ifndef COMMON_H
#define COMMON_H

#include <ctype.h>
#include <stdint.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <pthread.h>
#include <gmp.h>

#define KMAX 12			// the maximum allowable value for k. 
#define BLOCKSIZE 131072	// the size of a sieve block. Entries are 1 byte.

struct relation;	// this is going to be a rel_t. We have to forward-declare it here; there was a
			// cyclical dependency.

/* This struct defines a group of polynomials that share a common 'A' value. */
typedef struct {
	mpz_t a;		// the value of 'A'
	mpz_t *bvals;		// a list of all possible values of 'b' for this group.
	uint32_t gvals[KMAX];	// a list of the primes g_i that were used to produce 'A'
	uint32_t *ainverses;	// the values of a^-1 (mod p) for each p in the factor base. This stays 
				// the same for the whole group. Precomputing these saves a lot of time 
				// in computing the sieve offsets (the x_0 for which p | Q(x_0)).

	struct relation *victim;	// the selected victim for this group.

	/* We store a list of (pointers to) the relations we've accumulated from sieving polynomials 
	 * associated with this group here. Storing them here instead of as global state is helpful for
	 * multithreading concurrency problems, and also organizes things nicely: at the end of sieving
	 * the block, we select our victim and multiply everything else by the victim (which ends up 
	 * just being concatenating factor lists at this point), then add them to matrel_t's in the nsieve_t.
	*/
	uint32_t nrels;		// how many relations we've stored inside this polygroup.
	struct relation **relns;
} poly_group_t;

/* This struct defines a single polynomial. */
typedef struct {
	poly_group_t *group;	// a pointer to the group. This is very handy.
	mpz_t a;
	mpz_t b;
	mpz_t c;
	uint32_t *bmodp;	// for each prime p in the factor base, we compute b mod p. This accelerates the 
				// sieving, since we may do multiple blocks on each polynomial. This array can 
				// (and should) be freed after the polynomial is sieved.

	uint32_t M;		// the number of blocks to sieve for this polynomial. 
} poly_t;

/* This struct defines state for selecting the values of 'g' that are multiplied together to produce the 
 * polynomial group 'A' values. The optimal value for A is about sqrt(2N/M), so if for simplicity we choose
 * g_i of roughly equal size, each one should be about A^(1/k). Thus, we generate a list of potential g 
 * (primes such that (n/g) = 1) that are near this A^(1/k) value. Sometimes these values will be inside
 * the factor base, sometimes not. Both are OK, but we have to be a little careful to avoid them in the sieve.
 *
 * Once we have selected the values, we need a way to iterate through the various combinations of g values.
 * See poly.c for an explanation of the frog-hopping algorithm. 
*/
typedef struct {
	uint32_t *gpool;	// a list of allowable values of G. The idea here is to pick a certain number
				// of primes g such that (n/g) = 1 on either side of center. These will be 
				// close enough to center that whatever subset of k of them we decide to 
				// multiply together to produce A, this will be close enough to the ideal 
				// value of A = sqrt(2N)/M. 
				
	uint32_t ng;		// number of values in the pool
	uint32_t k;		// the same k as in the nsieve_t
	uint32_t *frogs;	// the last used set of g values. A call to advance_gpool will get the next set of values. 
} poly_gpool_t;


/* A simple linked-list structure for defining a list of factors. Used in storing the factor base indices
 * that divided a certain relation. */

typedef struct flentry {
	uint32_t fac;	// this is actually an index into the factor base.
	struct flentry *next;
} fl_entry_t;

/* This structure defines a relation (either full or partial; NOT a 'combined' relation constructed
 * from two partials. Relations in a form suitable for using in the matrix (fulls or combined) are
 * stored in a matrel_t. The matrix building step is effectively the process of converting rel_t's
 * to matrel_t's.
*/
typedef struct relation {
	poly_t  *poly;		// keep track of the polynomial
	int32_t  x;		// where to evaluate the polynomial
	uint32_t cofactor;	// the part of poly(x) that didn't factor over the factor base. This should
				// either be 1 (for full relations) or a prime between fb_bound and lp_bound.
	fl_entry_t *factors;	
	/* storing the factors of the relation is very desirable for several reasons. First it allows 
	 * us to free the temporaries (bmodp, ainverses) after we're done sieving with them; otherwise 
	 * for our accelerated tdiv code we'd have to hold on to all of them. This is as INSANE amount 
	 * of memory for semi-large factorizations (it will eat 80 MB per second if you let it). This is
	 * the 'bad memory usage problem' mentioned in the commit log. Second, we don't have to re-do the 
	 * trial division when we add to the matrix, which was a substantial time drain, especially for partials.
	*/
} rel_t;

/* One of these structs gets associated with each row of the matrix. There are two pointers to rel_t's;
 * r2 is only non-NULL if this is a partial relation. 'row' contains the actual bits of this row of the
 * matrix, packed into 64-bit ints. If we ever implemented double large primes, this struct would have
 * to be modified to use an arbitrary-sized list of relations.
*/
typedef struct {
	rel_t *r1;
	rel_t *r2;
	uint64_t *row;
} matrel_t;

/* Now we have some structs that define a separate-chaining hashtable for storing the partial relations. */

typedef struct htentry{	
	rel_t *rel;
	struct htentry *next;
} ht_entry_t;

typedef struct {
	uint32_t nbuckets;
	uint32_t nentries;	// just for keeping track of how much junk we've stuffed in the hashtable; mostly for observing how fast partials are accumulating.
	ht_entry_t **buckets;
} hashtable_t;


/* Organizes storage of timing data into one place */
typedef struct {
	long init_time;
	long sieve_time;
	long filter_time;
	long matsolve_time;
	long facdeduct_time;
	long total_time;
} time_data_t;

/* The ubiquitous nsieve_t ("ns->" or "nsieve_t" occur on over 300 lines) - contains lots of global
 * data regarding the current factorization. Only one copy of this is ever made; the same one is
 * passed to all of the different sieving threads, so care must be taken not to modify these values
 * without proper protection. 
*/
typedef struct {
	mpz_t N;		// the number to factor
	unsigned char k;	// the number of distinct primes to use to construct polynomial 'A' values
	unsigned short bvals;	// the number of distinct values for 'B' - given by 2^(k-1).
	unsigned int  M;	// the sieve length. 
	float T;		// The sieve threshold will be T * log(lp_bound).

	unsigned int  fb_bound;	// upper bound for the primes in the factor base
	unsigned int  fb_len;	// number of primes in the factor base
	int           extra_rels;	// number of relations more than the factor base size to collect.
	uint32_t      rels_needed;	// this should be equal to fb_len + extra_rels; it is here for convenience.
	uint32_t *fb;		// the actual primes
	uint8_t  *fb_logs;	// approximations to log_2 (p)
	uint32_t *roots;	// the values of sqrt(n) mod p   for each p in the factor base. 
				// Computed once and for all at the beginning.
	
	uint32_t multiplier;	// instead of factoring N, factor kN, for a small squarefree k. This enables 
				// us to pick a factor base that's nicer (more small p s.t. (N/p) = 1).

	uint32_t nfull;		// running count of the number of full relations found thus far.
	uint32_t npartial;	// running count of the number of partial relations found so far. This will be 
				// updated by calling ht_count at the end of each batch of polynomials

	uint32_t row_len;	// number of 64-bit chunks in a row of the matrix.
	matrel_t *relns;	// this is the list of relations, and also constitutes the matrix. 

	uint32_t lp_bound;	// large prime bound. Only relations whose large prime cofactors are smaller
				// than this bound are admitted into the hashtable. 
	hashtable_t partials;	// hashtable for storing partial relations.

	int nthreads;		// number of sieving threads to use.
	pthread_t *threads;	// pointers to the sieving threads
	pthread_mutex_t lock;	// mutex to coordinate updates of global state (mostly adding things to the matrix)


	/* These fields keep track of various properties of the sieving/timing, for informational purposes */

	int info_npoly;		// how many polys we've sieved
	int info_npg;		// how many poly groups we've sieved
	uint32_t tdiv_ct;	// number of relations that underwent trial division (passed the sieve)
	uint64_t sieve_locs;	// how many positions were sieved.

	time_data_t timing;

} nsieve_t;

/* Each thread needs its own gpool; see nsieve.c for comments on how this works */
typedef struct {
	poly_gpool_t gpool;
	nsieve_t *ns;
} thread_data_t;

/* Matrix row functions */
void clear_row (uint64_t *, nsieve_t *);
void flip_bit (uint64_t *, int);
int  get_bit  (uint64_t *, int);
void xor_row  (uint64_t *res, uint64_t *op, int len);
int  rightmost_1 (uint64_t *, int max_i);
int  is_zero_vec (uint64_t *, int len);
void print_row (uint64_t *, int max_i);

/* Hashtable functions */

void ht_init (nsieve_t *);
uint32_t hash_partial (uint32_t cofactor);	
void ht_add (hashtable_t *ht, rel_t *rel);
uint32_t ht_count (hashtable_t *ht);

/* Generic auxillary functions */ 

uint32_t find_root (mpz_t a, uint32_t p);	// finds modular square root of a (mod p)
uint32_t mod (int x, uint32_t p);
uint64_t mpz_get_64 (mpz_t a);
int mpz_fits_64 (mpz_t a);

/* Factor list ops */
void fl_add (rel_t *, uint32_t);
void fl_free (rel_t *);
int  fl_check (rel_t *, nsieve_t *);
void fl_fillrow (rel_t *, uint64_t *row, nsieve_t *ns);
void fl_concat (rel_t *res, rel_t *victim);	// appends victim's list to res's list
void rel_free (rel_t *);
int  rel_check (rel_t *, nsieve_t *);

int  fb_lookup (uint32_t p, nsieve_t *);
#endif

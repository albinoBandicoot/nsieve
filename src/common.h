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

#define KMAX 12
#define BLOCKSIZE 131072

struct relation;	// this is going to be a rel_t.

typedef struct {
	mpz_t a;
	mpz_t *bvals;
	uint32_t gvals[KMAX];
	uint32_t *ainverses;	// the values of a^-1 (mod p) for each p in the factor base. This stays the same for the whole group.

	struct relation *victim;
	uint64_t *victim_factors;
	uint32_t nrels;
	struct relation **relns;	// we store the relations we've collected for this poly group here. This facilitates keeping track of which
				// rels we need to multiply by what, and also will make multithreading easier if we ever decide to do it.
} poly_group_t;

typedef struct {
	poly_group_t *group;
	mpz_t a;
	mpz_t b;
	mpz_t c;
	uint32_t *offsets1;	// the actual offsets of -M.
	uint32_t *offsets2;	// (for the other root)
	uint32_t *bmodp;	// for each prime p in the factor base, we compute b mod p. This accelerates the sieving, since we may do multiple
				// blocks on each polynomial. This array can (and should) be freed after the polynomial is sieved.
	uint32_t M;	// we compute this on a per-polynomial basis. The value we store is rounded to the nearest multiple of BLOCKSIZE.
//	it turns out we don't need to keep istart, since it is (-b/a) - M. This is because b is produced as the solution to a congruence equation mod a, so it must be < a, so (-b/a) is in [-1, 1]. Thus the sieve interval can be taken to always be [-M, M].
} poly_t;

/* This struct defines state for selecting the values of 'g' that are multiplied together to produce the polynomial group 'A' values. 
 * The optimal value for A is about sqrt(2N/M), so if for simplicity we choose g_i of roughly equal size, each one should be about A^(1/k).
 * Thus, we generate a list of potential g (primes such that (n/g) = 1) that are near this A^(1/k) value (called 'center' in the struct).
 * It is likely that these values will be inside the factor base, so no extra work needs to be done to find them.
 *
 * Once we have selected the values, we need a way to iterate through the various combinations of g values. There are |gpool| C k  of them.
 * We initialize k 'frogs' to be the first (lowest) k primes in the gpool. Then, whenever advance_gpool is called, the highest frog (in the
 * k'th position of frogs) hops one slot forward. If it can't (it's at the end of gpool), the next one back tries to hop forwards. This 
 * continues until a frog is found who can hop (Billy). Then, all the frogs higher than Billy (which were jammed up at the top of gpool) are
 * set immediately following Billy. A large enough gpool will be chosen so that we won't run out of combinations. 
*/
typedef struct {
	uint32_t center;	// (sqrt(2N)/M)^(1/k)
	uint32_t *gpool;	// a list of allowable values of G. The idea here is to pick a certain number of primes g such that (n/g) = 1 on
				// either side of center. These will be close enough to center that whatever subset of k of them we decide to multiply
				// together to produce A, this will be close enough to the ideal value of A = sqrt(2N)/M. 
				
	uint32_t ng;		// number of values in the pool
	uint32_t k;		// the same k as in the nsieve_t
	uint32_t *frogs;	// the last used set of g values. A call to advance_gpool will get the next set of values. 
} poly_gpool_t;

typedef struct relation {
	poly_t  *poly;
	int32_t  x;
	uint32_t cofactor;
} rel_t;

typedef struct {	// One of these gets associated with each row in the matrix. r2 is non-null if this row was construted by combining partials.
	rel_t *r1;	// This will have to be changed to a list should we implement the double large prime variant.
	rel_t *r2;
	uint64_t *row;	// the row of the matrix, packed into 64-bit ints.
} matrel_t;	

typedef struct htentry{	// struct for the linked list in the separate-chaining hashtable for storing partial relations
	rel_t *rel;	// the relation. This is a pointer because we'll need the pointers for constructing the matrel_t's.
	struct htentry *next;	// pointer to the next reln in this bucket; NULL if it's the last one.
} ht_entry_t;

typedef struct {
	uint32_t nbuckets;
	uint32_t nentries;	// just for keeping track of how much junk we've stuffed in the hashtable; mostly for observing how fast partials are accumulating.
	ht_entry_t **buckets;
} hashtable_t;

typedef struct {
	long init_time;
	long sieve_time;
	long filter_time;
	long matsolve_time;
	long facdeduct_time;
	long total_time;
} time_data_t;

typedef struct {
	mpz_t N;		// the number to factor
	unsigned char k;	// the number of distinct prime squares to use to construct polynomial 'A' values
	unsigned short bvals;	// the number of distinct values for 'B' - given by 2^(k-1).
	unsigned int  M;	// the sieve length. 
	float T;		// The sieve threshold will be T * log(lp_bound).

	unsigned int  fb_bound;	// upper bound for the primes in the factor base
	unsigned int  fb_len;	// number of primes in the factor base
	int           extra_rels;	// number of relations more than the factor base size to collect.
	uint32_t      rels_needed;	// this should be equal to fb_len + extra_rels; it is here for convenience.
	uint32_t *fb;		// the actual primes
	uint8_t  *fb_logs;	// approximations to log_2 (p)
	uint32_t *roots;	// the values of sqrt(n) mod p   for each p in the factor base. Computed once and for all at the beginning.
	
	uint32_t multiplier;	// instead of factoring N, factor kN, for a small squarefree k. This enables us to pick a factor base that's nicer (more smooth relns).

	uint32_t nfull;		// running count of the number of full relations found thus far.
	uint32_t npartial;	// running count of the number of partial relations found so far. This will be updated by calling ht_count, probably at the end of each batch of polynomials (or maybe after each poly, depending on how many relations are being recovered from each poly vs. group)).

	uint32_t tdiv_ct;	// number of relations that underwent trial division (passed the sieve)
	uint64_t sieve_locs;	// how many positions were sieved.
	uint32_t row_len;	// number of 64-bit chunks in a row of the matrix.
	matrel_t *relns;	// this is the list of relations, and also constitutes the matrix. 

	uint32_t lp_bound;	// large prime bound. Relations with one factor between the top of the factor base and lp_bound are admitted as partials.
	hashtable_t partials;	// hashtable for storing partial relations. This might have to change for double large primes.

	int nthreads;		// number of sieving threads to use.
	pthread_t *threads;	
	pthread_mutex_t lock;

	int info_npoly;
	int info_npg;

	time_data_t timing;

} nsieve_t;

typedef struct {
	poly_gpool_t gpool;	// thread's gpool. All of the gpools will have the same gvalues; the frogs will be staggered, and to get the next group,
				// the frogs will be advanced (nthreads) times. 
	nsieve_t *ns;
} thread_data_t;

/* Matrix row get/set */
void clear_row (uint64_t *, nsieve_t *);
void flip_bit (uint64_t *, int);
int  get_bit  (uint64_t *, int);
void xor_row  (uint64_t *res, uint64_t *op, int len);
int  rightmost_1 (uint64_t *, int max_i);
int  is_zero_vec (uint64_t *, int len);
void print_row (uint64_t *, int max_i);

/* Hashtable functions */

void ht_init (nsieve_t *);			// allocates space for and initializes the hashtable stored in the nsieve_t.
uint32_t hash_partial (uint32_t cofactor);	// hashes the cofactor of a partial relation for determining its bucket in the hashtable. 
void ht_add (hashtable_t *ht, rel_t *rel);	// add rel to the hashtable
uint32_t ht_count (hashtable_t *ht);		// counts the number of full relations that can be made from the partials in the hashtable.

/* Generic auxillary functions */ 

uint32_t find_root (mpz_t a, uint32_t p);	// finds modular square root of a (mod p)
uint32_t mod (int x, uint32_t p);
uint64_t mpz_get_64 (mpz_t a);
int mpz_fits_64 (mpz_t a);
//uint32_t find_root_ui (uint32_t a, uint32_t p);

#endif

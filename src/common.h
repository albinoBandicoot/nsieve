#ifndef COMMON_H
#define COMMON_H

typedef struct {	// represents a relation, poly(istart+x) = 
	poly_t  *poly;
	uint32_t  x;	// offset from poly->istart; this can be unsigned since istart is the lower bound, and all offsets will be >=0. 
	uint32_t cofactor;
} rel_t;

typedef struct {	// One of these gets associated with each row in the matrix. r2 is non-null if this row was construted by combining partials.
	rel_t *r1;	// This will have to be changed to a list should we implement the double large prime variant.
	rel_t *r2;
	uint64_t *row;	// the row of the matrix, packed into 64-bit ints.
} matrel_t;	

typedef struct {	// struct for the linked list in the separate-chaining hashtable for storing partial relations
	rel_t *rel;	// the relation. This is a pointer because we'll need the pointers for constructing the matrel_t's.
	ht_entry_t *next;	// pointer to the next reln in this bucket; NULL if it's the last one.
} ht_entry_t;

typedef struct {
	uint32_t nbuckets;
	ht_entry_t *buckets;
} hashtable_t;

typedef struct {
	mpz_t N;		// the number to factor
	unsigned char k;	// the number of distinct prime squares to use to construct polynomial 'A' values
	unsigned short bvals;	// the number of distinct values for 'B' - given by 2^(k-1).
	unsigned int  M;	// the sieve length. 

	unsigned int  fb_bound;	// upper bound for the primes in the factor base
	unsigned int  fb_len;	// number of primes in the factor base
	int           extra_rels;	// number of relations more than the factor base size to collect.
	uint32_t      rels_needed;	// this should be equal to fb_len + extra_rels; it is here for convenience.
	uint32_t *fb;		// the actual primes
	uint32_t *roots;	// the values of sqrt(n) mod p   for each p in the factor base. Computed once and for all at the beginning.

	uint32_t nfull;		// running count of the number of full relations found thus far.
	uint32_t npartial;	// running count of the number of partial relations found so far. This will be updated by calling ht_count, probably at the end of each batch of polynomials (or maybe after each poly, depending on how many relations are being recovered from each poly vs. group)).

	matrel_t *relns;	// this is the list of relations, and also constitutes the matrix. 
	hashtable_t partials;	// hashtable for storing partial relations. This might have to change for double large primes.

} nsieve_t;


/* Hashtable functions */

void ht_init (nsieve_t *);			// allocates space for and initializes the hashtable stored in the nsieve_t.
uint32_t hash_partial (uint32_t cofactor);	// hashes the cofactor of a partial relation for determining its bucket in the hashtable. 
void ht_add (hashtable_t *ht, rel_t *rel);	// add rel to the hashtable
uint32_t ht_count (hashtable_t *ht);		// counts the number of full relations that can be made from the partials in the hashtable.

/* Generic auxillary functions */ 

uint32_t find_root (mpz_t a, uint32_t p);	// finds modular square root of a (mod p)
uint32_t find_root_ui (uint32_t a, uint32_t p);

#endif

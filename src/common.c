#include "common.h"

/* Matrix row get/set */

// important note: the matrix positions are like this: [63, 62, ... 1, 0][127, 126, ... 65, 64] etc.
// This is so that position p can be found in block p/64 at bit (1 << p % 64).

void clear_row (matrel_t *m, nsieve_t *ns){
	for (int i=0; i < ns->row_len; i++){
		m->row[i] = 0;
	}
}

int get_bit (matrel_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1 << (pos % 64);
	return (m->row[block] & mask) > 0 ? 1 : 0;
}

void flip_bit (matrel_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1 << (pos % 64);
	m->row[block] = m->row[block] ^ mask;
}

/* Hashtable functions */

void ht_init (nsieve_t *ns){			// allocates space for and initializes the hashtable stored in the nsieve_t.
	ns->partials.nbuckets = 5483;	// make some educated choice about this in the future. It should be prime.
	ns->partials.buckets = (ht_entry_t **)(calloc(5483, sizeof(ht_entry_t*)));	// allocate and clear them so they're all NULL to start off.
}

uint32_t hash_partial (uint32_t cofactor){	// hashes the cofactor of a partial relation for determining its bucket in the hashtable. 
	return (cofactor * 263633281) + 135666227;	// random primes.
}

/* Something to think about: do we want to keep the lists sorted? This might be a good idea, esp. since they will be relatively short. */
void ht_add (hashtable_t *ht, rel_t *rel){	// add rel to the hashtable
	uint32_t slot = hash_partial (rel->cofactor) % ht->nbuckets;
	if (ht->buckets[slot] == NULL){
		ht->buckets[slot]->rel = rel;
		ht->buckets[slot]->next = NULL;
		return;
	}
	// this is a work in progress
	ht_entry_t *rover = ht->buckets[slot];
}

uint32_t ht_count (hashtable_t *ht){		// counts the number of full relations that can be made from the partials in the hashtable.
	return 0;
}

/* Generic auxillary functions */ 

uint32_t find_root (mpz_t a, uint32_t p){	// finds modular square root of a (mod p)
	// for now we just do the stupid brute force thing. This will be replaced later with some of the modular
	// exponentiation algorithms (for certain cases), and perhaps the Tonelli-Shanks algorithm for the remaining case.
	mpz_t temp;
	mpz_init (temp);
	uint32_t k = (uint32_t) mpz_mod_ui(temp, a, p);
	mpz_clear (temp);
	return find_root_ui (k, p);
}

uint32_t find_root_ui (uint32_t a, uint32_t p){
	a %= p;
	uint32_t t = 0;
	while ( t*t % p != a && t < p/2 + 1 ){
		t ++;
	}
	if (t*t%p == a) return t;
	return -1;	// this should not happen, since we should only be calling this once we've confirmed (a/p) = 1.
}

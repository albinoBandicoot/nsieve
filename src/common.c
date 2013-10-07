#include "common.h"

/* Matrix row get/set */

// important note: the matrix positions are like this: [63, 62, ... 1, 0][127, 126, ... 65, 64] etc.
// This is so that position p can be found in block p/64 at bit (1 << p % 64).

void clear_row (uint64_t *m, nsieve_t *ns){
	for (int i=0; i < ns->row_len; i++){
		m[i] = 0;
	}
}

int get_bit (uint64_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1 << (pos % 64);
	return (m[block] & mask) > 0 ? 1 : 0;
}

void flip_bit (uint64_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1 << (pos % 64);
	m[block] = m[block] ^ mask;
}

void xor_row (uint64_t *res, uint64_t *op, int len){
	for (int i=0; i<len; i++){
		res[i] = res[i] ^ op[i];
	}
}

int rightmost_1 (uint64_t *m, int max_i){
	for (int i = max_i; i >= 0; i--){
		if (get_bit(m, i) == 1)	return i;
	}
	return -1;
}

int is_zero_vec (uint64_t *m, int len){
	for (int i=0; i<len; i++){
		if (m[i] != 0) return 0;
	}
	return 1;
}
/* Hashtable functions */

void ht_init (nsieve_t *ns){			// allocates space for and initializes the hashtable stored in the nsieve_t.
	ns->partials.nbuckets = 5483;	// make some educated choice about this in the future. It should be prime.
	ns->partials.buckets = (ht_entry_t **)(calloc(5483, sizeof(ht_entry_t*)));	// allocate and clear them so they're all NULL to start off.
}

uint32_t hash_partial (uint32_t cofactor){	// hashes the cofactor of a partial relation for determining its bucket in the hashtable. 
	return (cofactor * 263633281) + 135666227;	// random primes.
}

/* We want to keep the lists sorted. This makes the counting much easier, and is not much extra trouble or computational expense */
void ht_add (hashtable_t *ht, rel_t *rel){	// add rel to the hashtable
	uint32_t slot = hash_partial (rel->cofactor) % ht->nbuckets;
	if (ht->buckets[slot] == NULL){
		ht->buckets[slot]->rel = rel;
		ht->buckets[slot]->next = NULL;
		return;
	}
	// this is a work in progress
	ht_entry_t *rover = ht->buckets[slot];
	ht_entry_t *trailer = rover;
	while (rover != NULL && rel->cofactor < rover->rel->cofactor){
		trailer = rover;
		rover = rover->next;
	}
	// we want to insert right after trailer and before rover.
	ht_entry_t *newentry = (ht_entry_t *) malloc (sizeof(ht_entry_t));
	newentry -> rel = rel;
	newentry -> next = rover;	// this even works if rover is NULL (we are appending to the list)
	trailer->next = newentry;
}

uint32_t ht_count_bucket (ht_entry_t *h){
	if (h == NULL) return 0;
	uint32_t res = 0;
	uint32_t ct = 0;
	uint32_t curr_cofactor = h->rel->cofactor;
	while (1){
		while (h != NULL && h->rel->cofactor == curr_cofactor){
			ct ++;
			h = h->next;
		}
		res += ct - 1;
		if (h == NULL){
			return res;
		}
		curr_cofactor = h->rel->cofactor;
		ct = 0;
	}
	return 0;	// this should never happen.
}

uint32_t ht_count (hashtable_t *ht){		// counts the number of full relations that can be made from the partials in the hashtable.
	uint32_t res = 0;
	for (int i=0; i < ht->nbuckets; i++){
		res += ht_count_bucket (ht->buckets[i]);
	}
	return res;
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

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
	uint64_t mask = 1ull << (pos % 64);
	return (m[block] & mask) > 0 ? 1 : 0;
}

void flip_bit (uint64_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1ull << (pos % 64);
//	printf("pos = %d; block = %d; mask = %llx\n", pos, block, mask);
	m[block] = m[block] ^ mask;
}

void xor_row (uint64_t *res, uint64_t *op, int len){
	for (int i=0; i<len; i++){
		res[i] = res[i] ^ op[i];
	}
}

int rightmost_1 (uint64_t *m, int max_i){
	int block = max_i / 64;
	// first skip over all of the leading 0 blocks.
	while (block >= 0 && m[block] == 0)  block--;
	if (m[block] == 0) return -1;

	// now we need to determine the rightmost 1 in m[block].
	uint64_t x = m[block];
	int pos = 0;
	while ((x & (1ull << 63)) == 0){
		pos++;
		x = x << 1;
	}
	int res =  (63 - pos) + 64 * block;
	return res;
}

int is_zero_vec (uint64_t *m, int len){
	for (int i=0; i<len; i++){
		if (m[i] != 0) return 0;
	}
	return 1;
}

void print_row (uint64_t *m, int max_i){
	for (int i=0; i<max_i; i++){
		if (get_bit(m, i) == 1){
			fprintf(stderr, "1");
		} else {
			fprintf(stderr, "0");
		}
	}
	fprintf(stderr, "\n");
}

/* Hashtable functions */

void ht_init (nsieve_t *ns){			// allocates space for and initializes the hashtable stored in the nsieve_t.
	ns->partials.nbuckets = 5483;	// make some educated choice about this in the future. It should be prime.
	ns->partials.buckets = (ht_entry_t **)(calloc(5483, sizeof(ht_entry_t*)));	// allocate and clear them so they're all NULL to start off.
	ns->partials.nentries = 0;
}

uint32_t hash_partial (uint32_t cofactor){	// hashes the cofactor of a partial relation for determining its bucket in the hashtable. 
	return (cofactor * 263633281) + 135666227;	// random primes.
}

/* We want to keep the lists sorted. This makes the counting much easier, and is not much extra trouble or computational expense */
void ht_add (hashtable_t *ht, rel_t *rel){	// add rel to the hashtable
	uint32_t slot = hash_partial (rel->cofactor) % ht->nbuckets;
//	if (ht->buckets[slot] != NULL) printf("Adding rel to bucket %d with cofactor = %d to bucket. bucket[0] cofactor is %d\n", slot, rel->cofactor, ht->buckets[slot]->rel->cofactor);
	ht->nentries ++;
	if (ht->buckets[slot] == NULL){
//		printf("-- making it as first element\n");
		ht_entry_t *newentry = (ht_entry_t *) malloc(sizeof(ht_entry_t));
		newentry -> rel = rel;
		newentry -> next = NULL;
		ht->buckets[slot] = newentry;
		return;
	}
	ht_entry_t *rover = ht->buckets[slot];
	if (rel->cofactor < rover->rel->cofactor){	// insert at beginning, which is a little different since we have to modify one of the buckets.
//		printf("-- adding to start\n");
		ht_entry_t *newentry = (ht_entry_t *) malloc(sizeof(ht_entry_t));
		newentry -> rel = rel;
		newentry -> next = rover;
		ht->buckets[slot] = newentry;
		return;
	}
	ht_entry_t *trailer = rover;
	rover = rover->next;
	while (rover != NULL && rel->cofactor > rover->rel->cofactor){
		trailer = rover;
		rover = rover->next;
		if (rover != NULL){
//			printf(" - advancing; will now compare against %d\n", rover->rel->cofactor);
		}
	}
	// we want to insert right after trailer and before rover.
//	printf(" -- adding\n");
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

uint32_t find_root (mpz_t k, uint32_t p){	// finds modular square root of k (mod p)
	// for now we just do the stupid brute force thing. This will be replaced later with some of the modular
	// exponentiation algorithms (for certain cases), and perhaps the Tonelli-Shanks algorithm for the remaining case.
	mpz_t pol, temp, temp2;	// the only reason the names are like this is because that's how they were in the place I stole this from.
	mpz_inits (pol, temp, temp2, NULL);

	mpz_mod_ui (temp, k, p);
	uint32_t a = mpz_get_ui (temp);
	uint32_t res = 0;
	// This Pocklington code is stolen from my previous quadratic sieve implementation.
	if (p % 4 == 3){	// use Case 1 of Pocklington's algorithm.
		uint32_t m = p/4;
		mpz_set_ui(pol, p);
		mpz_powm_ui (temp, temp, m+1, pol);
		return mpz_get_ui (temp);
	} else if (p % 8 == 5){
		uint32_t m = p/8;
		mpz_set_ui (pol, p);
		mpz_powm_ui (temp2, temp, 2*m+1, pol);
		if (mpz_cmp_ui(temp2, 1) == 0){
			mpz_powm_ui(temp2, temp, m+1, pol);
			res = mpz_get_ui (temp2);
			mpz_clears (pol, temp, temp2, NULL);
			return res;
		} else {
			mpz_mul_ui(temp, temp, 4);
			mpz_powm_ui(temp2, temp, m+1, pol);
			if (mpz_divisible_ui_p(temp2, 2)){
				res = mpz_get_ui(temp2)/2;
				mpz_clears (pol, temp, temp2, NULL);
				return res;
			} else {
				res = (mpz_get_ui(temp2)+p)/2;
				mpz_clears (pol, temp, temp2, NULL);
				return res;
			}
		}
	}
	// otherwise, do brute force.

	mpz_clears (pol, temp, temp2, NULL);

	uint32_t t = 0;
	while ( t*t % p != a && t < p/2 + 1 ){
		t ++;
	}
	if (t*t%p == a) return t;
	return -1;	// this should not happen, since we should only be calling this once we've confirmed (a/p) = 1.
}

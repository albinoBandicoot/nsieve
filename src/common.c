#include "common.h"

/* Matrix row operations */

// important note: the matrix positions are like this: [63, 62, ... 1, 0][127, 126, ... 65, 64] etc.
// This is so that position p can be found in block p/64 at bit (1 << p % 64).

/* Clears out (zeroes) the row. */
void clear_row (uint64_t *m, nsieve_t *ns){
	for (int i=0; i < ns->row_len; i++){
		m[i] = 0;
	}
}

/* Gets the bit at position 'pos' of the row 'm' */
int get_bit (uint64_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1ull << (pos % 64);
	return (m[block] & mask) > 0 ? 1 : 0;
}

/* Flips the bit at position 'pos' of row 'm' */
void flip_bit (uint64_t *m, int pos){
	int block = pos/64;
	uint64_t mask = 1ull << (pos % 64);
	m[block] = m[block] ^ mask;
}

/* This will take two rows, and xor the second one into the first. 'len' is the number of 64-bit
 * ints in a row. There are assembly routines that take advantage of SSE instructions to do this
 * slightly more efficiently than this code, so this routine is only included in C if USE_ASM hasn't
 * been defined.
*/
#ifndef USE_ASM
void xor_row (uint64_t *res, uint64_t *op, int len){
	for (int i=0; i<len; i++){
		res[i] = res[i] ^ op[i];
	}
}
#endif

/* bitscan is another assembly routine that takes advantage of the BSF instruction, which stands for
 * 'bit scan forwards' - it finds the position of the first 1 bit in the number. The function contains
 * all of 2 or 3 instructions. It will only be called if USE_ASM has been defined.
*/
extern int bitscan (uint32_t x);

/* Find the rightmost set bit in a row of the matrix, starting at column max_i. A value smaller than
 * the length for max_i is specified in some places in the matrix solving code when it is already 
 * known that the rightmost 1 must be to the left of a certain place. */
int rightmost_1 (uint64_t *m, int max_i){
	int block = max_i / 64;
	// first skip over all of the leading 0 blocks.
	while (block > 0 && m[block] == 0)  block--;
	if (m[block] == 0) return -1;

	// now we need to determine the rightmost 1 in m[block].
#ifdef USE_ASM	// then we can use our bitscan() function. Unfortunately, it only works on 32-bit numbers.
	uint32_t low = m[block];
	uint32_t high = m[block] >> 32;
	if (high == 0){
		return (63 - 32 - bitscan(low)) + 64 * block;
	} else {
		return (63 - bitscan(high)) + 64 * block;
	}
#else
	// if we don't have asm, we do it ourselves in C.
	uint64_t x = m[block];
	
	int pos = 0;
	while ((x & (1ull << 63)) == 0){
		pos++;
		x = x << 1;
	}
	int res =  (63 - pos) + 64 * block;
	return res;
#endif
}

/* Tests whether this row is the zero vector. This will alert us to the presence of a linear dependency
 * in the matrix, and hence a congruence of squares */
int is_zero_vec (uint64_t *m, int len){
	for (int i=0; i<len; i++){
		if (m[i] != 0) return 0;
	}
	return 1;
}

/* Debugging routine for printing out a matrix row in binary */
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

/* A hashtable is used to store the partial relations. Partials are hashed based on their cofactor,
 * so partials with the same cofactor end up in the same bin. When we add relations to the buckets,
 * we make sure to add them in sorted order, so all of the partials with the same cofactor are clustered
 * together. This practice makes coalescing the partials into full relations, and counting how many full
 * relations can be made out of the partials much less painful.
*/

/* Allocate storage for the hashtable. Currently it's just a fixed (prime) size of 5483; in the future
 * some educated choice should be made about this size. However, it's actually not all that important
 * since not that much time gets spent doing hashtable operations. */
void ht_init (nsieve_t *ns){
	ns->partials.nbuckets = 5483;
	ns->partials.buckets = (ht_entry_t **)(calloc(5483, sizeof(ht_entry_t*)));
	ns->partials.nentries = 0;
}

// hashes the cofactor of a partial relation for determining its bucket in the hashtable. 
uint32_t hash_partial (uint32_t cofactor){
	return (cofactor * 263633281) + 135666227;
}

/* We want to keep the lists sorted. This makes the counting much easier, and is not much extra trouble or computational expense */
void ht_add (hashtable_t *ht, rel_t *rel){	// add rel to the hashtable
	uint32_t slot = hash_partial (rel->cofactor) % ht->nbuckets;
	ht->nentries ++;
	if (ht->buckets[slot] == NULL){
		ht_entry_t *newentry = (ht_entry_t *) malloc(sizeof(ht_entry_t));
		newentry -> rel = rel;
		newentry -> next = NULL;
		ht->buckets[slot] = newentry;
		return;
	}
	ht_entry_t *rover = ht->buckets[slot];
	if (rel->cofactor < rover->rel->cofactor){	// insert at beginning, which is a little different since we have to modify one of the buckets.
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
	}
	// we want to insert right after trailer and before rover.
	ht_entry_t *newentry = (ht_entry_t *) malloc (sizeof(ht_entry_t));
	newentry -> rel = rel;
	newentry -> next = rover;	// this even works if rover is NULL (we are appending to the list)
	trailer->next = newentry;
}

/* Counts the number of full relations that can be made out of the partials in this bucket. Since one
 * partial with each cofactor must be sacrificed to the factoring gods as a victim so that the others
 * may ascend into full relation status, only d-1 full relations can be produced from d partials that
 * share a cofactor.
*/
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

/* counts the number of full relations that can be made from the partials in the hashtable. This is
 * just the sum of all of the bucket counts. */
uint32_t ht_count (hashtable_t *ht){		
	uint32_t res = 0;
	for (int i=0; i < ht->nbuckets; i++){
		res += ht_count_bucket (ht->buckets[i]);
	}
	return res;
}

/* Generic auxillary functions */

#define POCKLINGTON
#define TONELLI_SHANKS

/* Finds the modular square root of k (mod p). This will use either Pocklington's algorithm, if
 * p is congruent to 3 (mod 4) or 5 (mod 8), and the Tonelli-Shanks algorithm otherwise. What 
 * algorithms are used can be changed by the #defines above; if one of those symbols is not defined,
 * that algorithm is not used. There is a brute-force basecase implemented, so something that
 * produces a correct result (eventually) will always be there.
*/
uint32_t find_root (mpz_t k, uint32_t p){	// finds modular square root of k (mod p)
	mpz_t pol, temp, temp2;	// the only reason the names are like this is because that's how they 
				// were in the place I stole this from (my old QS code).
	mpz_inits (pol, temp, temp2, NULL);

	mpz_mod_ui (temp, k, p);
	uint32_t a = mpz_get_ui (temp);
	if (p == 2){
		return a;
	}
	uint32_t res = 0;
#ifdef POCKLINGTON
	// This Pocklington code is stolen from my previous quadratic sieve implementation.
	if (p % 4 == 3){	// use Case 1 of Pocklington's algorithm.
		uint32_t m = p/4;
		mpz_set_ui(pol, p);
		mpz_powm_ui (temp, temp, m+1, pol);
		return mpz_get_ui (temp);
	} else if (p % 8 == 5){	// Case 2
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

#endif

#ifdef TONELLI_SHANKS
	/* Based on algorithm description on programmingpraxis.com/2012/11/23/tonelli-shanks-algorithm/ */
	/* Write p as s * 2^e + 1, for the largest possible such e. We seek sqrt(a) (mod p). 
	 * Find n such that n^((p-1)/2) ~= -1 (mod p)	(do this by trial and error until one is found, starting with n=2)
	 * Then set x = a^((s+1)/2) % p,
	 * 	    b = a^s % p
	 * 	    g = n^s % p
	 * 	    r = e.
	 * (1) Now find the smallest m >= 0 such that b^2^m = 1 (mod p). 
	 * If m = 0, the algorithm terminates: output x. Otherwise set
	 * 	x = x * g^2^(r-m-1) (mod p)
	 * 	b = b * g^2^(r-m)   (mod p)
	 * 	g = g^2^(r-m)	    (mod p)
	 * 	r = m
	 * and goto (1).
	*/
	uint32_t s = p-1;
	uint32_t e = 0;
	while (s % 2 == 0){
		s /= 2;
		e ++;
	}
	mpz_set_ui (pol, p);	// this will just hold a mpz version of p.
	// find a suitable n...
	mpz_set_ui(temp2, 2);
	mpz_powm_ui (temp, temp2, (p-1)/2, pol);
	while (mpz_cmp_ui (temp, p-1) != 0){
		mpz_add_ui (temp2, temp2, 1);
		mpz_powm_ui (temp, temp2, (p-1)/2, pol);
	}
	// temp2 now contains a suitable n.
	
	mpz_t x, b, g;
	mpz_inits (x, b, g, NULL);
	mpz_set_ui (temp, a);
	mpz_powm_ui (x, temp, (s+1)/2, pol);	// x = a^((s+1)/2) mod p
	mpz_powm_ui (b, temp, s, pol);		// b = a^s mod p
	mpz_powm_ui (g, temp2, s, pol);		// g = n^s mod p
	uint32_t r = e;
	
	// now the main loop begins. We will use temp2 to hold b^2^m, since we don't need n anymore.
	uint32_t m;
tshanks_mainloop:
	m = 0;
	mpz_set (temp, b);	// temp = b^2^0 = b^1 = b
	while (mpz_cmp_ui (temp, 1) != 0){
		mpz_mul (temp, temp, temp);
		mpz_mod (temp, temp, pol);
		m++;
	}
	if (m == 0){
		// we're done; output x.
		uint32_t res = mpz_get_ui (x);
		mpz_clears (x, b, g, pol, temp, temp2, NULL);
		return res;
	} else {
		// x *= g^2^(r-m-1)
		// b *= g^2^(r-m)
		// g  = g^2^(r-m)
		// r  = m
		mpz_ui_pow_ui (temp, 2, r-m-1);
		mpz_powm (temp, g, temp, pol);	// temp = g^2^(r-m-1)
		mpz_mul (x, x, temp);		// multiply in to x
		mpz_mod (x, x, pol);
		mpz_mul (temp, temp, temp);	// temp now = g^2^(r-m)
		mpz_mod (temp, temp, pol);	// reduce mod p
		mpz_mul (b, b, temp);		// multiply in to b
		mpz_mod (b, b, pol);
		mpz_set (g, temp);		// set g
		r = m;
	}
	goto tshanks_mainloop;

#endif
	mpz_clears (pol, temp, temp2, NULL);

/* If we're here, none of the other algorithms were applicable / enabled for this value of p, so 
 * we proceed with the brute force method: try each t < p/2; if t*t == a (mod p), output t. */
	uint32_t t = 0;
	while ( t*t % p != a && t < p/2 + 1 ){
		t ++;
	}
	if (t*t%p == a) return t;
	return -1;	// this should not happen, since we should only be calling this once we've confirmed (a/p) = 1.
}

// compute x (mod p), according to the mathematical definition.
inline uint32_t mod (int32_t x, uint32_t p){	
	if (x > 0){
		return x % p;
	} else {
		uint32_t z = p + (x % p);
		if (z == p) return 0;
		return z;
	}
}

/* Since a long int is not 64 bits on all systems, we implement two routines that mimic GMP syntax
 * to test if a mpz_t fits in 64 bits, and to return a uint64_t if it does. This is useful in
 * accelerating the trial division code, because we would like to switch to 64-bit arithmetic once
 * the values become small enough.
*/

int mpz_fits_64 (mpz_t a){
	return mpz_sizeinbase (a, 2) < 63;	// play this safe.
}

uint64_t mpz_get_64 (mpz_t a){
	if (mpz_fits_ulong_p (a)){
		return (uint64_t) mpz_get_ui(a);
	}
	uint64_t res = 0;
	if (sizeof(mp_limb_t) == 4){
		res = mpz_getlimbn (a, 1);
		res = (res << 32) | mpz_getlimbn(a, 0);
	} else {
		res = mpz_getlimbn (a, 0);
	}
	return res;
}

/* Factor list ops */

/* Add a prime to the factor list of 'rel' */
void fl_add (rel_t *rel, uint32_t p){
	// prepend, because that is more efficient. It does mean that the list will be in reverse
	// order, which should not matter at all.
	fl_entry_t *entry = (fl_entry_t *) malloc(sizeof (fl_entry_t));
	entry->fac = p;
	entry->next = rel->factors;
	rel->factors = entry;
}

/* Free a factor list */
void fl_free (rel_t * rel){
	fl_entry_t *entry = rel->factors;
	while (entry != NULL){
		fl_entry_t *next = entry->next;
		free (entry);
		entry = next;
	}
}

/* A safety check for factor lists. Sometimes, for some reason, primes less than the factor base
 * bound but not in the factor base divide a couple of sieve values. Mathematically, this should
 * be impossible. If this happens, a very large value is stored in the factor list to represent it 
 * (it is actually -p, interpreted as unsigned).
*/
int fl_check (rel_t *rel, nsieve_t *ns){
	fl_entry_t *entry = rel->factors;
	while (entry != NULL){
		if (entry->fac > ns->fb_len){
			return 0;
		}
		entry = entry->next;
	}
	return 1;
}

/* Walk down the factor list, flipping the bit of the matrix row 'row' corresponding to the entries
 * in the list. This is used to build the matrix from the relations. Note that it clears the row 
 * first, so multiplying in additional rows should be done via xors.
*/
void fl_fillrow (rel_t *rel, uint64_t *row, nsieve_t *ns){
	clear_row (row, ns);
	fl_entry_t *entry = rel->factors;
	while (entry != NULL){
		if (entry->fac > ns->fb_len+1){	// the same check as in fl_check, repeated here.
			printf ("WARNING - bad factor: %d\n", entry->fac);
		}
		flip_bit (row, entry->fac);
		entry = entry->next;
	}
}

/* Concatenate two factor lists. Very useful for multiplying relations by victims: we can just stick
 * the two factor lists together.
*/
void fl_concat (rel_t *res, rel_t *victim){
	fl_entry_t *entry = res->factors;
	if (entry == NULL){
		res->factors = victim->factors;
	} else {
		while (entry->next != NULL){
			entry = entry->next;
		}
		entry->next = victim->factors;
	}
}

void rel_free (rel_t *rel){
	fl_free (rel);
	free (rel);
}

/* Do a binary search on the factor base for the prime 'p.' The index returned is 1 more than the
 * index of the prime in the factor base - it is actually returning a matrix row position. Since
 * -1 occupies the first position, everything is shifted.
 *
 * If the prime is not found, it returns -p, which, since it will be stored as unsigned, will be a
 * very large number. This can be used to check the factor lists for validity later, while giving
 * us some info about what went wrong.
*/
int fb_lookup (uint32_t p, nsieve_t *ns){
	// returns the index of p in the factor base, plus 1 (for indexing into the matrix rows).
	int low = 0;
	int high = ns->fb_len - 1;
	int guess = (low + high)/2;
	while (high - low > 1){
		if (p < ns->fb[guess]){
			high = guess;
			guess = (high + low)/2;
		} else if (p > ns->fb[guess]){
			low = guess;
			guess = (high + low)/2;
		} else {
			return guess + 1;
		}
	}
	for (int i=low; i<high; i++){
		if (p == ns->fb[i]){
			return i+1;
		}
	}
	return -p;	// this is so we have some info about what went wrong.
}


#include "filter.h"

void build_matrix (nsieve_t * ns){
}

/* Whenever D partial relations share a cofactor, we can build D-1 full relations from them as follows:
 * We pick one of them, and multiply the others by it (does this sound familiar?) - we end up with the
 * matrix row side of the equation having a cofactor^2 term. Since we want to put this side in the matrix, but
 * keep the info about the cofactor being there, we multiply the LHS by (cofactor ^ -1 (mod N))^2. 
 *
 * Note that at this point, the 'victim' that's used to cancel the polynomial A values  has not been folded 
 * into the partial relations yet, so we have to do that as well.
 *
*/

void combine_bucket (ht_entry_t *h, nsieve_t *ns){
	if (h == NULL) return;
	while (1){
		rel_t *base_rel = h->rel;
		uint64_t base_factors[ns->row_len];
		fb_factor_rel (base_rel, &base_factors[0], ns);

		h = h->next;	// skip past the base rel.
		while (h != NULL && base_rel->cofactor == h->rel->cofactor){
//			printf("Making combined rel with base cof = %d and h cof = %d\n", base_rel->cofactor, h->rel->cofactor);
			// any h inside this loop creates another combined relation.
			matrel_t *m = &ns->relns[ns->nfull];
			m -> r1 = base_rel;
			m -> r2 = h -> rel;
			fb_factor_rel (h->rel, m->row, ns);
			xor_row (m->row, &base_factors[0], ns->row_len);	// multiply the factorizations together.
			xor_row (m->row, base_rel->poly->group->victim_factors, ns->row_len);	// and since this hasn't been done yet, multiply
			xor_row (m->row, m->r2->poly->group->victim_factors, ns->row_len);		// in their corresponding victims.
			ns->nfull ++;
			if (ns -> nfull >= ns -> rels_needed){
				return;
			}

			h = h->next;
		}
		if (h == NULL){
		       return;
		}
	}
}

void combine_partials (nsieve_t *ns){
	for (int i=0; i < ns->partials.nbuckets; i++){
		if (i % 100 == 0){
			printf("Combining bucket %d\r", i);
			fflush(stdout);
		}
		if (ns->nfull >= ns->rels_needed){
			return;
		}
		combine_bucket (ns->partials.buckets[i], ns);
	}
}

void filter (nsieve_t *ns){
}

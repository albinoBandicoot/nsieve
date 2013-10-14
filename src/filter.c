#include "filter.h"

void build_matrix (nsieve_t * ns){
	for (int i=0; i < ns->nfull; i++){
		ns->relns[i].row = (uint64_t *) calloc (ns->row_len, 8);
		fl_fillrow (ns->relns[i].r1, ns->relns[i].row, ns);
	}
	combine_partials (ns);
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
	while (1){
		if (h == NULL) return;
		if (h->rel->poly->group->victim == NULL){
			h = h->next;
			continue;
		}
		rel_t *base_rel = h->rel;
		uint64_t base_factors[ns->row_len];
		fl_concat (base_rel, base_rel->poly->group->victim);
		fl_fillrow (base_rel, &base_factors[0], ns);
//		fb_factor_rel (base_rel, &base_factors[0], ns);

		h = h->next;	// skip past the base rel.
		while (h != NULL && base_rel->cofactor == h->rel->cofactor){
//			printf("Making combined rel with base cof = %d and h cof = %d\n", base_rel->cofactor, h->rel->cofactor);
			// any h inside this loop creates another combined relation.
			if (h->rel->poly->group->victim == NULL){
				h = h->next;
				continue;	// we can't use this relation; there were no full relations in its pg to use as the victim.
			}
			matrel_t *m = &ns->relns[ns->nfull];
			m->row = (uint64_t *) calloc (ns->row_len, 8);
			m -> r1 = h -> rel;
			m -> r2 = base_rel;
			fl_concat (h->rel, h->rel->poly->group->victim);
			fl_fillrow (h->rel, m->row, ns);
			xor_row (m->row, &base_factors[0], ns->row_len);	// multiply the factorizations together.
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
	long start = clock();
	for (int i=0; i < ns->partials.nbuckets; i++){
		if (i % 100 == 0){
			printf("Combining bucket %d\r", i);
			fflush(stdout);
		}
		if (ns->nfull >= ns->rels_needed){
			ns->timing.filter_time = clock() - start;
			return;
		}
		combine_bucket (ns->partials.buckets[i], ns);
	}
	ns->timing.filter_time = clock() - start;
}

void filter (nsieve_t *ns){
}

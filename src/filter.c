#include "filter.h"

/* Allocates the actual matrix rows, and fills them in. Also makes a call to combine_partials, which,
 * predictably, combines the partials and adds them into the matrix */
void build_matrix (nsieve_t * ns){
	for (int i=0; i < ns->nfull; i++){
		ns->relns[i].row = (uint64_t *) calloc (ns->row_len, 8);
		fl_fillrow (ns->relns[i].r1, ns->relns[i].row, ns);
	}
	combine_partials (ns);
}

/* Whenever D partial relations share a cofactor, we can build D-1 full relations from them by picking
 * one as a victim (this is a similar concept but a different application of what we did to cancel out
 * polynomial 'A' values - each partial also has those A-victims associated with it) and multiply all
 * of the others by it. The cofactor then becomes a square in all of the resulting relations.
 *
*/

/* This routine will combine all of the relations it can in one bucket. */
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

		h = h->next;	// skip past the base rel.
		while (h != NULL && base_rel->cofactor == h->rel->cofactor){
			// any h inside this loop creates another combined relation.
			if (h->rel->poly->group->victim == NULL){
				h = h->next;
				continue;	// we can't use this relation; there were no full relations in its pg to use as the victim.
						// such 'orphaned partials' should not occur; they should never be added
						// to the hashtable, so this is more of a sanity check than anything. It
						// would probably segfault later if one somehow crept in.
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

/* Combine all of the partials. This just loops through the buckets, combining partials until
 * the matrix is full. */
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

/* Filtering goes here */
void filter (nsieve_t *ns){
}

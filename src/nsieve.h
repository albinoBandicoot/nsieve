#ifndef NSIEVE_H
#define NSIEVE_H
#include "gmp.h"
#include <stdlib.h>
#include "common.h"

void generate_fb (nsieve_t *);	// fills in 'fb' and 'roots'

void nsieve_init (nsieve_t *, mpz_t n);		// initialize all of the other parameters, given only N. 
void factor (nsieve_t *);		// the main top-level routine.

#endif
#endif

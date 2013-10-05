#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

int main (int argc, const char *argv[]){
	gmp_randstate_t rand;
	gmp_randinit_default (rand);
	int asize, bsize;
	asize = atoi(argv[1]);
	bsize = atoi(argv[2]);
	mpz_t a, b;
	mpz_inits (a, b, NULL);
	mpz_urandomb (a, rand, asize);
	mpz_urandomb (b, rand, bsize);
	mpz_mul (a, a, b);
	mpz_out_str (stdout, 10, a);
	printf("\n");
}

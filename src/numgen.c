#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

void getprime (mpz_t res, gmp_randstate_t rand, int bits){
	mpz_urandomb (res, rand, bits);
	while (mpz_sizeinbase (res, 2) != bits){
		mpz_urandomb (res, rand, bits);
	}
	mpz_nextprime (res, res);
}

int main (int argc, const char *argv[]){
	if (argc == 1){
		printf("! Error: expecting command line options (list of integers representing bit lengths of the primes to multiply together)\n");
		return 1;
	}
	mpz_t res, temp;
	mpz_inits (res, temp, NULL);
	mpz_set_ui (res, 1);

	gmp_randstate_t rand;
	gmp_randinit_default (rand);

	int arg =1;
	while (arg < argc){
		getprime (temp, rand, atoi(argv[arg]));
		mpz_mul (res, res, temp);
		arg++;
	}
	mpz_out_str (stdout, 10, res);
	printf("\n");
}

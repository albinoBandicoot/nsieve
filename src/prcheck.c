#include <stdio.h>
#include <stdlib.h>
#include <gmp.h>

int main (int argc, const char *argv[]){
	mpz_t n;
	mpz_init(n);
	if (argc >= 2){
		mpz_set_str(n, argv[1], 10);
	} else {
		mpz_inp_str(n, stdin, 10);
	}

	if (mpz_probab_prime_p (n, 30)){
		printf("N is probably prime\n");
	} else {
		printf("N is composite\n");
	}
	mpz_clear(n);
}

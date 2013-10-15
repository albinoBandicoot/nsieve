#include "rho.h"

int main (int argc, const char *argv[]){
	mpz_t n;
	mpz_init (n);
	if (argc >= 2){
		mpz_set_str (n, argv[1], 10);
	} else {
		mpz_inp_str (n, stdin, 10);
	}
	tdiv (n, 5000);
	if (mpz_probab_prime_p (n, 14)){
		mpz_out_str (stdout, 10, n);
		printf("\n");
		return 0;
	}
	rho  (n, 0xffffffff, 1);
	mpz_clear (n);
}


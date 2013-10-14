#include <stdio.h>
#include <gmp.h>

int main (int argc, const char *argv[]){
	mpz_t n;
	mpz_init(n);
	if (argc >= 2){
		mpz_set_str(n, argv[1], 10);
	} else {
		mpz_inp_str(n, stdin, 10);
	}
	while (mpz_divisible_ui_p(n, 2)){
		mpz_divexact_ui(n, n, 2);
		printf("2\n");
	}
	long int t = 3;
	if (mpz_probab_prime_p(n, 12)){
		mpz_out_str(stdout, 10, n);
		printf("\n");
		return 0;
	}
	while (mpz_cmp_ui(n, 1) != 0){
		while (mpz_divisible_ui_p(n, t)){
			mpz_divexact_ui(n, n, t);
			printf("%ld\n", t);
			if (mpz_probab_prime_p(n, 12)){
				mpz_out_str(stdout, 10, n);
				printf("\n");
				return 0;
			}
		}
		t += 2;
	}
}


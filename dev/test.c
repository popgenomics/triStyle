#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>
#include <gsl/gsl_statistics.h>
#include <gsl/gsl_permutation.h>

int main(int argc, char* argv[]){
        const gsl_rng_type *T;
        gsl_rng *r;
        gsl_rng_env_setup();

        T = gsl_rng_default;
        r = gsl_rng_alloc(T);
        gsl_rng_set(r, 13);

	int i = 0;
	int j = 0;
	int N = 10;
	int* table = NULL;
	table = malloc( N * sizeof(int));
	
	for(i=0; i<N; ++i){
		table[i] = i;
	}	

	for(i=0; i<1000; i++){
		j=gsl_rng_uniform_int(r, N);
		printf("%d\n",table[j]);
	}
	free(table);
}


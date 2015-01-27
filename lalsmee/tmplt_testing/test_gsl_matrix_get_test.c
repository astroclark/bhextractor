#include <gsl/gsl_matrix.h>
#include <gsl/gsl_vector.h>
#include <gsl/gsl_blas.h>
#include <complex.h>
#include <stdio.h>

static gsl_matrix_complex* get_matrix_from_file(const char *file_name, int M, int N){

	FILE *fp;

	gsl_matrix_complex *matrix_from_file = gsl_matrix_complex_calloc(M, N);

	fp = fopen(file_name, "rb");

	gsl_matrix_complex_fread(fp, matrix_from_file); //NOTE: matrix_from_file must be in binary format

	fprintf(stderr, "read in matrix\n");

	fclose(fp);

	return matrix_from_file;	
}

int main(int argc, char **argv){


    int i, M, N;
    char *pcfile = NULL;
    gsl_complex result;
	gsl_matrix_complex *pcs_from_file = NULL;

    pcfile = argv[1];
    M = atoi(argv[2]);
    N = atoi(argv[3]);


	pcs_from_file = get_matrix_from_file(pcfile, M, N);

	/* do some simple operations on columns and rows */	
    for (i=0; i<10; i++){

        result = gsl_matrix_complex_get(pcs_from_file, i, 0);

        fprintf(stdout, "%e +i%e\n", GSL_REAL(result), GSL_IMAG(result));
    }

	gsl_matrix_complex_free(pcs_from_file);

	return 0;
}




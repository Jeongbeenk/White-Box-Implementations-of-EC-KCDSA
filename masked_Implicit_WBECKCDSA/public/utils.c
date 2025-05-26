#include "utils.h"

void init_array(mpz_t *arr, int len) {
    for( int i = 0; i < len; i++ ) {
        mpz_init(arr[i]);
    }
}

void clear_array(mpz_t *arr, int len) {
    for( int i = 0; i < len; i++ ) {
        mpz_clear(arr[i]);
    }
}

void init_table(mpz_t *table, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            mpz_init(table[i * col + j]);
        }
    }
}

void clear_table(mpz_t *table, int row, int col) {
    for( int i = 0; i < row; i++ ) {
        for( int j = 0; j < col; j++ ) {
            mpz_clear(table[i * col + j]);
        }
    }
}

void init_table_3(mpz_t *table, int dim1, int dim2, int dim3) {
    for( int i = 0; i < dim1; i++ ) {
        for( int j = 0; j < dim2; j++ ) {
            for( int k = 0; k < dim3; k++ ) {
                mpz_init(table[i * dim2 * dim3 + j * dim3 + k]);
            }
        }
    }
}

void clear_table_3(mpz_t *table, int dim1, int dim2, int dim3) {
    for( int i = 0; i < dim1; i++ ) {
        for( int j = 0; j < dim2; j++ ) {
            for( int k = 0; k < dim3; k++ ) {
                mpz_clear(table[i * dim2 * dim3 + j * dim3 + k]);
            }
        }
    }
}

int get_two_bits(mpz_t num, int bit_position) {
    int bit1 = mpz_tstbit(num, bit_position);
    int bit2 = mpz_tstbit(num, bit_position+1);

    return (bit2 << 1) | bit1;
}


void print_matrix(mpz_t* Mat, int row, int col) {
    for (int i = 0; i < row; i++) {
        for (int j = 0; j < col; j++) {
            gmp_printf("%064ZX ", Mat[i * col + j]);
        }
        printf("\n");
    }
    printf("\n");
}


void print_array(mpz_t* arr, int len) {
    for( int i=0; i < len; i++ ) {
        gmp_printf("%Zx\n", arr[i]);
    }
    puts("");
}


void init_eq_coeff_table(struct EQ_COEFF_TABLE *table){
    for( int i=0; i<127; i++ ) {
        for( int j=0; j<4; j++ ) {
            for( int k=0; k<12; k++ ) {
                mpz_init(table->coeff_m[i][j][k]);
            }
        }
    }

    for( int i=0; i<3; i++ ) {
        for( int j=0; j<44; j++ ) {
            mpz_init(table->coeff_m_last[i][j]);
        }
    }
}

void clear_eq_coeff_table(struct EQ_COEFF_TABLE *table) {
    for( int i=0; i<127; i++ ) {
        for( int j=0; j<4; j++ ) {
            for( int k=0; k<12; k++ ) {
                mpz_clear(table->coeff_m[i][j][k]);
            }
        }
    }

    for( int i=0; i<3; i++ ) {
        for( int j=0; j<44; j++ ) {
            mpz_clear(table->coeff_m_last[i][j]);
        }
    }
}


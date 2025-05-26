#ifndef ENCODE_GEN_H
#define ENCODE_GEN_H

#include "../public/params.h"
#include "secret_params.h"
#include "../public/EC.h"
#include "../public/utils.h"


mpz_t** init_matrix(int size);
void clear_matrix(mpz_t **Matrix, int size);
// void print_matrix(mpz_t Matrix[4][4]);

#define SIZE 4
void cal_determinant_4x4(mpz_t det, mpz_t matrix[4][4], int n, mpz_t mod);

void generate_random_matrix_3x3(mpz_t matrix[3][3], gmp_randstate_t state, mpz_t mod);
void generate_random_matrix_4x4(mpz_t matrix[4][4], gmp_randstate_t state, mpz_t mod);
void gen_random_invertible_matrix_3x3(mpz_t matrix[3][3], gmp_randstate_t state, mpz_t mod);
void gen_random_invertible_matrix_4x4(mpz_t matrix[4][4], gmp_randstate_t state, mpz_t mod);
void get_inverse_matrix_4x4(mpz_t inv[4][4], mpz_t matrix[4][4], mpz_t mod);
void get_inverse_matrix_3x3(mpz_t inv[3][3], mpz_t matrix[3][3], mpz_t mod);

void gen_random_vector_3(mpz_t vec[3], gmp_randstate_t state, const mpz_t mod);
void gen_random_vector_4(mpz_t vec[4], gmp_randstate_t state, const mpz_t mod);

void init_identity_matrix_3x3(mpz_t Matrix[3][3]);
void init_identity_matrix_4x4(mpz_t Matrix[4][4]);


void tmp_k_d(mpz_t d_i_table[128][4],mpz_t k_i_table[128][4]);
void make_and_write_ki_di_alpha_beta_table(const char* filename, mpz_t d_i_table[128][4], mpz_t alpha_i_table[128][4], mpz_t beta_i_table[128][4],mpz_t k_i_table[128][4]);
void cal_ki_G_with_table(const char* ki_str[128][4], Point G, mpz_t e, mpz_t table[128][4][2]);
void write_ki_G_table(const char* filename, mpz_t table[128][4][2]);

void gen_eq_coeff_table( mpz_t m[128][4][4], mpz_t a[128][4][4], mpz_t b[128][4], mpz_t pre_a[128][4][4], mpz_t pre_b[128][4], 
                        mpz_t m_last[3][3], mpz_t a_last[3][3], mpz_t b_last[3],
                        mpz_t ki_table[128][4], mpz_t alpha_i_table[128][4], mpz_t beta_i_table[128][4], mpz_t di_table[128],
                        struct EQ_COEFF_TABLE *eq_coeff_table,  mpz_t n);

void write_eq_coeff_table(const char* filename, struct EQ_COEFF_TABLE *eq_coeff_table);
void write_last_encoding_A(const char* filename, mpz_t table[3][3]);
void write_last_encoding_B(const char* filename, mpz_t array[3]);
void write_u_input(const char* filename, mpz_t array[4]);

#endif

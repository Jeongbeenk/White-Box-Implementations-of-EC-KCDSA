#ifndef GE_H
#define GE_H

#include "utils.h"
#include "params.h"


void gauss_elimination_3x3_mod(mpz_t Mat[3][4], mpz_t y[3], const mpz_t p);
void gauss_elimination_4x4_mod(mpz_t Mat[4][5], mpz_t y[4], const mpz_t p);

void matrix_update_from_masking(mpz_t ge_matrix[4][5], mpz_t u[4], mpz_t coeff_m_table[4][12], mpz_t hash_e, mpz_t n, int round);
void matrix_update_from_masking_128(mpz_t ge_matrix[3][4], mpz_t u[4], mpz_t coeff_m_table[3][44], mpz_t hash_e, mpz_t rho, mpz_t n);

void encrypt(mpz_t u[4], mpz_t v[4], mpz_t res[3], struct EQ_COEFF_TABLE eq_coeff_table, mpz_t hash_e, mpz_t rho, mpz_t n);

#endif

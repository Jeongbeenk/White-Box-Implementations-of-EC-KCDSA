#ifndef UTILS_H
#define UTILS_H

#include <stdio.h>
#include <gmp.h>
#include <stdlib.h>
#include "params.h"

int get_two_bits(mpz_t num, int bit_position);

void init_array(mpz_t *arr, int len);
void clear_array(mpz_t *arr, int len);

void init_table(mpz_t *table, int row, int col);
void clear_table(mpz_t *table, int row, int col);

void init_table_3(mpz_t *table, int dim1, int dim2, int dim3);
void clear_table_3(mpz_t *table, int dim1, int dim2, int dim3);

void print_matrix(mpz_t* Mat, int row, int col);
void print_array(mpz_t* arr, int len);

void init_eq_coeff_table(struct EQ_COEFF_TABLE *table);
void clear_eq_coeff_table(struct EQ_COEFF_TABLE *table);

#endif

#ifndef EC_H
#define EC_H
/** 타원 곡선에 대한 연산을 정의 **/

// #include "params.h"
#include <stdio.h>
#include <gmp.h>

typedef struct point {
    mpz_t x, y;
}__Point;
typedef __Point Point[1];

void point_init(Point P);
void point_clear(Point P);

void point_init_set_str(Point P, const char * x_str, const char * y_str, int base);

void point_init_infinity(Point P);
int point_is_infinity(Point P);
int point_equal(Point P, Point Q);
int point_is_inverse(Point P, Point Q);

void point_set(Point R, Point P);

void point_add(Point R, Point P, Point Q, mpz_t z, mpz_t p);

void point_scalar(Point R, Point P, mpz_t scalar, mp_bitcnt_t num_bits, mpz_t a, mpz_t p);

#endif
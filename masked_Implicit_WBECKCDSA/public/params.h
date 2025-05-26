#ifndef PARAMS_H
#define PARAMS_H

#include <gmp.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>

/* ECKCDSA 고정 파라미터 */
extern const char*  p_str,
                 *  n_str,
                 * Gx_str,
                 * Gy_str,
                 *  a_str,
                 *  b_str;

extern mpz_t p,
             n,
             Gx,
             Gy,
             a,
             b,
             e,  
             Qx,
             Qy;

/* 메세지 M */
extern const char* M_str;
extern const char* e_str;



extern const char* Qx_str,
                 * Qy_str;

extern const char* ki_G_str[128][4][2];

struct EQ_COEFF_TABLE {
    mpz_t coeff_m[127][4][12];
    mpz_t coeff_m_last[3][44];
};

extern const char* u_input[4];
extern const char* coeff_m_str[127][4][12];
extern const char* coeff_m_last_str[3][44];

extern const char* encoding_A_last[3][3];
extern const char* encoding_B_last[3];

extern const char* inv_6_str;
extern const char* inv_36_str;
extern mpz_t inv_6;
extern mpz_t inv_36;

#define ENC_ROUD 128

#endif
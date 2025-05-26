#ifndef SECRET_PARAMS
#define SECRET_PARAMS

#include <stdio.h>
#include <gmp.h>

//Test Vector를 맞추기 위해 생성
extern const char* k_str;
extern const char* k_values_str[128];
extern const char* E_str;

extern mpz_t k;


/*************(서명 변수 테이블)*************************/
extern const char* d_str,
                 * inv_d_str;
extern const char* d_values_str[128];
extern const char* k_values_str[128];

extern const char* k_i_str[128][4];
extern const char* alpha_i_str[128][4];
extern const char* beta_i_str[128][4];
extern const char* d_i_str[128][4];
#endif
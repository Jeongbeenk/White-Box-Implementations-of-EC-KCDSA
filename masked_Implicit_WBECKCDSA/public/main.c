#include "params.h"
#include "EC.h"
#include "ge.h"
#include "utils.h"
#include "sha256.h"

mpz_t p, n, Gx, Gy, a, b, e, M, Qx, Qy;
Point G;
mpz_t inv_6, inv_36;
mpz_t kGx,kGy;
mpz_t r;


void set_params()
{
    mpz_init_set_str( p,  p_str, 16); mpz_init_set_str( n,  n_str, 16); 
    mpz_init_set_str(Gx, Gx_str, 16); mpz_init_set_str(Gy, Gy_str, 16); 
    mpz_init_set_str( a,  a_str, 16); mpz_init_set_str( b,  b_str, 16);
    mpz_init_set_str(Qx, Qx_str, 16); mpz_init_set_str(Qy, Qy_str, 16);
    mpz_init_set_str(inv_6, inv_6_str, 16); mpz_init_set_str(inv_36, inv_36_str, 16);
    point_init_set_str(G, Gx_str, Gy_str, 16);
}

void clear_params()
{
    mpz_clear(p); mpz_clear(n); mpz_clear(Gx); mpz_clear(Gy); 
    mpz_clear(a); mpz_clear(b); 
    mpz_clear(Qx); mpz_clear(Qy);
    mpz_clear(inv_6); mpz_clear(inv_36);
    point_clear(G);
}

void mat_mul_3x3(mpz_t matrix[3][3], mpz_t vec[3], mpz_t res[3]) {
    mpz_t temp1, temp2;
    mpz_init(temp1); mpz_init(temp2);

    for( int i=0; i<3; i++ ) {
        mpz_set_ui(temp2, 0);

        for( int j=0; j<3; j++ ) {
            mpz_mul(temp1, matrix[i][j], vec[j]);
            mpz_add(temp2, temp2, temp1);
            mpz_mod(temp2, temp2, n);
        }
        mpz_set(res[i], temp2);
    }

    mpz_clear(temp1); mpz_clear(temp2);
}

void get_sha256_mpz(mpz_t input, mpz_t output) 
{
    size_t count;
    BYTE* byte_array = (BYTE*)mpz_export(NULL, &count, 1, sizeof(BYTE), 1, 0, input);

    SHA256_CTX ctx;
    BYTE hash_bytes[SHA256_BLOCK_SIZE];

    sha256_init(&ctx);

    sha256_update(&ctx, byte_array, count);

    sha256_final(&ctx, hash_bytes);

    mpz_import(output, SHA256_BLOCK_SIZE, 1, sizeof(BYTE), 1, 0, hash_bytes);

    free(byte_array);
}

void get_sha256_byte(BYTE* in_byte, size_t in_byte_length, mpz_t output)
{
    SHA256_CTX ctx;
    BYTE hash_bytes[SHA256_BLOCK_SIZE];

    sha256_init(&ctx);

    sha256_update(&ctx, in_byte, in_byte_length);

    sha256_final(&ctx, hash_bytes);

    mpz_import(output, SHA256_BLOCK_SIZE, 1, sizeof(BYTE), 1, 0, hash_bytes);
}

void Gen_Sign_r(mpz_t r)
{
    mpz_t ki_G_table[128][4][2]; init_table_3((mpz_t*)ki_G_table, 128, 4, 2);
    
    mpz_t Rx; mpz_init(Rx);

    for( int i=0; i<128; i++ ) {
        for( int j = 0; j<4; j++ ) {
            for( int k = 0; k<2; k++ ) {
                mpz_set_str(ki_G_table[i][j][k], ki_G_str[i][j][k], 16);
            }
        }
    }

    Point kG; point_init(kG);

    Point t; point_init(t);
    Point kiG; point_init(kiG);
    for( int i=0; i<128; i++ ) {
        int bit = get_two_bits(e, (i << 1));
        mpz_set(kiG->x, ki_G_table[i][bit][0]);
        mpz_set(kiG->y, ki_G_table[i][bit][1]);
        point_add(kG, t, kiG, a, p);
        point_set(t, kG);
    }

    gmp_printf("kG_x = %064ZX\n", kG->x);
    gmp_printf("kG_y = %064ZX\n", kG->y);

    mpz_set(kGx, kG->x);
    mpz_set(kGy, kG->y);

    /****************(예진 추가)*******************/
    mpz_set(Rx, kG->x);
    get_sha256_mpz(Rx, r);

    gmp_printf("r = %064ZX\n", r);
    /********************************************/

    point_clear(kG);
    point_clear(t);
    point_clear(kiG);
    mpz_clear(Rx);
    

}

void FUNC_COEFF(mpz_t rho)
{
    mpz_t u[4]; init_array((mpz_t*)u, 4);
    mpz_t v[4]; init_array((mpz_t*)v, 4);
    mpz_t res[3]; init_array((mpz_t*)res, 3);
    struct EQ_COEFF_TABLE eq_coeff_table; init_eq_coeff_table(&eq_coeff_table);

    for( int i=0; i<127; i++ ) {
        for( int j=0; j<4; j++ ) {
            for( int k=0; k<12; k++ ) {
                mpz_set_str(eq_coeff_table.coeff_m[i][j][k], coeff_m_str[i][j][k], 16);
            }
        }
    }
    for( int j=0; j<3; j++ ) {
        for( int k=0; k<44; k++ ) {
            mpz_set_str(eq_coeff_table.coeff_m_last[j][k], coeff_m_last_str[j][k], 16);
        }
    }

/*-------------(정빈 수정 : \secret\main (pre_A[0])^{-1}(-pre_b[0])를 연산한 값을 입력으로 받음 )--------------*/
    for(int i=0; i<4; i++){
        mpz_set_str(u[i], u_input[i], 16);
    }
/*----------------------------------------( 수정 끝 )-----------------------------------------------------*/

/*----------------------------------------( 예진 수정 시작 : rho 연결)-----------------------------------------------------*/
    
    // mpz_t rho; mpz_init_set_str(rho, "66C17858FE9968F5DEC19445C9B2DC4C21B255E74079329D62B2F928A308F4B7", 16);
    // mpz_t rho; mpz_init(rho);
    // size_t count_M;
    
    // // Qx와 Qy를 32 바이트로 고정하여 export
    // BYTE *byte_array_Qx = (BYTE *)malloc(32);
    // mpz_export(byte_array_Qx, NULL, 1, sizeof(BYTE), 1, 0, Qx);
    
    // BYTE *byte_array_Qy = (BYTE *)malloc(32);
    // mpz_export(byte_array_Qy, NULL, 1, sizeof(BYTE), 1, 0, Qy);
    
    // // M은 가변 길이로 export
    // BYTE *byte_array_M = (BYTE *)mpz_export(NULL, &count_M, 1, sizeof(BYTE), 1, 0, M);
    
    // // Qx||Qy||M 바이트 배열 결합
    // size_t total_length = 32 + 32 + count_M;
    // BYTE *concatenated_bytes = (BYTE *)malloc(total_length);
    
    // memcpy(concatenated_bytes, byte_array_Qx, 32);
    // memcpy(concatenated_bytes + 32, byte_array_Qy, 32);
    // memcpy(concatenated_bytes + 64, byte_array_M, count_M);

    // // 해시값 계산 - rho = H(Qx||Qy||M)
    // get_sha256_byte(concatenated_bytes, total_length, rho);


    // // 메모리 해제
    // free(byte_array_Qx);
    // free(byte_array_Qy);
    // free(byte_array_M);
    // free(concatenated_bytes);

    // // rho = r ^ H(Qx||Qy||M) mod n
    // mpz_xor(rho, r, rho); mpz_mod(rho, rho, n);

    // gmp_printf("rho = % 064ZX\n", rho);
/*----------------------------------------( 수정 끝 )-----------------------------------------------------*/

    encrypt(u, v, res, eq_coeff_table, e, rho, n);
//=======================================================================================================================================
    mpz_t tmp[3], tmp2[3]; init_array((mpz_t*)tmp, 3); init_array((mpz_t*)tmp2, 3);
    mpz_t A_last_inv[3][3]; init_table((mpz_t*)A_last_inv, 3, 3);
    mpz_t B_last[3]; init_array((mpz_t*)B_last, 3);
    mpz_t A_last[3][3];init_table((mpz_t*)A_last,3,3);
/*-----------------------( 정빈 수정 : \secret\main A_last[4][4], B_last 를 사용할 수 있도록 함수를 생성하고, 파라미터에 작성 )-------------*/
    for( int i=0; i<3; i++ ) {
        mpz_set_str(B_last[i],encoding_B_last[i], 16);
        for( int j=0; j<3; j++ ) {
            mpz_set_str(A_last[i][j],encoding_A_last[i][j], 16);
            
        }
    }
    mat_mul_3x3(A_last, res, tmp2);
    mpz_add(tmp[0], tmp2[0], B_last[0]); mpz_mod(tmp[0], tmp[0], n);
    mpz_add(tmp[1], tmp2[1], B_last[1]); mpz_mod(tmp[1], tmp[1], n);
    mpz_add(tmp[2], tmp2[2], B_last[2]); mpz_mod(tmp[2], tmp[2], n);

    gmp_printf("tmp[0] = %064ZX\n", tmp[0]);
    gmp_printf("tmp[1] = %064ZX\n", tmp[1]);
    gmp_printf("tmp[2] = %064ZX\n", tmp[2]);

/*----------------------------------------( 수정 끝 )-----------------------------------------------------*/
    mpz_t result2; mpz_init(result2);
    mpz_t s; mpz_init(s);
    mpz_sub(result2, tmp[0], rho); mpz_mod(result2,result2,n);
    mpz_mul(result2, result2, tmp[1]); mpz_mod(result2,result2,n);
    mpz_add(result2,result2,tmp[2]);mpz_mod(result2,result2,n);

    gmp_printf(" 서명 (r) : %064ZX\n",r);
    gmp_printf(" 서명 (s) : %064ZX\n",result2);
    
    mpz_set_str(r,"EC3847B0CA52038A823D023014546B414946EF0A6EE09228389484595F30E26C",16);
    mpz_set_str(s,"9B333457661C7CF741BDDBC0835553DFBB37EE74F53DB699E0A17780C7B6F1D0",16);
    if(mpz_cmp(kGx,r)!=0){
        printf("r = invalid!\n");
    }
    else{
        printf("r = valid!\n");
    }
  if(mpz_cmp(result2,s)!=0){
        printf("s = invalid!\n");
    }
    else{
        printf("s = valid!\n");
    }
    clear_array((mpz_t*)tmp, 3); clear_array((mpz_t*)tmp2, 3); clear_table((mpz_t*)A_last_inv, 3, 3); clear_array((mpz_t*)B_last, 3);

    clear_array((mpz_t*)u, 4); clear_array((mpz_t*)v, 4); clear_array((mpz_t*)res, 3);
    clear_eq_coeff_table(&eq_coeff_table);
    mpz_clear(result2);
    mpz_clear(s);

    /****(예진 추가)*****/
     mpz_clear(r); mpz_clear(rho);
    /****(예진 추가 끝)*****/
}

void Masking_Sign(mpz_t M)
{
    get_sha256_mpz(M,e);
    size_t count_M; 
    BYTE *byte_array_M = (BYTE *)mpz_export(NULL, &count_M, 1, sizeof(BYTE), 1, 0, M);
    gmp_printf(" 메시지 M의 해시비트 (e) : %064ZX \n",e);

    
    /*---------------------------( 필요한 부분이지만, 테스트벡터를 위해 주석처리  : e <- HASH{(int)M||(int)t} )-----------------------------------------------*/
    
    // gmp_randstate_t state; gmp_randinit_default(state); gmp_randseed_ui(state, 0);
    // mpz_t t; mpz_init(t);
    // mpz_urandomm(t, state, n);  
    // gmp_printf("t= %064ZX\n", t);
   
    // BYTE *byte_array_t = (BYTE *)malloc(32);
    // mpz_export(byte_array_t, NULL, 1, sizeof(BYTE), 1, 0, t);
    // size_t total_len_Mt = 32 + count_M;
    // BYTE *concatenated_bytes_Mt = (BYTE *)malloc(total_len_Mt);
    // memcpy(concatenated_bytes_Mt, byte_array_t,32);
    // memcpy(concatenated_bytes_Mt + 32, byte_array_M, count_M);
    // get_sha256_byte(concatenated_bytes_Mt, total_len_Mt, e);

    // gmp_printf(" 메시지 해시값 (e2) : %064ZX \n", e);
    // free(byte_array_t);
    // free(concatenated_bytes_Mt);
    // mpz_clear(t);
    /*------------------------------------------( 끝 : 테스트벡터를 위해 주석처리 )------------------------------------------------------*/

/*---------------------------( 1: 서명 r 생성부 )----------------------------------------------*/
    Gen_Sign_r(r);
/*---------------------------( 2: \rho 생성부 )----------------------------------------------*/
    mpz_t rho; mpz_init(rho);
    
    // Qx와 Qy를 32 바이트로 고정하여 export
    BYTE *byte_array_Qx = (BYTE *)malloc(32);
    mpz_export(byte_array_Qx, NULL, 1, sizeof(BYTE), 1, 0, Qx);
    
    BYTE *byte_array_Qy = (BYTE *)malloc(32);
    mpz_export(byte_array_Qy, NULL, 1, sizeof(BYTE), 1, 0, Qy);
    
    // Qx||Qy||M 바이트 배열 결합
    size_t total_length = 32 + 32 + count_M;
    BYTE *concatenated_bytes = (BYTE *)malloc(total_length);
    
    memcpy(concatenated_bytes, byte_array_Qx, 32);
    memcpy(concatenated_bytes + 32, byte_array_Qy, 32);
    memcpy(concatenated_bytes + 64, byte_array_M, count_M);

    // 해시값 계산 - rho = H(Qx||Qy||M)
    get_sha256_byte(concatenated_bytes, total_length, rho);
    mpz_xor(rho, r, rho); mpz_mod(rho, rho, n);
    gmp_printf("rho = % 064ZX\n", rho);

    // 메모리 해제
    free(byte_array_Qx);
    free(byte_array_Qy);
    free(byte_array_M);
    free(concatenated_bytes);
/*---------------------------( 3: 서명 s 생성부 )----------------------------------------------*/
    FUNC_COEFF(rho);
    mpz_clear(e); 
}

int main() 
{
    mpz_set_str( M,  "5468697320697320612073616D706C65206D65737361676520666F722045432D4B4344534120696D706C656D656E746174696F6E2076616C69646174696F6E2E", 16);
    
    set_params();
    Masking_Sign(M);
    clear_params();

    mpz_clear(M);
    return 0;
}


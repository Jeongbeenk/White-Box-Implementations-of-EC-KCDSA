#include "encode_gen.h"

mpz_t p, n, Gx, Gy, a, b, e;
Point G;
mpz_t k, d;

void set_params()
{
    mpz_set_str( p,  p_str, 16); mpz_set_str( n,  n_str, 16); 
    mpz_set_str(Gx, Gx_str, 16); mpz_set_str(Gy, Gy_str, 16); 
    mpz_set_str( a,  a_str, 16); mpz_set_str( b,  b_str, 16);
    mpz_set_str( e,  e_str, 16); mpz_set_str( k,  k_str, 16);
    mpz_set_str( d,  d_str, 16);

    point_init_set_str(G, Gx_str, Gy_str, 16);
}

void clear_params()
{
    mpz_clear(p); mpz_clear(n); mpz_clear(Gx); mpz_clear(Gy); mpz_clear(a); mpz_clear(b); mpz_clear(e); mpz_clear(k); mpz_clear(d);
    point_clear(G);
}

void mat_mul_4x4(mpz_t matrix[4][4], mpz_t vec[4], mpz_t res[4], mpz_t n) {
    mpz_t temp1, temp2;
    mpz_init(temp1);
    mpz_init(temp2);

    for (int i = 0; i < 4; i++) {
        mpz_set_ui(temp2, 0);  // temp2를 0으로 초기화하여 합산을 시작합니다.

        for (int j = 0; j < 4; j++) {
            mpz_mul(temp1, matrix[i][j], vec[j]);  // matrix[i][j] * vec[j]
            mpz_add(temp2, temp2, temp1);          // 결과를 temp2에 누적
        }

        mpz_mod(res[i], temp2, n);  // 각 행의 결과를 모듈로 n 연산하여 저장
    }

    mpz_clear(temp1);
    mpz_clear(temp2);
}

void Gen_Values(){
    gmp_randstate_t state; gmp_randinit_default(state);; gmp_randseed_ui(state, 0);
    unsigned int seed = 0;

    
    mpz_t d_i_table[128][4], alpha_i_table[128][4], beta_i_table[128][4], k_i_table[128][4];

    init_table((mpz_t*)alpha_i_table, 128, 4); init_table((mpz_t*)beta_i_table, 128, 4); init_table((mpz_t*)k_i_table, 128, 4); init_table((mpz_t*)d_i_table, 128, 4);

    make_and_write_ki_di_alpha_beta_table("secret_params.c", d_i_table, alpha_i_table, beta_i_table, k_i_table);

    // [ki]G 테이블 만들기
    mpz_t ki_G_table[128][4][2]; init_table_3((mpz_t*)ki_G_table, 128, 4, 2);
    cal_ki_G_with_table(k_i_str, G, e, ki_G_table);      // [ki]G 테이블 생성
    write_ki_G_table("../public/params.c", ki_G_table); // 생성한 테이블 파일에 쓰기


    // clear_table_3((mpz_t*)ki_G_table, 128, 4, 2);
    clear_table((mpz_t*)k_i_table, 128, 4);
    clear_table((mpz_t*)d_i_table, 128, 4);
    clear_table((mpz_t*)alpha_i_table, 128, 4);
    clear_table((mpz_t*)beta_i_table, 128, 4);
}

void Gen_Eq_Coeff_Table()
{   
    Gen_Values();
    mpz_t encoding_M[127][4][4], encoding_A[127][4][4], encoding_B[127][4];
    mpz_t pre_encoding_A[128][4][4], pre_encoding_B[128][4];
    mpz_t encoding_M_last[3][3], encoding_A_last[3][3], encoding_B_last[3];

    for (int i=0;i<127;i++){
        init_table((mpz_t*)encoding_A[i], 4, 4);
        init_table((mpz_t*)encoding_M[i], 4, 4);
        init_array((mpz_t*)encoding_B[i], 4);
    }
    for( int i = 0; i < 128; i++ ) {
        init_table((mpz_t*)pre_encoding_A[i], 4, 4);
        init_array((mpz_t*)pre_encoding_B[i], 4);
    }

    for (int i=0;i<127;i++){
        init_identity_matrix_4x4(encoding_M[i]);    
        init_identity_matrix_4x4(encoding_A[i]);  
    }
    
    for( int i=0 ; i<128; i++ ) {
        init_identity_matrix_4x4(pre_encoding_A[i]); 
    }

    init_identity_matrix_3x3(encoding_M_last);
    init_identity_matrix_3x3(encoding_A_last);
    init_array((mpz_t*)encoding_B_last, 3);

/*****************************************************( 1. 인코딩 관련 : 시작 )*********************************************************************/

    gmp_randstate_t state;
    gmp_randinit_default(state);
    gmp_randseed_ui(state, 0);

    mpz_t tmp,tmp2,tmp3,tmp4; mpz_init(tmp),mpz_init(tmp2),mpz_init(tmp3),mpz_init(tmp4);


/*----------------------( 랜덤 인코딩 생성 )-------------------------------------*/
    // 1~127 라운드 인코딩 : 4x4, 4x1
    for( int i=0; i<127; i++ ) {
        gen_random_invertible_matrix_4x4(encoding_A[i], state, n);
    }
    for(int i=0;i<127;i++){
        gen_random_invertible_matrix_4x4(encoding_M[i], state, n);

    }
    for( int i=0; i<127; i++ ) gen_random_vector_4(encoding_B[i], state, n);

    //128 라운드 인코딩 : 3x3, 3x1
    gen_random_invertible_matrix_3x3(encoding_M_last, state, n);
    gen_random_invertible_matrix_3x3(encoding_A_last, state, n);
    gen_random_vector_3(encoding_B_last, state, n);

/*-----------( 인코딩 상쇄 구조 : 이전 출력 인코딩 = 다음 입력 인코딩 )------------------*/
    for( int i=0; i<127; i++ ) {
        for( int j=0; j<4; j++ ) {
            for( int k=0; k<4; k++ ) {
                mpz_set(pre_encoding_A[i+1][j][k], encoding_A[i][j][k]);
            }
        }
    }
    for( int i=0; i<127; i++ ) {
        for( int j=0; j<4; j++ ) {
            mpz_set(pre_encoding_B[i+1][j], encoding_B[i][j]);
        }
    }

    gen_random_vector_4(pre_encoding_B[0],state,n);
    gen_random_invertible_matrix_4x4(pre_encoding_A[0],state,n); 


/*-----------( public u input 쓰기 :\secret\main (pre_A[0])^{-1}(-pre_b[0])를 연산한 값을 입력으로 받음)---------------*/

    mpz_t temp[4]; init_array(temp,4);
    mpz_t u[4];init_array((mpz_t *)u,4);
    mpz_t first_pre_A_inv[4][4]; init_table((mpz_t *)first_pre_A_inv,4,4);

    mpz_neg(tmp,pre_encoding_B[0][0]); mpz_mod(tmp,tmp,n);
    // gmp_printf("%064ZX\n",tmp);
    mpz_neg(tmp2,pre_encoding_B[0][1]); mpz_mod(tmp2,tmp2,n);
    // gmp_printf("%064ZX\n",tmp2);
    mpz_neg(tmp3,pre_encoding_B[0][2]); mpz_mod(tmp3,tmp3,n);
    // gmp_printf("%064ZX\n",tmp3);
    mpz_neg(tmp4,pre_encoding_B[0][3]); mpz_mod(tmp4,tmp4,n);
    // gmp_printf("%064ZX\n",tmp4);

    mpz_set(temp[0],tmp);
    mpz_set(temp[1],tmp2);
    mpz_set(temp[2],tmp3);
    mpz_set(temp[3],tmp4);

    get_inverse_matrix_4x4(first_pre_A_inv,pre_encoding_A[0],n);
    mat_mul_4x4(first_pre_A_inv,temp,u,n);
    write_u_input("../public/params.c",&u);

/*--------------------------( public last encoding 쓰기 )---------------------------------------------------------*/

    write_last_encoding_A("../public/params.c",&encoding_A_last);
    write_last_encoding_B("../public/params.c",&encoding_B_last);

    clear_table(first_pre_A_inv,4,4);
    clear_array(temp,4); clear_array(u,4);
    mpz_clear(tmp),mpz_clear(tmp2),mpz_clear(tmp3),mpz_clear(tmp4);
/*****************************************************( 인코딩 관련 : 끝 )*************************************************************************/

/**********************************************( 2. 변수 쓰기 및 방정식 계수 배열 쓰기 )***************************************************************/

    // ki_table 가져오기
    mpz_t k_i_table[128][4]; init_table((mpz_t*)k_i_table, 128, 4);
    for( int i=0; i<128; i++ ) {
        for( int j=0; j<4; j++ ) {
            mpz_set_str(k_i_table[i][j], k_i_str[i][j], 16);
        }
    }
    
    // alpha_i_table 가져오기
    mpz_t alpha_i_table[128][4]; init_table((mpz_t*)alpha_i_table, 128, 4);
    for( int i=0; i<128; i++ ) {
        for( int j=0; j<4; j++ ) {
            mpz_set_str(alpha_i_table[i][j], alpha_i_str[i][j], 16);
        }
    }

    // beta_i_table 가져오기
    mpz_t beta_i_table[128][4]; init_table((mpz_t*)beta_i_table, 128, 4);
    for( int i=0; i<128; i++ ) {
        for( int j=0; j<4; j++ ) {
            mpz_set_str(beta_i_table[i][j], beta_i_str[i][j], 16);
        }
    }

    // d_values 가져오기
    mpz_t d_values[128]; init_array((mpz_t*)d_values, 128);
    for( int i=0; i<128; i++ ) {
        mpz_set_str(d_values[i], d_values_str[i], 16);
    }

    // 방정식 계수 배열 선언
    struct EQ_COEFF_TABLE eq_coeff_table;
    init_eq_coeff_table(&eq_coeff_table);

    // 방정식 계수 배열 생성
    gen_eq_coeff_table(encoding_M, encoding_A, encoding_B, pre_encoding_A, pre_encoding_B, 
                        encoding_M_last, encoding_A_last, encoding_B_last, 
                        k_i_table, alpha_i_table, beta_i_table, d_values, &eq_coeff_table, n);

    write_eq_coeff_table("../public/params.c", &eq_coeff_table);

    

    /**************** 인코딩 행렬 및 테이블 메모리 해제 ****************/
    clear_eq_coeff_table(&eq_coeff_table);
    clear_table((mpz_t*)k_i_table, 128, 4);
    clear_table((mpz_t*)alpha_i_table, 128, 4);
    clear_table((mpz_t*)beta_i_table, 128, 4);
    clear_array((mpz_t*)d_values, 128);

    for( int i = 0; i < 127; i++ ) {
        clear_table((mpz_t*)encoding_M[i], 4, 4);
        clear_table((mpz_t*)encoding_A[i], 4, 4);
        clear_array((mpz_t*)encoding_B[i], 4);
    }
    for( int i = 0; i < 128; i++ ) {
        clear_table((mpz_t*)pre_encoding_A[i], 4, 4);
        clear_array((mpz_t*)pre_encoding_B[i], 4);
    }
    
}

int main()
{
    set_params();
    Gen_Eq_Coeff_Table();
    clear_params();
    return 0;
}






/*------(테스트 벡터 확인을 위해서 고정된 k와 d를 d[128],k[128]을 생성하는 코드)-------*/
/*
void FUNC_divide_values()
{   
    gmp_randstate_t state; gmp_randinit_default(state);gmp_randseed_ui(state, 0);
    unsigned int seed = 0;
    
//==========================================================
// 변수 및 테이블 선언
//==========================================================
    // k_values
    mpz_t k_values[4][128]; init_table((mpz_t*)k_values, 4, 128);

    // d_values
    mpz_t d_values[128]; init_array((mpz_t*)d_values, 128);

    // alpha_values
    mpz_t alpha_values[128]; init_array((mpz_t*)alpha_values, 128);
    mpz_t alpha; mpz_init(alpha);
    mpz_urandomm(alpha, state, n);
    
    // beta_values
    mpz_t beta_values[128]; init_array((mpz_t*)beta_values, 128);
    mpz_t beta; mpz_init(beta);
    mpz_urandomm(beta, state, n);


//==========================================================
// k, d, alpha, beta 나누기
//==========================================================
    do{
        for( int i = 0; i < 4; i++ ) 
            divide_values(k, k_values[i], seed++);
    }while( cal_max_of_ki(k_values) );

    divide_values(d, d_values, seed++);

    divide_values(alpha, alpha_values, seed++);

    divide_values(beta, beta_values, seed++);


//==========================================================
// 테이블 및 파일에 쓰기
//==========================================================
    // k_values, d_values 파일에 쓰기
    write_values("secret_params.c", k_values, d_values);

    // alpha, beta values 파일에 쓰기
    mpz_t alpha_i_table[128][4], beta_i_table[128][4], ki_table[128][4];
    init_table((mpz_t*)alpha_i_table, 128, 4); init_table((mpz_t*)beta_i_table, 128, 4); init_table((mpz_t*)ki_table, 128, 4);

    make_and_write_ki_alpha_beta_table("secret_params.c", k, alpha, beta, alpha_values, beta_values, k_values, alpha_i_table, beta_i_table, ki_table);

    // [ki]G 테이블 만들기
    mpz_t ki_G_table[128][4][2]; init_table_3((mpz_t*)ki_G_table, 128, 4, 2);
    cal_ki_G_with_table(ki_str, G, e, ki_G_table);      // [ki]G 테이블 생성
    write_ki_G_table("../public/params.c", ki_G_table); // 생성한 테이블 파일에 쓰기
    
    // Point kG; point_init(kG);

    // Point t; point_init(t);
    // Point kiG; point_init(kiG);

    // // puts("\n디버깅 kiG");
    // for( int i=0; i<128; i++ ) {
    //     int bit = get_two_bits(e, (i << 1));
    //     mpz_set(kiG->x, ki_G_table[i][bit][0]);
    //     mpz_set(kiG->y, ki_G_table[i][bit][1]);
    //     point_add(kG, t, kiG, a, p);
    //     point_set(t, kG);
    // }

    // gmp_printf("\nkG_x = %064ZX\n", kG->x);
    // gmp_printf("kG_y = %064ZX\n", kG->y);

    // // Point 메모리 해제
    // point_clear(kG);
    // point_clear(t);
    // point_clear(kiG);


//==========================================================
// 변수 및 테이블 메모리 해제
//==========================================================
    clear_table_3((mpz_t*)ki_G_table, 128, 4, 2);
    clear_table((mpz_t*)k_values, 4, 128);
    clear_array((mpz_t*)d_values, 128);
    clear_array((mpz_t*)alpha_values, 128);
    clear_array((mpz_t*)beta_values, 128);
    clear_table((mpz_t*)ki_table, 128, 4);
    clear_table((mpz_t*)alpha_i_table, 128, 4);
    clear_table((mpz_t*)beta_i_table, 128, 4);

    mpz_clear(alpha); mpz_clear(beta);
}
*/


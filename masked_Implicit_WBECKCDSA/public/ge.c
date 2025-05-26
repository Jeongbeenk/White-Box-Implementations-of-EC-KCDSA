#include "ge.h"


void gauss_elimination_3x3_mod(mpz_t Mat[3][4], mpz_t y[3], const mpz_t p) {
    mpz_t ratio, temp;
    mpz_inits(ratio, temp, NULL);

    // Step 1: Forward elimination
    for (int i = 0; i < 3; i++) {
        // Check for division by zero
        if (mpz_cmp_ui(Mat[i][i], 0) == 0) {
            printf("Error: Division by zero in row %d.\n", i);
            return;
        }

        // Normalize the pivot row
        mpz_invert(temp, Mat[i][i], p); // temp = 1 / Mat[i][i] (mod p)
        for (int j = 0; j < 4; j++) {
            mpz_mul(Mat[i][j], Mat[i][j], temp); // Mat[i][j] *= temp
            mpz_mod(Mat[i][j], Mat[i][j], p); // Ensure we take modulo p
        }

        // Eliminate the current variable from the rows below
        for (int j = i + 1; j < 3; j++) {
            mpz_set(ratio, Mat[j][i]); // ratio = Mat[j][i]

            for (int k = i; k < 4; k++) {
                // Mat[j][k] -= ratio * Mat[i][k] (mod p)
                mpz_mul(temp, ratio, Mat[i][k]); // temp = ratio * Mat[i][k]
                mpz_sub(Mat[j][k], Mat[j][k], temp); // Mat[j][k] -= temp
                mpz_mod(Mat[j][k], Mat[j][k], p); // Ensure we take modulo p
            }
        }
    }

    // Step 2: Backward substitution
    for (int i = 2; i >= 0; i--) {
        if (mpz_cmp_ui(Mat[i][i], 0) == 0) {
            printf("Error: Division by zero in row %d.\n", i);
            return;
        }

        // Solve for y[i]
        mpz_invert(temp, Mat[i][i], p); // temp = 1 / Mat[i][i] (mod p)
        mpz_mul(y[i], temp, Mat[i][3]); // y[i] = (Mat[i][3] / Mat[i][i]) % p
        mpz_mod(y[i], y[i], p); // Ensure we take modulo p

        // Back substitute to eliminate the value from the previous rows
        for (int j = 0; j < i; j++) {
            // Mat[j][3] -= Mat[j][i] * y[i] (mod p)
            mpz_mul(temp, Mat[j][i], y[i]); // temp = Mat[j][i] * y[i]
            mpz_sub(Mat[j][3], Mat[j][3], temp); // Mat[j][3] -= temp
            mpz_mod(Mat[j][3], Mat[j][3], p); // Ensure we take modulo p
        }
    }

    // Clear variables
    mpz_clears(ratio, temp, NULL);
}

void gauss_elimination_4x4_mod(mpz_t Mat[4][5], mpz_t y[4], const mpz_t p) {
    mpz_t ratio, temp;
    mpz_inits(ratio, temp, NULL);

/*     Ay+B=0 [A|-B]  */
    // for( int i = 0; i < 4; i++ ) {
    //     mpz_neg(Mat[i][4], Mat[i][4]);
    //     mpz_mod(Mat[i][4], Mat[i][4], p);
    //     // gmp_printf("Mat[%d][4] = %Zd\n", i, Mat[i][4]);
    // }

    // Step 1: Forward elimination
    for (int i = 0; i < 4; i++) {


        // Check for division by zero
        if (mpz_cmp_ui(Mat[i][i], 0) == 0) {
            printf("Error: Division by zero in row %d.\n", i);
            return;
        }

        // Normalize the pivot row by finding the inverse of Mat[i][i]
        mpz_invert(temp, Mat[i][i], p); // temp = 1 / Mat[i][i] (mod p)
        for (int j = 0; j < 5; j++) {
            mpz_mul(Mat[i][j], Mat[i][j], temp); // Mat[i][j] *= temp
            mpz_mod(Mat[i][j], Mat[i][j], p); // Ensure we take modulo p
        }

        // Eliminate the current variable from the rows below
        for (int j = i + 1; j < 4; j++) {
            mpz_set(ratio, Mat[j][i]); // ratio = Mat[j][i]

            for (int k = i; k < 5; k++) {
                // Mat[j][k] -= ratio * Mat[i][k] (mod p)
                mpz_mul(temp, ratio, Mat[i][k]); // temp = ratio * Mat[i][k]
                mpz_sub(Mat[j][k], Mat[j][k], temp); // Mat[j][k] -= temp
                mpz_mod(Mat[j][k], Mat[j][k], p); // Ensure we take modulo p
            }
        }
    }

    // Step 2: Backward substitution
    for (int i = 3; i >= 0; i--) {
        if (mpz_cmp_ui(Mat[i][i], 0) == 0) {
            printf("Error: Division by zero in row %d.\n", i);
            return;
        }

        // Solve for y[i]
        mpz_invert(temp, Mat[i][i], p); // temp = 1 / Mat[i][i] (mod p)
        mpz_mul(y[i], temp, Mat[i][4]); // y[i] = (Mat[i][4] / Mat[i][i]) % p
        mpz_mod(y[i], y[i], p); // Ensure we take modulo p

        // Back substitute to eliminate the value from the previous rows
        for (int j = 0; j < i; j++) {
            // Mat[j][4] -= Mat[j][i] * y[i] (mod p)
            mpz_mul(temp, Mat[j][i], y[i]); // temp = Mat[j][i] * y[i]
            mpz_sub(Mat[j][4], Mat[j][4], temp); // Mat[j][4] -= temp
            mpz_mod(Mat[j][4], Mat[j][4], p); // Ensure we take modulo p
        }
    }

    // Clear variables
    mpz_clears(ratio, temp, NULL);
}


/**************************************************(방정식 계수 배열 업데이트)*************************************************************************************************/
/*----------------( 1~127 라운드 )----------------*/
void matrix_update_from_masking(mpz_t ge_matrix[4][5], mpz_t u[4], mpz_t coeff_m_table[4][12], mpz_t hash_e, mpz_t n, int round)
{
    int i, j, ptr;

    mpz_t constant; mpz_init(constant);

    int bit = get_two_bits(hash_e, (round << 1));
    mpz_t ei[3];    init_array((mpz_t*)ei, 3);
    mpz_set_ui(ei[0], bit);
    mpz_set_ui(ei[1], bit*bit);
    mpz_set_ui(ei[2], bit*bit*bit);

    // gmp_printf("ei = %ZX\n", ei[0]);

    for( int it = 0; it < 4; it++ )
    {
        mpz_t temp1, temp2;   mpz_init(temp1);  mpz_init(temp2);
        mpz_set_ui(constant, 0);

        // 1) vj : 4개
        mpz_set(ge_matrix[it][0], coeff_m_table[it][0]);
        mpz_set(ge_matrix[it][1], coeff_m_table[it][1]);
        mpz_set(ge_matrix[it][2], coeff_m_table[it][2]);
        mpz_set(ge_matrix[it][3], coeff_m_table[it][3]);

        // 2) uj : 4개
        for( i = 0; i < 4; i++ ) {
            mpz_mul(temp1, u[i], coeff_m_table[it][4 + i]);         mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                     mpz_mod(constant, constant, n);
        }
        // gmp_printf("[%d] constant = %064ZX\n", it, constant);

        // 3~5) ei^3, ei^2, ei : 각 1개 -> 더한 후 6으로 나누어서 저장
        for( i = 0; i < 3; i++ ) {
            mpz_mul(temp1, coeff_m_table[it][8 + i], ei[2 - i]);    mpz_mod(temp1, temp1, n);
            mpz_add(temp2, temp2, temp1);                           mpz_mod(temp2, temp2, n);

            // gmp_printf("coeff[%d][%d] = %064ZX\n", it, 8+i, coeff_m_table[it][8 + i]);
            // gmp_printf("temp1 = %064ZX\n", temp1);
            // gmp_printf("temp2 = %064ZX\n", temp2);
        }
        mpz_mul(temp2, temp2, inv_6);                               mpz_mod(temp2, temp2, n);
        mpz_add(constant, constant, temp2);                         mpz_mod(constant, constant, n);
        // gmp_printf("[%d] constant = %064ZX\n", it, constant);

        // 6) constant : 1개
        mpz_add(constant, constant, coeff_m_table[it][11]);         mpz_mod(constant, constant, n);
        // gmp_printf("[%d] constant = %064ZX\n\n", it, constant);

        mpz_set(ge_matrix[it][4], constant);

        mpz_clear(temp1);   mpz_clear(temp2);
    }
    
    mpz_clear(constant);
    clear_array((mpz_t*)ei, 3);
}

/*----------------( 128 라운드 )----------------*/
void matrix_update_from_masking_128(mpz_t ge_matrix[3][4], mpz_t u[4], mpz_t coeff_m_table[3][44], mpz_t hash_e, mpz_t rho, mpz_t n)
{
    int i, j, ptr;

    mpz_t constant; mpz_init(constant);

    int bit = get_two_bits(hash_e, (127 << 1));
    mpz_t ei[3];    init_array((mpz_t*)ei, 3);
    mpz_set_ui(ei[0], bit);
    mpz_set_ui(ei[1], bit*bit);
    mpz_set_ui(ei[2], bit*bit*bit);

    // puts("입력 정보");
    // for(int i=0; i<4; i++) gmp_printf("u[%d] = %064ZX\n", i, u[i]);

    for( int it = 0; it < 3; it++ )
    {
        mpz_t temp1, temp2;   mpz_init(temp1);  mpz_init(temp2);
        mpz_set_ui(constant, 0);

        // 1) vj : 3개
        mpz_set(ge_matrix[it][0], coeff_m_table[it][0]);
        mpz_set(ge_matrix[it][1], coeff_m_table[it][1]);
        mpz_set(ge_matrix[it][2], coeff_m_table[it][2]);

        // 2) (uj)^2 : 4개
        for( i = 0; i < 4; i++ ) {
            mpz_pow_ui(temp1, u[i], 2);                         mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, coeff_m_table[it][3 + i]);    mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                 mpz_mod(constant, constant, n);
            // gmp_printf("uj^2 constant = %064ZX\n", constant);
        }

        // gmp_printf("constant = %064ZX\n", constant);
        
        // 3) ujuk : 6개
        ptr = 7;
        for( i = 0; i < 4; i++ ) {
            for( j = i+1; j < 4; j++ ) {
                mpz_mul(temp1, u[i], u[j]);                     mpz_mod(temp1, temp1, n);
                mpz_mul(temp1, temp1, coeff_m_table[it][ptr]);  mpz_mod(temp1, temp1, n);
                mpz_add(constant, constant, temp1);             mpz_mod(constant, constant, n);
                ptr++;
            }
        }
        // gmp_printf("constant = %064ZX\n", constant);

        // 4~6) (e_128)^3uj, (e_128)^2uj, e_128uj : 각 4개 -> 모두 계산해서 더한 후, 6으로 나눔
        for( i = 0; i < 4; i++ ) {
            mpz_mul(temp1, u[i], coeff_m_table[it][13 + i]);    mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, ei[2]);                       mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, inv_6);                       mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                 mpz_mod(constant, constant, n);
        }
        for( i = 0; i < 4; i++ ) {
            mpz_mul(temp1, u[i], coeff_m_table[it][17 + i]);    mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, ei[1]);                       mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, inv_6);                       mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                 mpz_mod(temp1, temp1, n);
        }
        for( i = 0; i < 4; i++ ) {
            mpz_mul(temp1, u[i], coeff_m_table[it][21 + i]);    mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, ei[0]);                       mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, inv_6);                       mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                 mpz_mod(temp1, temp1, n);
        }

        // gmp_printf("constant = %064ZX\n", constant);

        // 7) uj : 4개
        for( i = 0; i < 4; i++ ) {
            mpz_mul(temp1, u[i], coeff_m_table[it][25 + i]);        mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                     mpz_mod(constant, constant, n);
        }
        // gmp_printf("constant = %064ZX\n", constant);

        /***********************************(예진 수정)********************************************/
        // 7-1) rho * uj : 4개
        for( i = 0; i < 4; i++ ) {
            mpz_mul(temp1, u[i], coeff_m_table[it][29 + i]);        mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, rho);                             mpz_mod(temp1, temp1, n);
            mpz_add(constant, constant, temp1);                     mpz_mod(constant, constant, n);
        }
        /***********************************(예진 수정 끝)******************************************/

        // 8) (e_128)^6 : 1개
        mpz_mul(temp1, ei[2], ei[2]);                               mpz_mod(temp1, temp1, n); // (e_128)^6
        mpz_mul(temp1, temp1, coeff_m_table[it][33]);               mpz_mod(temp1, temp1, n);

        // 9) (e_128)^5 : 1개
        mpz_mul(temp2, ei[2], ei[1]);                               mpz_mod(temp2, temp2, n); // (e_128)^5
        mpz_mul(temp2, temp2, coeff_m_table[it][34]);               mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2

        // 10) (e_128)^4 : 1개
        mpz_mul(temp2, ei[1], ei[1]);                               mpz_mod(temp2, temp2, n); // (e_128)^4
        mpz_mul(temp2, temp2, coeff_m_table[it][35]);               mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2
        
        // 11) (e_128)^3 : 1개
        mpz_mul(temp2, ei[2], coeff_m_table[it][36]);               mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2

        /***********************************(예진 수정)********************************************/
        // 11-1) e_128^3 * rho : 1개
        mpz_mul(temp2, ei[2], coeff_m_table[it][37]);               mpz_mod(temp2, temp2, n);
        mpz_mul(temp2, temp2, rho);                                 mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2
        /***********************************(예진 수정 끝)******************************************/

        // 12) (e_128)^2 : 1개
        mpz_mul(temp2, ei[1], coeff_m_table[it][38]);               mpz_mod(temp2, temp2, n);
        
        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2

        /***********************************(예진 수정)********************************************/
        // 12-1) e_128^2 * rho : 1개
        mpz_mul(temp2, ei[1], coeff_m_table[it][39]);               mpz_mod(temp2, temp2, n);
        mpz_mul(temp2, temp2, rho);                                 mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2
        /***********************************(예진 수정 끝)******************************************/

        // 13) e_128 : 1개
        mpz_mul(temp2, ei[0], coeff_m_table[it][40]);               mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2

        /***********************************(예진 수정)********************************************/
        // 13-1) e_128 * rho : 1개
        mpz_mul(temp2, ei[0], coeff_m_table[it][41]);               mpz_mod(temp2, temp2, n);
        mpz_mul(temp2, temp2, rho);                                 mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n); // temp1 + temp2

        mpz_mul(temp1, temp1, inv_36);                              mpz_mod(temp1, temp1, n);
        mpz_add(constant, constant, temp1);                         mpz_mod(constant, constant, n);

        // 13-2) rho : 1개
        mpz_mul(temp1, rho, coeff_m_table[it][42]);                 mpz_mod(temp1, temp1, n);
        mpz_add(constant, constant, temp1);                         mpz_mod(temp1, temp1, n);
        /***********************************(예진 수정 끝)******************************************/

        // 14) constant : 1개
        mpz_add(constant, constant, coeff_m_table[it][43]);         mpz_mod(constant, constant, n);

        // gmp_printf("constant = %064ZX\n", constant);
        mpz_set(ge_matrix[it][3], constant);
        
        mpz_clear(temp1);   mpz_clear(temp2);
    }
    
    mpz_clear(constant);
}



/**************************************************(라운드 연산)*************************************************************************************************/
void encrypt(mpz_t u[4], mpz_t v[4], mpz_t res[3], struct EQ_COEFF_TABLE eq_coeff_table, mpz_t hash_e, mpz_t rho, mpz_t n)
{    
    mpz_t ge_matrix[4][5]; init_table((mpz_t*)ge_matrix, 4, 5);
    mpz_t ge_matrix_last[3][4]; init_table((mpz_t*)ge_matrix_last, 3, 4);
    int round = 0;
    for( round = 0; round < 127; round++ ) {
        matrix_update_from_masking(ge_matrix, u, eq_coeff_table.coeff_m[round], hash_e, n, round);
        gauss_elimination_4x4_mod(ge_matrix, v, n);
        /* 이전 라운드 출력을 다음 라운드 입력으로 넣어주기 */
        for( int i = 0; i < 4; i++ ) {
            mpz_set(u[i], v[i]);
        }
    }    
    matrix_update_from_masking_128(ge_matrix_last, u, eq_coeff_table.coeff_m_last, hash_e, rho, n);
    gauss_elimination_3x3_mod(ge_matrix_last, res, n);
    
    
    clear_table((mpz_t*)ge_matrix, 4, 5);   clear_table((mpz_t*)ge_matrix_last, 3, 4);
}



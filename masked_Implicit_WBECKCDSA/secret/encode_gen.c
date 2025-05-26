#include "encode_gen.h"
#include "secret_params.h"

void init_identity_matrix_3x3(mpz_t Matrix[3][3])
{
    for( int i = 0; i < 3; i++ ) {
        for( int j = 0; j < 3; j++ ) {
            if( i == j ) mpz_set_ui(Matrix[i][j], 1);
            else         mpz_set_ui(Matrix[i][j], 0);
        }
    }
}

void init_identity_matrix_4x4(mpz_t Matrix[4][4])
{
    for( int i = 0; i < 4; i++ ) {
        for( int j = 0; j < 4; j++ ) {
            if( i == j ) mpz_set_ui(Matrix[i][j], 1);
            else         mpz_set_ui(Matrix[i][j], 0);
        }
    }
}

void generate_random_matrix_3x3(mpz_t matrix[3][3], gmp_randstate_t state, mpz_t mod) {
    for( int i = 0; i < 3; i++ ) {
        for( int j = 0; j < 3; j++ ) {
            mpz_urandomm(matrix[i][j], state, mod);
        }
    }
}

void generate_random_matrix_4x4(mpz_t matrix[4][4], gmp_randstate_t state, mpz_t mod) {
    for (int i = 0; i < SIZE; i++) {
        for (int j = 0; j < SIZE; j++) {
            mpz_urandomm(matrix[i][j], state, mod);  // mod n으로 랜덤값 생성
        }
    }
}

void cal_determinant_4x4(mpz_t det, mpz_t matrix[4][4], int n, mpz_t mod) {
    if (n == 1) {
        mpz_mod(det, matrix[0][0], mod);  // mod n으로 계산
        return;
    }

    // 디버깅용 출력: 현재 재귀 깊이 및 행렬 크기
    // printf("Calculating determinant for size: %d\n", n); fflush(stdout);

    mpz_t sub_det, sign, term;
    mpz_init(sub_det);  // 초기화
    mpz_init_set_si(sign, 1);  // 행렬식의 기호
    mpz_init(term);  // 초기화

    // 고정된 크기의 임시 배열을 사용
    mpz_t temp[4][4];

    // temp 배열 초기화 (mpz_init)
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mpz_init(temp[i][j]);
        }
    }

    mpz_set_ui(det, 0);  // 행렬식을 계산하기 전에 초기화

    for (int i = 0; i < n; i++) {
        // 소행렬 생성
        int sub_i = 0;
        for (int row = 1; row < n; row++) {
            int sub_j = 0;
            for (int col = 0; col < n; col++) {
                if (col == i) continue;  // i번째 열 제외
                mpz_set(temp[sub_i][sub_j], matrix[row][col]);  // 값 복사
                sub_j++;
            }
            sub_i++;
        }

        // 소행렬식 재귀 호출
        cal_determinant_4x4(sub_det, temp, n - 1, mod);

        // 현재 항을 더하거나 뺌
        mpz_mul(term, sign, sub_det);
        mpz_mul(term, term, matrix[0][i]);
        mpz_mod(term, term, mod);  // mod 연산

        mpz_add(det, det, term);
        mpz_mod(det, det, mod);  // mod 연산

        // 기호 변경
        mpz_neg(sign, sign);
    }

    // temp 배열 해제 (mpz_clear)
    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mpz_clear(temp[i][j]);
        }
    }

    mpz_clear(sub_det);
    mpz_clear(sign);
    mpz_clear(term);

    // 디버깅용 출력: 계산 완료
    // printf("Finished calculating determinant for size: %d\n", n); fflush(stdout);
}

void cal_determinant_3x3(mpz_t det, mpz_t matrix[3][3], mpz_t mod) {
    mpz_t term1, term2, term3;
    mpz_inits(term1, term2, term3, NULL);

    // term1 = a * (e * i - f * h) % mod
    mpz_mul(term1, matrix[1][1], matrix[2][2]); // e * i
    mpz_submul(term1, matrix[1][2], matrix[2][1]); // e * i - f * h
    mpz_mul(term1, matrix[0][0], term1); // a * (e * i - f * h)
    mpz_mod(term1, term1, mod);

    // term2 = -b * (d * i - f * g) % mod
    mpz_mul(term2, matrix[1][0], matrix[2][2]); // d * i
    mpz_submul(term2, matrix[1][2], matrix[2][0]); // d * i - f * g
    mpz_mul(term2, matrix[0][1], term2); // b * (d * i - f * g)
    mpz_neg(term2, term2); // -b * (d * i - f * g)
    mpz_mod(term2, term2, mod);

    // term3 = c * (d * h - e * g) % mod
    mpz_mul(term3, matrix[1][0], matrix[2][1]); // d * h
    mpz_submul(term3, matrix[1][1], matrix[2][0]); // d * h - e * g
    mpz_mul(term3, matrix[0][2], term3); // c * (d * h - e * g)
    mpz_mod(term3, term3, mod);

    // det = (term1 + term2 + term3) % mod
    mpz_add(det, term1, term2);
    mpz_add(det, det, term3);
    mpz_mod(det, det, mod);

    mpz_clears(term1, term2, term3, NULL);
}

void gen_random_invertible_matrix_3x3(mpz_t matrix[3][3], gmp_randstate_t state, mpz_t mod) {
    mpz_t det;
    mpz_init(det);
    int is_invertible = 0;

    while( !is_invertible ) {
        // 랜덤 행렬 생성 (mod n)
        generate_random_matrix_3x3(matrix, state, mod);

        // 행렬식 계산 (mod n)
        mpz_set_ui(det, 0);
        cal_determinant_3x3(det, matrix, mod);

        // 행렬식이 0이 아니면 역행렬이 존재함
        if( mpz_sgn(det) != 0 ) {
            is_invertible = 1;
        }
    }

    mpz_clear(det);
}

void gen_random_invertible_matrix_4x4(mpz_t matrix[4][4], gmp_randstate_t state, mpz_t mod) {
    mpz_t det;
    mpz_init(det);
    int is_invertible = 0;

    while (!is_invertible) {
        // 랜덤 행렬 생성 (mod n)
        generate_random_matrix_4x4(matrix, state, mod);

        // 행렬식 계산 (mod n)
        mpz_set_ui(det, 0);  // 행렬식을 계산하기 전에 초기화
        cal_determinant_4x4(det, matrix, 4, mod);

        // 행렬식이 0이 아니면 역행렬이 존재함
        if (mpz_sgn(det) != 0) {
            is_invertible = 1;
        }
    }

    mpz_clear(det);
}

void gen_random_vector_3(mpz_t vec[3], gmp_randstate_t state, const mpz_t mod)
{
    mpz_t rnd; mpz_init(rnd);
    for( int i=0; i<3; i++ ) {
        mpz_urandomm(rnd, state, mod);
        if( mpz_cmp_ui(rnd, 0) )
            mpz_set(vec[i], rnd);
    }
    mpz_clear(rnd);
}

void gen_random_vector_4(mpz_t vec[4], gmp_randstate_t state, const mpz_t mod)
{
    mpz_t rnd; mpz_init(rnd);
    for( int i=0; i<4; i++ ) {
        mpz_urandomm(rnd, state, mod);
        if( mpz_cmp_ui(rnd, 0) )
            mpz_set(vec[i], rnd);
    }
    mpz_clear(rnd);
}

void get_inverse_matrix_4x4(mpz_t inv[4][4], mpz_t matrix[4][4], mpz_t mod)
{   
    mpz_t det;
    mpz_init(det);

    // 행렬식 계산
    cal_determinant_4x4(det, matrix, 4, mod);

    // 행렬식이 0이면 역행렬이 존재하지 않음
    if (mpz_sgn(det) == 0) {
        printf("Inverse does not exist (determinant is zero).\n");
        mpz_clear(det);
        return;
    }

    // 여인자 행렬을 계산
    mpz_t cofactor[4][4], sub_det;
    mpz_init(sub_det);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mpz_init(cofactor[i][j]);
        }
    }

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            // 소행렬 생성
            int sub_i = 0, sub_j;
            for (int row = 0; row < 4; row++) {
                if (row == i) continue;
                sub_j = 0;
                for (int col = 0; col < 4; col++) {
                    if (col == j) continue;
                    mpz_set(cofactor[sub_i][sub_j], matrix[row][col]);
                    sub_j++;
                }
                sub_i++;
            }

            // 소행렬의 행렬식 계산
            cal_determinant_4x4(sub_det, cofactor, 3, mod);
            // 여인자 계산
            mpz_set(inv[j][i], sub_det);  // transposing
            if ((i + j) % 2 == 1) {
                mpz_neg(inv[j][i], inv[j][i]);  // cofactor의 부호 조정
            }
        }
    }

    // 역행렬 계산
    mpz_t det_inv;
    mpz_init(det_inv);
    mpz_invert(det_inv, det, mod);  // det의 역수 계산

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mpz_mul(inv[i][j], inv[i][j], det_inv);  // (1/det) * cofactor
            mpz_mod(inv[i][j], inv[i][j], mod);  // mod 연산
        }
    }

    mpz_clear(det);
    mpz_clear(det_inv);

    for (int i = 0; i < 4; i++) {
        for (int j = 0; j < 4; j++) {
            mpz_clear(cofactor[i][j]);
        }
    }
}


/*----------------------------------( 삭제하면 안됨 !!! 시작 : 테스트 벡터를 위해서 일시적으로 생성한 함수 )-----------------------------*/
void divide_values(mpz_t num, mpz_t num_values[128], unsigned int seed)
{
    mpz_t average_value, remaining, delta, min_v, max_v, temp_min, temp_max;
    mpz_inits(average_value, remaining, delta, temp_min, temp_max, min_v, max_v, NULL);

    gmp_randstate_t state; gmp_randinit_default(state); gmp_randseed_ui(state, seed);


    mpz_fdiv_q_ui(average_value, num, 128); // 평균값 계산
    
    // min_v = average_value - (average_value//10)
    mpz_fdiv_q_ui(min_v, average_value, 10);
    mpz_sub(min_v, average_value, min_v);
    // max_v = average_value + (average_value//10)
    mpz_fdiv_q_ui(max_v, average_value, 10);
    mpz_add(max_v, average_value, max_v);

    while( 1 ) {
        mpz_set(remaining, num);  // remaining = num

        for( int i = 0; i < 127; i++ ) {
            mpz_fdiv_q_ui(delta, average_value, 10);
            mpz_urandomm(delta, state, delta);

            // delta가 양수 혹은 음수가 되도록 설정
            if( rand() % 2 == 0 ) {
                mpz_neg(delta, delta);
            }

            // num_values[i] = average_value + delta
            mpz_add(num_values[i], average_value, delta);
            mpz_sub(remaining, remaining, num_values[i]);

            // ki가 0이하인 경우 1로 설정
            if( mpz_sgn(num_values[i]) <= 0 ) {
                mpz_set_ui(num_values[i], 1);
                mpz_sub_ui(remaining, remaining, 1);
            }
        }
        // 마지막 값 설정
        mpz_set(num_values[127], remaining);

        // 조건 체크: min(num_values) >= min_v/2 && max(num_values) >= max_v*2
        int valid = 1;
        for( int i = 0; i < 128; i++ ) {
            mpz_fdiv_q_ui(temp_min, min_v, 2);
            mpz_mul_ui(temp_max, max_v, 2);

            if( mpz_cmp(num_values[i], temp_min) < 0 || mpz_cmp(num_values[i], temp_max) > 0 ) {
                valid = 0;
                break;
            }
        }
        if( valid ) break;
    }

    mpz_clears(average_value, remaining, delta, temp_min, temp_max, min_v, max_v, NULL);
}


void write_values(const char* filename, mpz_t k_values[128], mpz_t d_values[128])
{
    FILE *file = fopen(filename, "r+");
    if( file == NULL ) {
        perror("파일 열기 실패\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    long insert_pos = -1;
    while( fgets(line, sizeof(line), file) ) {
        if( strstr(line, "const char* k_values_str[128] = ")) {
            insert_pos = ftell(file);
            break;
        }
    }

    if( insert_pos == -1 ) {
        fprintf(stderr, "k_values_str 선언을 파일에서 찾을 수 없음.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fseek(file, insert_pos, SEEK_SET);

    fprintf(file, "\n{");
    for( int i=0; i<128; i++ ) {
        if( i % 8 == 0) fprintf(file, "\n\t");
        char* str1 = mpz_get_str(NULL, 16, k_values[i]);
        fprintf(file, "\"%064s\", ", str1);
        free(str1);
        // gmp_fprintf(file, "\"%064ZX\", ", k_values[i]);
    }
    fprintf(file, "\n};");

    // d_values 쓰기
    fprintf(file, "\nconst char* d_values_str[128] =\n");
    fprintf(file, "\n{");
    for( int i=0; i<128; i++ ) {
        if( i % 8 == 0) fprintf(file, "\n\t");
        char* str1 = mpz_get_str(NULL, 16, d_values[i]);
        fprintf(file, "\"%064s\", ", str1);
        free(str1);
        // gmp_printf("%064ZX\n", d_values[i]);
    }
    fprintf(file, "\n};");

    fclose(file);
}
/*----------------------------------( 끝 : 테스트 벡터를 위해서 일시적으로 생성한 함수 )-----------------------------*/


void cal_ki_G_with_table(const char* k_i_str[128][4], Point G, mpz_t e, mpz_t table[128][4][2])
{   
    int i;

    /* k_i mpz_t 자료형으로 변환 */
    mpz_t ki[128][4]; init_table((mpz_t *)ki,128,4);
        /* Convert ki_str to mpz_t values */
    for (i = 0; i < 128; i++) {
        for (int j = 0; j < 4; j++) {
            mpz_set_str(ki[i][j], k_i_str[i][j], 16);
        }
    }
    /*
    int bit;
    for( i = 0; i < 128; i++ ) {
        bit = get_two_bits(e, (i << 1));
        mpz_set_str(ki[i], ki_str[i][bit], 16);
    }
    */
    /* 포인트 계산을 위한 변수 설정 */
    Point ki_G[128][4]; /* for( int i = 0; i < 256; i++ ) point_init(ki_G[i]); */

    /* ki_G[i] = [k_i]G mod p 계산 */
    for( i = 0; i < 128; i++ ) {
        for(int j=0;j<4;j++){
        point_init(ki_G[i][j]);
        point_scalar(ki_G[i][j], G, ki[i][j], 256, a, p); /* ki_G[i] = [k_i]G mod p */
        }
        
    }


    
    for( i = 0; i < 128; i++ ) {
        for(int j=0;j<4;j++){
        // bit = get_two_bits(e, (i << 1));
        mpz_set(table[i][j][0], ki_G[i][j]->x);
        mpz_set(table[i][j][1], ki_G[i][j]->y);
        }
    }

    clear_table((mpz_t *)ki,128,4);
}

void write_ki_G_table(const char* filename, mpz_t table[128][4][2])
{
    FILE *file = fopen(filename, "r+");
    if( file == NULL ) {
        perror("파일 열기 실패\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    long insert_pos = -1;
    while (fgets(line, sizeof(line), file)) {
        if( strstr(line, "const char* ki_G_str[128][4][2] =")) {
            insert_pos = ftell(file);
            break;
        }
    }

    if( insert_pos == -1 ) {
        fprintf(stderr, "ki_G_str 선언을 파일에서 찾을 수 없음.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fseek(file, insert_pos, SEEK_SET);

    fprintf(file, "\n{\n");
    for( int i=0; i<128; i++ ) {
        fprintf(file, "\t{ ");
        for( int j=0; j<4; j++ ) {
            char* num_str1 = mpz_get_str(NULL, 16, table[i][j][0]);
            char* num_str2 = mpz_get_str(NULL, 16, table[i][j][1]);
            fprintf(file, "{ \"%064s\", ", num_str1);
            fprintf(file, "\"%064s\" }, ", num_str2);
            free(num_str1); free(num_str2);
        }
        fprintf(file, " }, \n");
    }
    fprintf(file, "};\n");

    fclose(file);
}




void write_last_encoding_A(const char* filename, mpz_t table[3][3])
{
    FILE *file = fopen(filename, "r+");
    if (file == NULL) {
        perror("파일 열기 실패\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    long insert_pos = -1;
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "const char* encoding_A_last[3][3] =")) {
            insert_pos = ftell(file);
            break;
        }
    }

    if (insert_pos == -1) {
        fprintf(stderr, "encoding_A_last 선언을 파일에서 찾을 수 없음.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fseek(file, insert_pos, SEEK_SET);

    
    fprintf(file, "\n{\n");
    for (int i = 0; i < 3; i++) {
        fprintf(file, "\t{ ");
        for (int j = 0; j < 3; j++) {
            // mpz 값을 16진수 문자열로 변환
            char* num_str = mpz_get_str(NULL, 16, table[i][j]);
            fprintf(file, "\"%064s\"", num_str);

            if (j < 3) {  // 마지막 열 이후에만 쉼표를 추가
                fprintf(file, ", ");
            }
            free(num_str);
        }
        fprintf(file, " }");
        if (i < 3) {  // 마지막 행 이후에만 쉼표와 줄바꿈을 추가
            fprintf(file, ",\n");
        }
    }
    fprintf(file, "\n};\n");

    fclose(file);
}

void write_last_encoding_B(const char* filename, mpz_t array[3])
{
    FILE *file = fopen(filename, "r+");
    if (file == NULL) {
        perror("파일 열기 실패\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    long insert_pos = -1;
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "const char* encoding_B_last[3] =")) {
            insert_pos = ftell(file);
            break;
        }
    }

    if (insert_pos == -1) {
        fprintf(stderr, "encoding_B_last 선언을 파일에서 찾을 수 없음.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fseek(file, insert_pos, SEEK_SET);

    fprintf(file, "\n{\n");
    for (int i = 0; i < 3; i++) {
        char* num_str = mpz_get_str(NULL, 16, array[i]);
        fprintf(file, "\t\"%064s\"", num_str); // 16진수로 변환된 값 출력

        free(num_str); // 메모리 해제

        if (i < 2) { // 마지막 항목 뒤에는 쉼표 추가하지 않음
            fprintf(file, ",");
        }
        fprintf(file, "\n");
    }
    fprintf(file, "};\n");

    fclose(file);
}

void write_u_input(const char* filename, mpz_t array[4])
{
    FILE *file = fopen(filename, "r+");
    if (file == NULL) {
        perror("파일 열기 실패\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    long insert_pos = -1;
    while (fgets(line, sizeof(line), file)) {
        if (strstr(line, "const char* u_input[4] =")) {
            insert_pos = ftell(file);
            break;
        }
    }

    if (insert_pos == -1) {
        fprintf(stderr, "u_input[4] 선언을 파일에서 찾을 수 없음.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fseek(file, insert_pos, SEEK_SET);

    fprintf(file, "\n{\n");
    for (int i = 0; i < 4; i++) {
        char* num_str = mpz_get_str(NULL, 16, array[i]);
        fprintf(file, "\t\"%064s\"", num_str); // 16진수로 변환된 값 출력

        free(num_str); // 메모리 해제

        if (i < 4) { // 마지막 항목 뒤에는 쉼표 추가하지 않음
            fprintf(file, ",");
        }
        fprintf(file, "\n");
    }
    fprintf(file, "};\n");

    fclose(file);
}

void tmp_k_d(mpz_t d_i_table[128][4],mpz_t k_i_table[128][4]){
    int bit;
    
    mpz_t k_value_table[128]; init_array((mpz_t*)k_value_table, 128);
    for( int i=0; i<128; i++ ) {
        mpz_set_str(k_value_table[i], k_values_str[i], 16);
    }
    for( int i=0; i<128; i++ ) {
        bit = get_two_bits(e, (i << 1));
        // printf("i = %d, bit = %d\n", i, bit);
        for( int j=0; j<4; j++ ) {
            if( j == bit ) {
                mpz_set_str(   k_i_table[i][bit],k_values_str[i], 16);
                mpz_set_str(d_i_table[i][bit],d_values_str[i], 16);
            }
            
    // printf("ddd\n"); fflush(stdout);
        }
    }
    clear_array((mpz_t*)k_value_table, 128);
}

void make_and_write_ki_di_alpha_beta_table(const char* filename, mpz_t d_i_table[128][4], mpz_t alpha_i_table[128][4], mpz_t beta_i_table[128][4],mpz_t k_i_table[128][4]){
    
    gmp_randstate_t state; gmp_randinit_default(state); gmp_randseed_ui(state, 0);   // 현재 시간을 시드로 설정
    for(int i=0;i<128;i++){
        for(int j=0;j<4;j++){
            mpz_urandomm(k_i_table[i][j], state, n);
            mpz_urandomm(alpha_i_table[i][j], state, n);
            mpz_urandomm(beta_i_table[i][j], state, n);
            mpz_urandomm(d_i_table[i][j],state,n);
        }
    }
    /*test vector 확인을 위해 일시적으로 추가되는 코드*/
    tmp_k_d(d_i_table, k_i_table);
    
    
    FILE *file = fopen(filename, "r+");
        if( file == NULL ) {
            perror("파일 열기 실패\n");
            exit(EXIT_FAILURE);
        }

        char line[256];
        long insert_pos = -1;
        while( fgets(line, sizeof(line), file) ) {
            if( strstr(line, "const char* k_i_str[128][4] =") ) {
                insert_pos = ftell(file);
                break;
            }
        }

        if( insert_pos == -1 ) {
            fprintf(stderr, "k_i_str 선언을 찾을 수 없음.\n");
            fclose(file);
            exit(EXIT_FAILURE);
        }

        fseek(file, insert_pos, SEEK_SET);

        fprintf(file, "\n{\n");
        for( int i=0; i<128; i++ ) {
            fprintf(file, "\t{ ");
            for( int j=0; j<4; j++ ) {
                char* str1 = mpz_get_str(NULL, 16, k_i_table[i][j]);
                fprintf(file, "\"%064s\", ", str1);
                free(str1);
            }
            fprintf(file, " }, \n");
        }
        fprintf(file, "};\n");


        fprintf(file, "\n\nconst char* d_i_str[128][4] =");
        fprintf(file, "\n{\n");
        for( int i=0; i<128; i++ ) {
            fprintf(file, "\t{ ");
            for( int j=0; j<4; j++ ) {
                char* str1 = mpz_get_str(NULL, 16, d_i_table[i][j]);
                fprintf(file, "\"%064s\", ", str1);
                free(str1);
            }
            fprintf(file, " }, \n");
        }
        fprintf(file, "};\n");




        fprintf(file, "\n\nconst char* alpha_i_str[128][4] =");
        fprintf(file, "\n{\n");
        for( int i=0; i<128; i++ ) {
            fprintf(file, "\t{ ");
            for( int j=0; j<4; j++ ) {
                char* str1 = mpz_get_str(NULL, 16, alpha_i_table[i][j]);
                fprintf(file, "\"%064s\", ", str1);
                free(str1);
            }
            fprintf(file, " }, \n");
        }
        fprintf(file, "};\n");

        fprintf(file, "\n\nconst char* beta_i_str[128][4] =");
        fprintf(file, "\n{\n");
        for( int i=0; i<128; i++ ) {
            fprintf(file, "\t{ ");
            for( int j=0; j<4; j++ ) {
                char* str1 = mpz_get_str(NULL, 16, beta_i_table[i][j]);
                fprintf(file, "\"%064s\", ", str1);
                free(str1);
            }
            fprintf(file, " }, \n");
        }
        fprintf(file, "};\n");

        fclose(file);





}



void cal_6_x_e1(mpz_t xi_table[4], mpz_t x_e1)
{   
    // 6*x_e1 = -11*x_i0 + 18*x_i1 - 9*x_i2 + 2*x_i3
    mpz_t temp1, temp2, temp3, temp4; 
    mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);

    mpz_mul_ui(temp1, xi_table[1], 18); mpz_mod(temp1, temp1, n); // 18*x_i1
    mpz_mul_ui(temp2, xi_table[0], 11); mpz_mod(temp2, temp2, n); // 11*x_i0
    mpz_sub(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n); // -11*x_i0 + 18*x_i1

    mpz_mul_ui(temp3, xi_table[3], 2);  mpz_mod(temp3, temp3, n); // 2*x_i3
    mpz_mul_ui(temp4, xi_table[2], 9);  mpz_mod(temp4, temp4, n); // 9*x_i2
    mpz_sub(temp3, temp3, temp4);       mpz_mod(temp3, temp3, n); // -9*x_i2 + 2*x_i3

    mpz_add(x_e1, temp1, temp3);        mpz_mod(x_e1, x_e1, n);    // -11*x_i0 + 18*x_i1 - 9*x_i2 + 2*x_i3

    mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); mpz_clear(temp4);
}

void cal_6_x_e2(mpz_t xi_table[4], mpz_t x_e2)
{   
    // 6*x_e2 = 6*x_i0 - 15*x_i1 + 12*x_i2 - 3*x_i3
    mpz_t temp1, temp2, temp3, temp4; 
    mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4);

    mpz_mul_ui(temp1, xi_table[0], 6);  mpz_mod(temp1, temp1, n);   // 6*x_i0
    mpz_mul_ui(temp2, xi_table[1], 15); mpz_mod(temp2, temp2, n);   // 15*x_i1
    mpz_sub(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);   // 6*x_i0 - 15*x_i1

    mpz_mul_ui(temp3, xi_table[2], 12); mpz_mod(temp3, temp3, n);   // 12*x_i2
    mpz_mul_ui(temp4, xi_table[3], 3);  mpz_mod(temp4, temp4, n);   // 3*x_i3
    mpz_sub(temp3, temp3, temp4);       mpz_mod(temp3, temp3, n);   // 12*x_i2 - 3*x_i3

    mpz_add(x_e2, temp1, temp3);        mpz_mod(x_e2, x_e2, n);     // 6*x_i0 - 15*x_i1 +  12*x_i2 - 3*x_i3

    mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); mpz_clear(temp4);
}

void cal_6_x_e3(mpz_t xi_table[4], mpz_t x_e3)
{   
    // 6*x_e3 = -x_i0 + 3*x_i1 - 3*x_i2 + x_i3
    mpz_t temp1, temp2; 
    mpz_init(temp1); mpz_init(temp2);

    mpz_mul_ui(temp1, xi_table[1], 3);  mpz_mod(temp1, temp1, n);   // 3*x_i1
    mpz_sub(temp1, temp1, xi_table[0]); mpz_mod(temp1, temp1, n);   // -x_i0 + 3*x_i1

    mpz_mul_ui(temp2, xi_table[2], 3);  mpz_mod(temp2, temp2, n);   // 3*x_i2
    mpz_sub(temp2, xi_table[3], temp2); mpz_mod(temp2, temp2, n);   // -3*x_i2 + x_i3

    mpz_add(x_e3, temp1, temp2);        mpz_mod(x_e3, x_e3, n);     // -x_i0 + 3*x_i1 - 3*x_i2 + x_i3

    mpz_clear(temp1); mpz_clear(temp2);
}


void gen_eq_coeff_table( mpz_t m[128][4][4], mpz_t a[128][4][4], mpz_t b[128][4], mpz_t pre_a[128][4][4], mpz_t pre_b[128][4], 
                        mpz_t m_last[3][3], mpz_t a_last[3][3], mpz_t b_last[3],
                        mpz_t ki_table[128][4], mpz_t alpha_i_table[128][4], mpz_t beta_i_table[128][4], mpz_t di_table[128],
                        struct EQ_COEFF_TABLE *eq_coeff_table,mpz_t n)
{   
    int round, it, i, j, k, \
        bit;
    /*
        round -> 총 128라운드
        it -> 4개의 식
        i, j ... -> 필요할 떄 사용
    */
    mpz_t temp1, temp2, temp3, temp4, temp5, temp6, temp7, res;
    mpz_init(temp1); mpz_init(temp2); mpz_init(temp3); mpz_init(temp4); mpz_init(temp5); mpz_init(temp6); mpz_init(temp7); mpz_init(res);

    for( round = 0; round < 127; round++ ) 
    {   
        // bit = get_two_bits(hash_e, (round << 1));
        // printf("\nround = %d\n", round);
        for( it = 0; it < 4; it++ ) 
        {   
            // 1) vj : 4개
            for( i = 0; i < 4; i++ ) {
                mpz_mul(temp1, m[round][it][0], a[round][0][i]);  mpz_mod(temp1, temp1, n);
                mpz_mul(temp2, m[round][it][1], a[round][1][i]);  mpz_mod(temp2, temp2, n);
                mpz_mul(temp3, m[round][it][2], a[round][2][i]);  mpz_mod(temp3, temp3, n);
                mpz_mul(temp4, m[round][it][3], a[round][3][i]);  mpz_mod(temp4, temp4, n);

                // gmp_printf("m_%d0 * a_0%d = %064ZX\n", it, i, temp1);
                // gmp_printf("m_%d1 * a_1%d = %064ZX\n", it, i, temp2);
                // gmp_printf("m_%d2 * a_2%d = %064ZX\n", it, i, temp3);
                // gmp_printf("m_%d3 * a_3%d = %064ZX\n", it, i, temp4);


                mpz_add(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);
                mpz_add(temp3, temp3, temp4);       mpz_mod(temp3, temp3, n);
                mpz_add(res, temp1, temp3);         mpz_mod(res, res, n);
                // mpz_neg(res, res);                  mpz_mod(res, res, n);

                // gmp_printf("v_%d의 계수 = %064ZX\n\n", i, res);

                mpz_set(eq_coeff_table->coeff_m[round][it][i], res);
            }
            
            // 2) uj : 4개
            for( i = 0; i < 4; i++ ) {
                mpz_mul(temp1, m[round][it][0], pre_a[round][0][i]);  mpz_mod(temp1, temp1, n);
                mpz_mul(temp2, m[round][it][1], pre_a[round][1][i]);  mpz_mod(temp2, temp2, n);
                mpz_mul(temp3, m[round][it][2], pre_a[round][2][i]);  mpz_mod(temp3, temp3, n);
                mpz_mul(temp4, m[round][it][3], pre_a[round][3][i]);  mpz_mod(temp4, temp4, n);

                mpz_add(temp1, temp1, temp2);           mpz_mod(temp1, temp1, n);
                mpz_add(temp3, temp3, temp4);           mpz_mod(temp3, temp3, n);
                mpz_add(res, temp1, temp3);             mpz_mod(res, res, n);

                mpz_set(eq_coeff_table->coeff_m[round][it][4 + i], res);
            }

            // 3) (ei)^3 : 1개 -> 6을 곱해서 저장
            // 저장식 : m_j1(-k_i0 + 3k_i1 - 3k_i2 + k_i3) + m_j3(-\alpha_i0 + 3\alpha_i1 - 3\alpha_i2 + \alpha_i3) + m_j4(-\beta_i0 + 3\beta_i1 - 3\beta_i2 + \beta_i3)
            mpz_mul_ui(temp1, ki_table[round][1], 3);       mpz_mod(temp1, temp1, n);   // temp1 = 3k_i1
            mpz_sub(temp1, temp1, ki_table[round][0]);      mpz_mod(temp1, temp1, n);   // temp1 = 3k_i1 - k_i0
            mpz_mul_ui(temp2, ki_table[round][2], 3);       mpz_mod(temp2, temp2, n);   // temp2 = 3k_i2
            mpz_sub(temp2, ki_table[round][3], temp2);      mpz_mod(temp2, temp2, n);   // temp2 = k_i3 - 3k_i2
            mpz_add(temp1, temp1, temp2);                   mpz_mod(temp1, temp1, n);   // temp1 = 3k_i1 - k_i0 + k_i3 - 3k_i2
            mpz_mul(temp1, temp1, m[round][it][0]);         mpz_mod(temp1, temp1, n);   // temp1 = m_j1(3k_i1 - k_i0 + k_i3 - 3k_i2)

            mpz_mul_ui(temp2, alpha_i_table[round][1], 3);  mpz_mod(temp2, temp2, n);
            mpz_sub(temp2, temp2, alpha_i_table[round][0]); mpz_mod(temp2, temp2, n);
            mpz_mul_ui(temp3, alpha_i_table[round][2], 3);  mpz_mod(temp3, temp3, n);
            mpz_sub(temp3, alpha_i_table[round][3], temp3); mpz_mod(temp3, temp3, n);
            mpz_add(temp2, temp2, temp3);                   mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, temp2, m[round][it][2]);         mpz_mod(temp2, temp2, n);

            mpz_mul_ui(temp3, beta_i_table[round][1], 3);   mpz_mod(temp3, temp3, n);
            mpz_sub(temp3, temp3, beta_i_table[round][0]);  mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, beta_i_table[round][2], 3);   mpz_mod(temp4, temp4, n);
            mpz_sub(temp4, beta_i_table[round][3], temp4);  mpz_mod(temp4, temp4, n);
            mpz_add(temp3, temp3, temp4);                   mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, temp3, m[round][it][3]);         mpz_mod(temp3, temp3, n);

            mpz_add(temp1, temp1, temp2);                   mpz_mod(temp1, temp1, n);
            mpz_add(res, temp1, temp3);                     mpz_mod(res, res, n);

            mpz_set(eq_coeff_table->coeff_m[round][it][8], res);

            // 4) (ei)^2 : 1개 -> 6을 곱해서 저장
            // 저장식 : m_j1(6k_i0 - 15k_i1 + 12k_i2 - 3k_i3) + 
            //        m_j3(6\alpha_i0 - 15\alpha_i1 + 12\alpha_i2 - 3\alpha_i3) + 
            //        m_j4(6\beta_i0 - 15\beta_i1 + 12\beta_i2 - 3\beta_i3)
            mpz_mul_ui(temp1, ki_table[round][1], 15);  mpz_mod(temp1, temp1, n);   // temp1 = 15k_i1
            mpz_mul_ui(temp2, ki_table[round][0], 6);   mpz_mod(temp2, temp2, n);   // temp2 = 6k_i0
            mpz_sub(temp1, temp2, temp1);               mpz_mod(temp1, temp1, n);   // temp1 = 6k_i0 - 15k_i1
            mpz_mul_ui(temp2, ki_table[round][2], 12);  mpz_mod(temp2, temp2, n);   // temp2 = 12k_i2
            mpz_mul_ui(temp3, ki_table[round][3], 3);   mpz_mod(temp3, temp3, n);   // temp3 = 3k_i3
            mpz_sub(temp2, temp2, temp3);               mpz_mod(temp2, temp2, n);   // temp2 = 12k_i2 - 3k_i0
            mpz_add(temp1, temp1, temp2);               mpz_mod(temp1, temp1, n);   // temp1 = 6k_i0 - 15k_i1 + 12k_i2 - 3k_i0
            mpz_mul(temp1, m[round][it][0], temp1);     mpz_mod(temp1, temp1, n);   // temp1 = m_j1(6k_i0 - 15k_i1 + 12k_i2 - 3k_i0)

            mpz_mul_ui(temp2, alpha_i_table[round][1], 15); mpz_mod(temp2, temp2, n);
            mpz_mul_ui(temp3, alpha_i_table[round][0], 6);  mpz_mod(temp3, temp3, n);
            mpz_sub(temp2, temp3, temp2);                   mpz_mod(temp2, temp2, n);
            mpz_mul_ui(temp3, alpha_i_table[round][2], 12); mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, alpha_i_table[round][3], 3);  mpz_mod(temp4, temp4, n);
            mpz_sub(temp3, temp3, temp4);                   mpz_mod(temp3, temp3, n);
            mpz_add(temp2, temp2, temp3);                   mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, m[round][it][2], temp2);         mpz_mod(temp2, temp2, n);
            
            mpz_mul_ui(temp3, beta_i_table[round][1], 15);  mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, beta_i_table[round][0], 6);   mpz_mod(temp4, temp4, n);
            mpz_sub(temp3, temp4, temp3);                   mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, beta_i_table[round][2], 12);  mpz_mod(temp4, temp4, n);
            mpz_mul_ui(temp5, beta_i_table[round][3], 3);   mpz_mod(temp5, temp5, n);
            mpz_sub(temp4, temp4, temp5);                   mpz_mod(temp4, temp4, n);
            mpz_add(temp3, temp3, temp4);                   mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, m[round][it][3], temp3);         mpz_mod(temp3, temp3, n);

            mpz_add(temp1, temp1, temp2);                   mpz_mod(temp1, temp1, n);
            mpz_add(res, temp1, temp3);                     mpz_mod(res, res, n);
            
            mpz_set(eq_coeff_table->coeff_m[round][it][9], res);
            // gmp_printf("result[1] = %064ZX\n", eq_coeff_table->coeff_m[round][it][9]);

            // 5) ei : 1개 -> 6을 곱해서 저장
            // 저장식 : m_j1(-11k_i0 + 18k_i1 - 9k_i2 + 2k_i3) + 
            //        m_j3(-11\alpha_i0 + 18\alpha_i1 - 9\alpha_i2 + 2\alpha_i3) +
            //        m_j4(-11\beta_i0 + 18\beta_i1 - 9\beta_i2 + 2\beta_i3)
            mpz_mul_ui(temp1, ki_table[round][0], 11);             mpz_mod(temp1, temp1, n);   // temp1 = 11k_i0
            mpz_mul_ui(temp2, ki_table[round][1], 18);             mpz_mod(temp2, temp2, n);   // temp2 = 18k_i1
            mpz_sub(temp1, temp2, temp1);                          mpz_mod(temp1, temp1, n);   // temp1 = 18k_i1 - 11k_i0
            mpz_mul_ui(temp2, ki_table[round][2], 9);              mpz_mod(temp2, temp2, n);   // temp2 = 9k_i2
            mpz_mul_ui(temp3, ki_table[round][3], 2);              mpz_mod(temp3, temp3, n);   // temp3 = 2k_i3
            mpz_sub(temp2, temp3, temp2);                          mpz_mod(temp2, temp2, n);   // temp2 = 2k_i3 - 9k_i2
            mpz_add(temp1, temp1, temp2);                          mpz_mod(temp1, temp1, n);   // temp1 = 18k_i1 - 11k_i0 + 2k_i3 - 9k_i2
            mpz_mul(temp1, m[round][it][0], temp1);                mpz_mod(temp1, temp1, n);   // temp1 = m_j1(18k_i1 - 11k_i0 + 2k_i3 - 9k_i2)

            mpz_mul_ui(temp2, alpha_i_table[round][0], 11);        mpz_mod(temp2, temp2, n);
            mpz_mul_ui(temp3, alpha_i_table[round][1], 18);        mpz_mod(temp3, temp3, n);
            mpz_sub(temp2, temp3, temp2);                          mpz_mod(temp2, temp2, n);
            mpz_mul_ui(temp3, alpha_i_table[round][2], 9);         mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, alpha_i_table[round][3], 2);         mpz_mod(temp4, temp4, n);
            mpz_sub(temp3, temp4, temp3);                          mpz_mod(temp3, temp3, n);
            mpz_add(temp2, temp2, temp3);                          mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, m[round][it][2], temp2);                mpz_mod(temp2, temp2, n);

            mpz_mul_ui(temp3, beta_i_table[round][0], 11);         mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, beta_i_table[round][1], 18);         mpz_mod(temp4, temp4, n);
            mpz_sub(temp3, temp4, temp3);                          mpz_mod(temp3, temp3, n);
            mpz_mul_ui(temp4, beta_i_table[round][2], 9);          mpz_mod(temp4, temp4, n);
            mpz_mul_ui(temp5, beta_i_table[round][3], 2);          mpz_mod(temp5, temp5, n);
            mpz_sub(temp4, temp5, temp4);                          mpz_mod(temp4, temp4, n);
            mpz_add(temp3, temp3, temp4);                          mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, m[round][it][3], temp3);                mpz_mod(temp3, temp3, n);

            mpz_add(temp1, temp1, temp2);                          mpz_mod(temp1, temp1, n);
            mpz_add(res, temp1, temp3);                            mpz_mod(res, res, n);

            mpz_set(eq_coeff_table->coeff_m[round][it][10], res);
            // gmp_printf("result[2] = %064ZX\n", eq_coeff_table->coeff_m[round][it][10]);

            // 6) constant : 1개
            
            mpz_add(temp1, ki_table[round][0], pre_b[round][0]);        mpz_mod(temp1, temp1, n);
            /**(been 수정)**/
            mpz_sub(temp1, temp1, b[round][0]);                         mpz_mod(temp1, temp1, n);
            /**(끝)**/
            mpz_mul(temp1, temp1, m[round][it][0]);                     mpz_mod(temp1, temp1, n);

            // gmp_printf("pre_b[%d][0] = %064ZX\n", round, pre_b[round][0]);
            // gmp_printf("b[%d][0] = %064ZX\n", round, b[round][0]);
            // gmp_printf("m_%d0 = %064ZX\n", it, temp1);

            mpz_add(temp2, pre_b[round][1], di_table[round]);           mpz_mod(temp2, temp2, n);
            mpz_sub(temp2, temp2, b[round][1]);                         mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, temp2, m[round][it][1]);                     mpz_mod(temp2, temp2, n);
        
            // gmp_printf("pre_b[%d][1] = %064ZX\n", round, pre_b[round][1]);
            // gmp_printf("b[%d][1] = %064ZX\n", round, b[round][1]);
            // gmp_printf("m_%d1 = %064ZX\n", it, temp2);
            
            mpz_add(temp3, alpha_i_table[round][0], pre_b[round][2]);   mpz_mod(temp3, temp3, n);
            mpz_sub(temp3, temp3, b[round][2]);                         mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, temp3, m[round][it][2]);                     mpz_mod(temp3, temp3, n);
            
            // gmp_printf("pre_b[%d][2] = %064ZX\n", round, pre_b[round][2]);
            // gmp_printf("b[%d][2] = %064ZX\n", round, b[round][2]);
            // gmp_printf("m_%d2 = %064ZX\n", it, temp3);

            mpz_add(temp4, beta_i_table[round][0], pre_b[round][3]);    mpz_mod(temp4, temp4, n);
            mpz_sub(temp4, temp4, b[round][3]);                         mpz_mod(temp4, temp4, n);
            mpz_mul(temp4, temp4, m[round][it][3]);                     mpz_mod(temp4, temp4, n);

            // gmp_printf("m_%d3 = %064ZX\n", it, temp4);

            mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n);
            mpz_add(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);
            mpz_add(res, temp1, temp3);                                 mpz_mod(res, res, n);

            mpz_set(eq_coeff_table->coeff_m[round][it][11], res);

            
        }
    }

    //==================================================================================================
    // 중간 디버깅 출력
    //==================================================================================================
    // for( int i=0; i<127; i++ ) {
    //     printf("round = %d\n", i);
    //     for( int j=0; j<4; j++ ) {
    //         for( int k=0; k<12; k++ ) {
    //             gmp_printf("coeff[%d][%d][%d] = %ZX\n", i, j, k, eq_coeff_table->coeff_m[i][j][k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("==============\n");
    // }
    //==================================================================================================
/****************************************************************************************************************/
    // 마지막 128 라운드
    // 필요한 변수 계산
    mpz_t k_e1, k_e2, k_e3, alpha_e1, alpha_e2, alpha_e3, beta_e1, beta_e2, beta_e3;
    mpz_init(    k_e1); mpz_init(    k_e2); mpz_init(    k_e3);
    mpz_init(alpha_e1); mpz_init(alpha_e2); mpz_init(alpha_e3);
    mpz_init( beta_e1); mpz_init( beta_e2); mpz_init( beta_e3);
    
    cal_6_x_e1(     ki_table[round],     k_e1); cal_6_x_e2(     ki_table[round],     k_e2); cal_6_x_e3(     ki_table[round],     k_e3);
    cal_6_x_e1(alpha_i_table[round], alpha_e1); cal_6_x_e2(alpha_i_table[round], alpha_e2); cal_6_x_e3(alpha_i_table[round], alpha_e3);
    cal_6_x_e1( beta_i_table[round],  beta_e1); cal_6_x_e2( beta_i_table[round],  beta_e2); cal_6_x_e3( beta_i_table[round],  beta_e3);
    
    for( it = 0; it < 3; it++ ) {
        // 1) vj : 3개
        for( i = 0; i < 3; i++ ) {
            mpz_mul(temp1, m_last[it][0], a_last[0][i]);    mpz_mod(temp1, temp1, n);
            mpz_mul(temp2, m_last[it][1], a_last[1][i]);    mpz_mod(temp2, temp2, n);
            mpz_mul(temp3, m_last[it][2], a_last[2][i]);    mpz_mod(temp3, temp3, n);
            mpz_add(temp1, temp1, temp2);                   mpz_mod(temp1, temp1, n);
            mpz_add(res, temp1, temp3);                     mpz_mod(res, res, n);

            mpz_set(eq_coeff_table->coeff_m_last[it][i], res);
        }

        // 2) (uj)^2 : 4개
    
        for( i = 0 ; i < 4; i++ ) {
            mpz_mul(temp1, pre_a[round][1][i], pre_a[round][2][i]); mpz_neg(temp1, temp1);      // -pre_a_2i * pre_a_3i

            mpz_mul(temp2, pre_a[round][0][i], pre_a[round][3][i]);
            mpz_sub(temp1, temp1, temp2);                           mpz_mod(temp1, temp1, n);

            mpz_mul(temp2, pre_a[round][2][i], pre_a[round][3][i]); mpz_mod(temp2, temp2, n);

            mpz_sub(temp1, temp1, temp2);                           mpz_mod(temp1, temp1, n);
            
            mpz_mul(res, m_last[it][2], temp1);                     mpz_mod(res, res, n);

            mpz_set(eq_coeff_table->coeff_m_last[it][3 + i], res);
        }

        // 3) ujuk : 6개
        int ptr = 7;
        for( i = 0; i < 3; i++ ) {
            for( j = i; j < 4; j++ ) {
                if( i != j ) {
                    mpz_mul(temp1, pre_a[round][1][i], pre_a[round][2][j]); mpz_neg(temp1, temp1); mpz_mod(temp1, temp1, n);
                    mpz_mul(temp2, pre_a[round][1][j], pre_a[round][2][i]); mpz_neg(temp2, temp2); mpz_mod(temp2, temp2, n);
                    mpz_mul(temp3, pre_a[round][0][i], pre_a[round][3][j]); mpz_neg(temp3, temp3); mpz_mod(temp3, temp3, n);
                    mpz_mul(temp4, pre_a[round][0][j], pre_a[round][3][i]); mpz_neg(temp4, temp4); mpz_mod(temp4, temp4, n);

                    mpz_add(temp1, temp1, temp2); mpz_mod(temp1, temp1, n);
                    mpz_add(temp3, temp3, temp4); mpz_mod(temp3, temp3, n);
                    mpz_add(temp1, temp1, temp3); mpz_mod(temp1, temp1, n); // temp1 = -pre_a_21*pre_a_32 - pre_a_22*pre_a_31 - pre_a_11*pre_a_42 - pre_a_12*pre_a_41

                    mpz_mul(temp2, pre_a[round][2][i], pre_a[round][3][j]); mpz_mod(temp2, temp2, n);
                    mpz_mul(temp3, pre_a[round][2][j], pre_a[round][3][i]); mpz_mod(temp3, temp3, n);

                    mpz_sub(temp1, temp1, temp2); mpz_mod(temp1, temp1, n); // 수정(add->sub)
                    mpz_sub(temp1, temp1, temp3); mpz_mod(temp1, temp1, n);//확인!

                    mpz_mul(res, m_last[it][2], temp1); mpz_mod(res, res, n);
                    mpz_set(eq_coeff_table->coeff_m_last[it][ptr], res);
                    ptr++;
                }
                
            }
        }

        // 4) (e_128)^3uj : 4개
        // 저장식 : m_j3[ (-\alpha_i0 + 3\alpha_i1 - 3\alpha_i2 + \alpha_i3)(-pre_a_21 - pre_a_41) +
        //              (-\beta_i0 + 3\beta_i1 - 3\beta_i2 + \beta_i3)(-pre_a_11 - pre_a_31) +
        //              (-k_i0 + 3k_i1 - 3k_i2 + k_i3)(-pre_a_41)   ]
        //      = m_j3[\alpha_e3*(-pre_a_21 - pre_a_41) + beta_e3*(-pre_a_11 - pre_a_31)  + k_e3*(-pre_a_41)]
        for( i = 0; i < 4; i++ ) 
        {
            mpz_add(temp1, pre_a[round][3][i], pre_a[round][1][i]); mpz_mod(temp1, temp1, n);
            mpz_neg(temp1, temp1);                                  mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, alpha_e3, temp1);                        mpz_mod(temp1, temp1, n);

            mpz_add(temp2, pre_a[round][2][i], pre_a[round][0][i]); mpz_mod(temp2, temp2, n);
            mpz_neg(temp2, temp2);                                  mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, beta_e3, temp2);                         mpz_mod(temp2, temp2, n);

            mpz_neg(temp3, pre_a[round][3][i]);                     mpz_mod(temp3, temp3, n);

            mpz_add(temp1, temp1, temp2);                           mpz_mod(temp1, temp1, n);
            mpz_add(temp1, temp1, temp3);                           mpz_mod(temp1, temp1, n);

            mpz_mul(res, m_last[it][2], temp1);                     mpz_mod(res, res, n);
            mpz_set(eq_coeff_table->coeff_m_last[it][13 + i], res);
        }

        // 5) (e_128)^2uj : 4개
        // 저장식 : m_j3[(6\alpha_i0 - 15\alpha_i1 + 12\alpha_i2 - 3\alpha_i3)(-pre_a_21 - pre_a_41) +
        //              (6\beta_i0 - 15\beta_i1 + 12\beta_i2 - 3\beta_i3)(-pre_a_11 - pre_a_31) +
        //              (6k_i0 - 15k_i1 + 12k_i2 - 3k_i3)(-pre_a_41) ]
        //      = m_j3[alpha_e2*(-pre_a_21 - pre_a_41) + beta_e2*(-pre_a_11 - pre_a_31) + k_e2*(-pre_a_41)]
        for( i = 0; i < 4; i++ ) 
        {
            mpz_add(temp1, pre_a[round][3][i], pre_a[round][1][i]); mpz_mod(temp1, temp1, n);
            mpz_neg(temp1, temp1);                                  mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, alpha_e2, temp1);                        mpz_mod(temp1, temp1, n);

            mpz_add(temp2, pre_a[round][2][i], pre_a[round][0][i]); mpz_mod(temp2, temp2, n);
            mpz_neg(temp2, temp2);                                  mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, beta_e2, temp2);                         mpz_mod(temp2, temp2, n);

            mpz_neg(temp3, pre_a[round][3][i]);                     mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, k_e2, temp3);                            mpz_mod(temp3, temp3, n);

            mpz_add(temp1, temp1, temp2);                           mpz_mod(temp1, temp1, n);
            mpz_add(temp1, temp1, temp3);                           mpz_mod(temp1, temp1, n);

            mpz_mul(res, m_last[it][2], temp1);                     mpz_mod(res, res, n);
            mpz_set(eq_coeff_table->coeff_m_last[it][17 + i], res);
        }

        // 6) e_128uj : 4개
        // 저장식 : m_j3[(-11\alpha_i0 + 18\alpha_i1 - 9\alpha_i2 + 2\alpha_i3)(-pre_a_21 - pre_a_41) +
        //              (-11\beta_i0 + 18\beta_i1 - 9\beta_i2 + 2\beta_i3)(-pre_a_11 - pre_a_31) +
        //              (-11k_i0 + 18k_i1 - 9k_i2 + 2k_i3)(-pre_a_41)   ]
        //      = m_j3[alpha_e1*(-pre_a_21 - pre_a_41) + beta_e1*(-pre_a_11 - pre_a_31) + k_e1*(-pre_a_41)]
        for( i = 0; i < 4; i++ )
        {
            mpz_add(temp1, pre_a[round][3][i], pre_a[round][1][i]); mpz_mod(temp1, temp1, n);
            mpz_neg(temp1, temp1);                                  mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, alpha_e1, temp1);                        mpz_mod(temp1, temp1, n);

            mpz_add(temp2, pre_a[round][2][i], pre_a[round][0][i]); mpz_mod(temp2, temp2, n);
            mpz_neg(temp2, temp2);                                  mpz_mod(temp2, temp2, n);
            mpz_mul(temp2, beta_e1, temp2);                         mpz_mod(temp2, temp2, n);

            mpz_neg(temp3, pre_a[round][3][i]);                     mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, k_e1, temp3);                            mpz_mod(temp3, temp3, n);

            mpz_add(temp1, temp1, temp2);                           mpz_mod(temp1, temp1, n);
            mpz_add(temp1, temp1, temp3);                           mpz_mod(temp1, temp1, n);

            mpz_mul(res, temp1, m_last[it][2]);                     mpz_mod(res, res, n);
            mpz_set(eq_coeff_table->coeff_m_last[it][21 + i], res);
        }

        // 7) uj : 4개
        for( i = 0; i < 4; i++ ) {
            mpz_add(temp1, pre_a[round][0][i], pre_a[round][2][i]); mpz_mod(temp1, temp1, n);
            mpz_mul(temp1, temp1, m_last[it][0]);                   mpz_mod(temp1, temp1, n);
            
            mpz_add(temp2, pre_a[round][1][i], pre_a[round][3][i]); mpz_mod(temp2, temp2, n);   // been수정*(sub->add)
            mpz_mul(temp2, temp2, m_last[it][1]);                   mpz_mod(temp2, temp2, n);

            mpz_add(temp3, pre_b[round][2], alpha_i_table[round][0]);   mpz_mod(temp3, temp3, n);
            mpz_add(temp4, pre_a[round][3][i], pre_a[round][1][i]);     mpz_mod(temp4, temp4, n);
            mpz_neg(temp4, temp4);                                      mpz_mod(temp4, temp4, n);
            mpz_mul(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

            mpz_add(temp4, pre_b[round][1], di_table[round]);           mpz_mod(temp4, temp4, n);
            mpz_neg(temp5, pre_a[round][2][i]);                         mpz_mod(temp5, temp5, n);
            mpz_mul(temp4, temp4, temp5);                               mpz_mod(temp4, temp4, n);

            mpz_add(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

            mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);
            mpz_add(temp5, pre_a[round][2][i], pre_a[round][0][i]);     mpz_mod(temp5, temp5, n);
            mpz_neg(temp5, temp5);                                      mpz_mod(temp5, temp5, n);   // 수정(sub->add, neg 추가)
            mpz_mul(temp4, temp4, temp5);                               mpz_mod(temp4, temp4, n);

            mpz_add(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

            /************************************(예진 수정)**************************************/
            mpz_add(temp4, ki_table[round][0], pre_b[round][0]);        mpz_mod(temp4, temp4, n);
            mpz_neg(temp4, temp4);                                      mpz_mod(temp4, temp4, n);
            /**********************************(예진 수정 끝)**************************************/
            mpz_mul(temp4, temp4, pre_a[round][3][i]);                  mpz_mod(temp4, temp4, n);

            mpz_add(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);
            mpz_mul(temp3, temp3, m_last[it][2]);                       mpz_mod(temp3, temp3, n);

            // temp1 + temp2 + temp3
            mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n);
            mpz_add(res, temp1, temp3);                                 mpz_mod(res, res, n);

            mpz_set(eq_coeff_table->coeff_m_last[it][25 + i], res);            
        }

        /************************************(예진 수정)**************************************/
        // 7-1) rho*uj : 4개
        for( i = 0; i < 4; i++ ) {
            mpz_mul(res, m_last[it][2], pre_a[round][3][i]);    mpz_mod(res, res, n);
            mpz_set(eq_coeff_table->coeff_m_last[it][29 + i], res);
        }
        /**********************************(예진 수정 끝)**************************************/

        // 8) (e_128)^6 : 1개
        mpz_mul(temp1,     k_e3, beta_e3);  mpz_mod(temp1, temp1, n);   // k_e3*\beta_e3
        mpz_mul(temp2, alpha_e3, beta_e3);  mpz_mod(temp2, temp2, n);   // \alpha_e3*\beta_e3
        mpz_add(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);   // k_e3*\beta_e3 + \alpha_e3*\beta_e3
        mpz_neg(temp1, temp1);              mpz_mod(temp1, temp1, n);   // -k_e3*\beta_e3 - \alpha_e3*\beta_e3)
        
        mpz_mul(res, temp1, m_last[it][2]); mpz_mod(res, res, n);
        mpz_set(eq_coeff_table->coeff_m_last[it][33], res);

        // 9) (e_128)^5 : 1개
        mpz_mul(temp1, k_e2, beta_e3);      mpz_mod(temp1, temp1, n);
        mpz_mul(temp2, k_e3, beta_e2);      mpz_mod(temp2, temp2, n);
        mpz_add(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);

        mpz_mul(temp2, alpha_e3, beta_e2);  mpz_mod(temp2, temp2, n);
        mpz_mul(temp3, alpha_e2, beta_e3);  mpz_mod(temp3, temp3, n);
        mpz_add(temp2, temp2, temp3);       mpz_mod(temp2, temp2, n);
        
        mpz_add(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);
        mpz_neg(temp1, temp1);              mpz_mod(temp1, temp1, n);

        mpz_mul(res, temp1, m_last[it][2]); mpz_mod(res, res, n);
        mpz_set(eq_coeff_table->coeff_m_last[it][34], res);

        // 10) (e_128)^4 : 1개
        mpz_mul(temp1, k_e1, beta_e3);      mpz_mod(temp1, temp1, n);
        mpz_mul(temp2, k_e2, beta_e2);      mpz_mod(temp2, temp2, n);
        mpz_mul(temp3, k_e3, beta_e1);      mpz_mod(temp3, temp3, n);
        mpz_add(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);
        mpz_add(temp1, temp1, temp3);       mpz_mod(temp1, temp1, n);

        mpz_mul(temp2, alpha_e3, beta_e1);  mpz_mod(temp2, temp2, n);
        mpz_mul(temp3, alpha_e2, beta_e2);  mpz_mod(temp3, temp3, n);
        mpz_mul(temp4, alpha_e1, beta_e3);  mpz_mod(temp4, temp4, n);
        mpz_add(temp2, temp2, temp3);       mpz_mod(temp2, temp2, n);
        mpz_add(temp2, temp2, temp4);       mpz_mod(temp2, temp2, n);

        mpz_add(temp1, temp1, temp2);       mpz_mod(temp1, temp1, n);
        mpz_neg(temp1, temp1);              mpz_mod(temp1, temp1, n);

        mpz_mul(res, temp1, m_last[it][2]); mpz_mod(res, res, n);
        mpz_set(eq_coeff_table->coeff_m_last[it][35], res);

        // 11) (e_128)^3 : 1개
        //  m_j1[6*(-k_i0 + 3k_i1 - 3k_i2 + k_i3 - \alpha_i0 + 3\alpha_i1 - 3\alpha_i2 + \alpha_i3)]
        //= m_j1[6*(k_e3 + alpha_e3)]
        mpz_add(temp1, k_e3, alpha_e3);         mpz_mod(temp1, temp1, n);
        mpz_mul_ui(temp1, temp1, 6);            mpz_mod(temp1, temp1, n);
        mpz_mul(temp1, temp1, m_last[it][0]);   mpz_mod(temp1, temp1, n);

        // m_j2[6*(-\beta_i0 + 3\beta_i1 - 3\beta_i2 + \beta_i3)]
        // = m_j2 * 6*beta_e3
        mpz_mul_ui(temp2, beta_e3, 6);          mpz_mod(temp2 ,temp2, n);
        mpz_mul(temp2, temp2, m_last[it][1]);   mpz_mod(temp2, temp2, n);

        // m_j3[ -6*(pre_b_2 + d_128)*alpha_e3 - 6*k_e3(pre_b_4 + \beta_128,0) - k_e2\beta_e2 - 6*\beta_e3(pre_b_1 + k_128,0)
        //          - 6*\alpha_e3(pre_b_4 + \beta_128,0) - \alpha_e2\beta_e1 - alpha_e1\beta_e2 - 6*\beta_e3(pre_b_3 + \alpha_128,0)]
        mpz_add(temp3, pre_b[round][1], di_table[round]);           mpz_mod(temp3, temp3, n);   // pre_b_2 + d_128
        mpz_mul_ui(temp3, temp3, 6);                                mpz_mod(temp3, temp3, n);   // 6*pre_b_2 + d_128
        mpz_neg(temp3, temp3);                                      mpz_mod(temp3, temp3, n);   // -6*(pre_b_2 + d_128)
        mpz_mul(temp3, temp3, alpha_e3);                            mpz_mod(temp3, temp3, n);   // -(pre_b_2 + d_128)*\alpha_e3

        mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);   // pre_b_4 + \beta_128,0
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);   // 6*(pre_b_4 + \beta_128,0)
        mpz_mul(temp4, temp4, k_e3);                                mpz_mod(temp4, temp4, n);   // 6*k_e3(pre_b_4 + \beta_128,0)

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);   // -6*(pre_b_2 + d_128)*\alpha_e3 - 6*k_e3(pre_b_4 + \beta_128,0)

        mpz_mul(temp4, k_e2, beta_e2);                              mpz_mod(temp4, temp4, n);   // k_e2\beta_e2
        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);   // -6*(pre_b_2 + d_128)*\alpha_e3 - 6*k_e3(pre_b_4 + \beta_128,0) - k_e2\beta_e2

        mpz_add(temp4, pre_b[round][0], ki_table[round][0]);        mpz_mod(temp4, temp4, n);   // pre_b_1 + k_128,0
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);   // 6*(pre_b_1 + k_128,0)
        mpz_mul(temp4, temp4, beta_e3);                             mpz_mod(temp4, temp4, n);   // 6*\beta_e3(pre_b_1 + k_128,0)

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);   // pre_b_4 + \beta_128,0
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);   // 6*(pre_b_4 + \beta_128,0)
        mpz_mul(temp4, temp4, alpha_e3);                            mpz_mod(temp4, temp4, n);   // 6*\alpha_e3(pre_b_4 + \beta_128,0)

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_mul(temp4, alpha_e2, beta_e1);                          mpz_mod(temp4, temp4, n);
        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_mul(temp4, alpha_e1, beta_e2);                          mpz_mod(temp4, temp4, n);
        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][2], alpha_i_table[round][0]);   mpz_mod(temp4, temp4, n);   // pre_b_3 + \alpha_128,0
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);   // 6*(pre_b_3 + \alpha_128,0)
        mpz_mul(temp4, temp4, beta_e3);                             mpz_mod(temp4, temp4, n);   // 6*\beta_e3(pre_b_3 + \alpha_128,0)

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_mul(temp3, temp3, m_last[it][2]);                       mpz_mod(temp3, temp3, n);
        
        // temp1 + temp2 + temp3
        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n);
        mpz_add(res, temp1, temp3);                                 mpz_mod(res, res, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][36], res);

        /************************************(예진 수정)**************************************/
        // 11-1) (e_128^3 * rho) : 1개
        // 6*m_j3*\beta_e3
        mpz_mul(temp1, m_last[it][2], beta_e3);         mpz_mod(temp1, temp1, n);
        mpz_mul_ui(res, temp1, 6);                      mpz_mod(res, res, n);
        
        mpz_set(eq_coeff_table->coeff_m_last[it][37], res);
        /**********************************(예진 수정 끝)**************************************/

        // 12) (e_128)^2 : 1개
        // m_j1[6(k_e2 + \alpha_e2)]
        mpz_add(temp1, k_e2, alpha_e2);                 mpz_mod(temp1, temp1, n);
        mpz_mul_ui(temp1, temp1, 6);                    mpz_mod(temp1, temp1, n);
        mpz_mul(temp1, temp1, m_last[it][0]);           mpz_mod(temp1, temp1, n);

        // m_j2[6\beta_e2]
        mpz_mul_ui(temp2, beta_e2, 6);                  mpz_mod(temp2, temp2, n);
        mpz_mul(temp2, temp2, m_last[it][1]);           mpz_mod(temp2, temp2, n);

        // m_j3[-6(pre_b_2 + d_128)\alpha_e3 - 6*k_e2(pre_b_4 + \beta_128,0) - k_e1\beta_e1
        //      - 6*\beta_e2(pre_b_1 + k_128,0) - 6*\alpha_e2(pre_b_4 + \beta_128,0) - \alpha_e1\beta_e1 - 6*\beta_e2(pre_b_3 + \alpha_128,0)]
        mpz_add(temp3, pre_b[round][1], di_table[round]);   mpz_mod(temp3, temp3, n);
        mpz_mul(temp3, temp3, alpha_e3);                    mpz_mod(temp3, temp3, n);
        mpz_mul_si(temp3, temp3, -6);                       mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, k_e2);                                mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);
        
        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_mul(temp4, k_e1, beta_e1);                              mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][0], ki_table[round][0]);        mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, beta_e2);                             mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, alpha_e2);                            mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_mul(temp4, alpha_e1, beta_e1);                          mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][2], alpha_i_table[round][0]);   mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, beta_e2);                             mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                                mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);


        mpz_mul(temp3, temp3, m_last[it][2]);                       mpz_mod(temp3, temp3, n);

        // temp1 + temp2 + temp3    
        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n);
        mpz_add(res, temp1, temp3);                                 mpz_mod(res, res, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][38], res);

        /************************************(예진 수정)**************************************/
        // 12-1) e_128^2 * rho : 1개
        // 6*m_j3\beta_e2
        mpz_mul(temp1, m_last[it][2], beta_e2);                 mpz_mod(temp1, temp1, n);
        mpz_mul_ui(res, temp1, 6);                              mpz_mod(res, res, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][39], res);
        /**********************************(예진 수정 끝)**************************************/

        // 13) e_128 : 1개
        // m_j1[6(k_e1 + \alpha_e1)]
        mpz_add(temp1, k_e1, alpha_e1);                         mpz_mod(temp1, temp1, n);
        mpz_mul_ui(temp1, temp1, 6);                            mpz_mod(temp1, temp1, n);
        mpz_mul(temp1, temp1, m_last[it][0]);                   mpz_mod(temp1, temp1, n);

        // m_j2[6\beta_e1]
        mpz_mul_ui(temp2, beta_e1, 6);                          mpz_mod(temp2, temp2, n);
        mpz_mul(temp2, temp2, m_last[it][1]);                   mpz_mod(temp2, temp2, n);

        // m_j3[ -6*(pre_b_2 + d_128)\alpha_e1 - 6*k_e1(pre_b_4 + \beta_128,0) - 6*\beta_e1*(pre_b_1 + k_128,0)
        //          - 6*\alpha_e1(pre_b_4 + \beta_128,0) - 6*\beta_e1(pre_b_3 + \alpha_128,0)]
        mpz_add(temp3, pre_b[round][1], di_table[round]);       mpz_mod(temp3, temp3, n);
        mpz_mul(temp3, temp3, alpha_e1);                        mpz_mod(temp3, temp3, n);
        mpz_mul_si(temp3, temp3, -6);                           mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, k_e1);                            mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                            mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                           mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][0], ki_table[round][0]);    mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, beta_e1);                         mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                            mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                           mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, alpha_e1);                        mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                            mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                           mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][2], alpha_i_table[round][0]);   mpz_mod(temp4, temp4, n);
        mpz_mul(temp4, temp4, beta_e1);                         mpz_mod(temp4, temp4, n);
        mpz_mul_ui(temp4, temp4, 6);                            mpz_mod(temp4, temp4, n);

        mpz_sub(temp3, temp3, temp4);                           mpz_mod(temp3, temp3, n);

        mpz_mul(temp3, temp3, m_last[it][2]);                   mpz_mod(temp3, temp3, n);

        // temp1 + temp2 + temp3
        mpz_add(temp1, temp1, temp2);                           mpz_mod(temp1, temp1, n);
        mpz_add(res, temp1, temp3);                             mpz_mod(res, temp1, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][40], res);

        /************************************(예진 수정)**************************************/
        // 13-1) e_128 * rho : 1개
        // 6*m_j3*beta_e1
        mpz_mul(temp1, m_last[it][2], beta_e1);                 mpz_mod(temp1, temp1, n);
        mpz_mul_ui(res, temp1, 6);                              mpz_mod(res, res, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][41], res);

        // 13-2) rho : 1개
        mpz_add(temp1, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp1, temp1, n);
        mpz_mul(res, m_last[it][2], temp1);                         mpz_mod(res, res, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][42], res);
        /**********************************(예진 수정 끝)**************************************/

        // 14) constant : 1개
        mpz_add(temp1, ki_table[round][0], pre_b[round][0]);        mpz_mod(temp1, temp1, n);
        mpz_add(temp2, alpha_i_table[round][0], pre_b[round][2]);   mpz_mod(temp2, temp2, n);
        /*(been)수정*/
        mpz_sub(temp2, temp2, b_last[0]);                           mpz_mod(temp2, temp2, n);
        /************/
        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n);
        mpz_mul(temp1, temp1, m_last[it][0]);                       mpz_mod(temp1, temp1, n);
        
        mpz_add(temp2, pre_b[round][1], di_table[round]);           mpz_mod(temp2, temp2, n);
        /*(been)수정*/
        mpz_sub(temp3, beta_i_table[round][0],b_last[1]);          mpz_mod(temp3, temp3, n);
        /************/
        mpz_add(temp3, temp3, pre_b[round][3]);                     mpz_mod(temp3, temp3, n);
        mpz_add(temp2, temp2, temp3);                               mpz_mod(temp2, temp2, n);
        mpz_mul(temp2, temp2, m_last[it][1]);                       mpz_mod(temp2, temp2, n);

        mpz_add(temp3, pre_b[round][1], di_table[round]);           mpz_mod(temp3, temp3, n);
        mpz_neg(temp3, temp3);                                      mpz_mod(temp3, temp3, n);
        mpz_add(temp4, pre_b[round][2], alpha_i_table[round][0]);   mpz_mod(temp4, temp4, n);
        mpz_mul(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);

        mpz_add(temp4, pre_b[round][2], alpha_i_table[round][0]);   mpz_mod(temp4, temp4, n);
        mpz_neg(temp4, temp4);
        mpz_add(temp5, pre_b[round][0], ki_table[round][0]);        mpz_mod(temp5, temp5, n);   // 예진 수정
        mpz_neg(temp5, temp5);                                     
        mpz_add(temp4, temp4, temp5);                               mpz_mod(temp4, temp4, n);

        mpz_add(temp5, pre_b[round][3], beta_i_table[round][0]);    mpz_mod(temp5, temp5, n);
        mpz_mul(temp4, temp4, temp5);                               mpz_mod(temp4, temp4, n);

        mpz_add(temp3, temp3, temp4);                               mpz_mod(temp3, temp3, n);
        mpz_sub(temp3, temp3, b_last[2]);                           mpz_mod(temp3, temp3, n);

        mpz_mul(temp3, temp3, m_last[it][2]);                       mpz_mod(temp3, temp3, n);

        // temp1 + temp2 + temp3
        mpz_add(temp1, temp1, temp2);                               mpz_mod(temp1, temp1, n);
        mpz_add(res, temp1, temp3);                                 mpz_mod(res, res, n);

        mpz_set(eq_coeff_table->coeff_m_last[it][43], res);
    }
    //==================================================================================================
    // 중간 디버깅 출력
    //==================================================================================================
    // for( int i=0; i<127; i++ ) {
    //     printf("round = %d\n", i);
    //     for( int j=0; j<4; j++ ) {
    //         for( int k=0; k<12; k++ ) {
    //             gmp_printf("coeff[%d][%d][%d] = %ZX\n", i, j, k, eq_coeff_table->coeff_m[i][j][k]);
    //         }
    //         printf("\n");
    //     }
    //     printf("==============\n");
    // }
    //==================================================================================================
    mpz_clear(temp1); mpz_clear(temp2); mpz_clear(temp3); mpz_clear(temp4); mpz_clear(temp5); mpz_clear(temp6); mpz_clear(temp7); mpz_clear(res);
    mpz_clear(k_e1); mpz_clear(k_e2); mpz_clear(k_e3); mpz_clear(alpha_e1); mpz_clear(alpha_e2); mpz_clear(alpha_e3); mpz_clear(beta_e1); mpz_clear(beta_e2); mpz_clear(beta_e3);
}


void write_eq_coeff_table(const char* filename, struct EQ_COEFF_TABLE *eq_coeff_table)
{
    FILE* file = fopen(filename, "r+");
    if( file == NULL ) {
        perror("파일 열기 실패\n");
        exit(EXIT_FAILURE);
    }

    char line[256];
    long insert_pos = -1;
    while( fgets(line, sizeof(line), file) ) {
        if( strstr(line, "const char* coeff_m_str[127][4][12] =") ) {
            insert_pos = ftell(file);
            break;
        }
    }

    if( insert_pos == -1 ) {
        fprintf(stderr, "coeff_m_str 선언을 파일에서 찾을 수 없음.\n");
        fclose(file);
        exit(EXIT_FAILURE);
    }

    fseek(file, insert_pos, SEEK_SET);

    fprintf(file, "{\n");
    for( int i=0; i<127; i++ ) {
        fprintf(file, "\t{\n");
        for( int j=0; j<4; j++ ) {
            fprintf(file, "\t\t{ ");
            for( int k=0; k<12; k++ ) {
                char* str1 = mpz_get_str(NULL, 16, eq_coeff_table->coeff_m[i][j][k]);
                fprintf(file, "\"%064s\", ", str1);
                free(str1);
            }
            fprintf(file, " }, \n");
        }
        fprintf(file, "\t}, \n");
    }
    fprintf(file, "};\n");

    fprintf(file, "const char* coeff_m_last_str[3][44] =\n{\n");
    for( int i=0; i<3; i++ ) {
        fprintf(file, "\t{ ");
        for( int j=0; j<44; j++ ) {
            char* str1 = mpz_get_str(NULL, 16, eq_coeff_table->coeff_m_last[i][j]);
            fprintf(file, "\"%064s\", ", str1);
            free(str1);
        }
        fprintf(file, " }, \n");
    }
    fprintf(file, "};\n");

    fclose(file);
}
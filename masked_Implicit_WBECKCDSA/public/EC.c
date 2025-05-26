#include "EC.h"

void point_init(Point P)  { mpz_inits( P->x, P->y, NULL); }
void point_clear(Point P) { mpz_clears(P->x, P->y, NULL); }


void point_init_set_str(Point P, const char* x_str, const char* y_str, int base)
{
    mpz_init_set_str(P->x, x_str, base);
    mpz_init_set_str(P->y, y_str, base);
}

void point_init_infinity(Point P)
{
    mpz_init_set_ui(P->x, 0);
    mpz_init_set_ui(P->y, 0);
}

int point_is_infinity(Point P)
{
  return (mpz_cmp_ui(P->x, 0) == 0) && (mpz_cmp_ui(P->y, 0) == 0);         // 두 값을 비교하고 왼쪽이 크면 양수, 같으면 0, 오른쪽이 크면 음수
}

int point_equal(Point P, Point Q)
{
  return (mpz_cmp(P->x, Q->x) == 0) && (mpz_cmp(P->y, Q->y) == 0);         // 두 값을 비교하고 왼쪽이 크면 양수, 같으면 0, 오른쪽이 크면 음수
}

int point_is_inverse(Point P, Point Q)
{
  int comp = mpz_cmp(P->x, Q->x) == 0;      // P랑 Q의 X좌표가 같을 때 true로 comp가 1, 다를때는 return comp = 0 , 즉 x좌표 같을때만 검사
  if (comp != 1) {
    return comp;
  }

  // compute negative
  mpz_t Q_y_neg;
  mpz_init(Q_y_neg);
  mpz_neg(Q_y_neg, Q->y);                   // 오른쪽 값에 음수취하고 왼쪽에 대입

  comp = mpz_cmp(P->y, Q_y_neg) == 0;       // P(x,y) = Q(x,-y)  일때만 참(1) return. 아니면 거짓(0) return
  mpz_clear(Q_y_neg);

  return comp;
}

void point_set(Point R, Point P)
{
    mpz_set(R->x, P->x);
    mpz_set(R->y, P->y);
}

void point_add(Point R, Point P, Point Q, mpz_t a, mpz_t p)
{

    if( point_is_infinity(P) ) {        /* Case 2) P: infi -> R = Q */
        // puts("Case) P : infinity");
        point_set(R, Q);
        return;
    } else if( point_is_infinity(Q) ) { /* Case 3) Q: infi -> R = P */
        // puts("Case) Q : infinity");
        point_set(R, P);
        return;
    }

/*
    Case: [P: !infinity] and [Q: !infinity] and Q = -P
        P + Q = P + (-P) = infi
*/
    if( point_is_inverse(P, Q) ) {      /* Q = -P -> R = infi */
        // puts("Case) Q = -P");
        point_init_infinity(R);
        return;
    }

/*
    Case: [P: !infinity] and [Q: !infinity] and Q != -P
        1) (Doubling: if P  = Q)
            lambda = ( 3*(P_x^2) + a )( (2*(P_y))^(-1) )
            
        2) (Addition: if P != Q)
            lambda = ( Q_y - P_y )( (Q_x - P_x)^(-1) )

        R = (R_x, R_y)
        R_x : lambda^2 - P_x - Q_x  mod p
        R_y : (P_x - R_x) * lambda - P_y  mod p
*/
    mpz_t lambda, denominator;
    mpz_inits(lambda, denominator, NULL);

    /* Case 1) Doubling */
    if( P == Q || point_equal(P, Q) ) 
    {   
        /* lambda = 3*(P_x^2) + a */
        mpz_powm_ui(lambda, P->x, 2, p); 
        mpz_mul_ui(lambda, lambda, 3);   
        mpz_add(lambda, lambda, a);     

        /* denominator = (2*(P_y))^(-1) */
        mpz_mul_ui(denominator, P->y, 2);
        mpz_invert(denominator, denominator, p);
    } 

    /* Case 2) Addition */
    else {
        /* lambda = Q_y - P_y */
        mpz_sub(lambda, Q->y, P->y);

        /* denominator = (Q_x - P_x)^(-1) */
        mpz_sub(denominator, Q->x, P->x);
        mpz_invert(denominator, denominator, p);
    }
    /* 
        Case
        - Doubling: lambda = ( 3*(P_x^2) + a )( (2*(P_y))^(-1) )
        - Addition: lambda = ( Q_y - P_y )( (Q_x - P_x)^(-1) )
    */
    mpz_mul(lambda, lambda, denominator);
    mpz_mod(lambda, lambda, p);

    /* R_x = lambda^2 - P_x - Q_x  mod p */
    mpz_powm_ui(R->x, lambda, 2, p);

    mpz_sub(R->x, R->x, P->x);
    mpz_sub(R->x, R->x, Q->x);
    mpz_mod(R->x, R->x, p);

    /* R_y : (P_x - R_x) * lambda - P_y  mod p */
    mpz_sub(R->y, P->x, R->x);
    mpz_mul(R->y, lambda, R->y);
    mpz_mod(R->y, R->y, p);
    mpz_sub(R->y, R->y, P->y);
    mpz_mod(R->y, R->y, p);

    // clear mpz
    mpz_clears(lambda, denominator, NULL);       // 사용끝난 변수 정리
}

void point_scalar(Point R, Point P, mpz_t scalar, mp_bitcnt_t num_bits, mpz_t a, mpz_t p)
{
    Point tmp;
    point_init(tmp);
    
    for (mp_bitcnt_t i = num_bits - 1; i >= 0 && i < num_bits; i--) /* scalar의 각 비트를 확인하며 이진 스칼라 곱 */
    {
        point_add(tmp, R, R, a, p); /* R <- 2R */
        
        if (mpz_tstbit(scalar, i) == 1) {
            point_add(R, tmp, P, a, p);  /* 비트가 1이면 tmp + P를 계산하여 R에 저장 (R = tmp + P) */
        } else {
            point_set(R, tmp);           /* 비트가 0이면 R에 tmp 값을 설정 (R = tmp) */
        }
    }

    point_clear(tmp);
}

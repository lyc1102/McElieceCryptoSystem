#include <iostream>
using namespace std;

typedef struct Polynomial{
    unsigned degree;
    unsigned coefficient[100];
}Poly_t;

typedef struct Field{
    unsigned max_ele;
    unsigned P_x;
    unsigned gen[1024];
    unsigned gen_inv[1024];
}Field_t;


unsigned Euclidean_add_c(unsigned oprand_1, unsigned oprand_2);

Poly_t Euclidean_add_p(Poly_t oprand_1, Poly_t oprand_2, const Field_t GF);

unsigned Euclidean_mult_cc(unsigned oprand_1, unsigned oprand_2, const Field_t GF);

Poly_t Euclidean_mult_pc(Poly_t oprand_1, unsigned oprand_2, const Field_t GF);

Poly_t Euclidean_mult_pz(Poly_t oprand, unsigned time, const Field_t GF);

Poly_t Euclidean_mult_pp(Poly_t oprand_1, Poly_t oprand_2, const Field_t GF);

Poly_t Euclidean_div_pp(Poly_t divider, Poly_t divisor, const Field_t GF);

Poly_t Euclidean_modp(Poly_t oprand, Poly_t modulo, const Field_t GF);

Poly_t Euclidean_pow(Poly_t oprand, Poly_t modulo, unsigned power, const Field_t GF);

Poly_t Euclidean_inv(Poly_t oprand, Poly_t modulo, const Field_t GF);

Poly_t Euclidean_gcd(Poly_t oprand, Poly_t modulo, const Field_t GF);
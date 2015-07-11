#include "Euclidean.h"
unsigned Euclidean_add_c(unsigned oprand_1, unsigned oprand_2){
    return oprand_1^oprand_2;
}

Poly_t Euclidean_add_p(Poly_t oprand_1, Poly_t oprand_2, const Field_t GF){
    if(oprand_1.degree > oprand_2.degree){
        for(unsigned i = oprand_2.degree + 1; i < oprand_1.degree + 1; i++){
            oprand_2.coefficient[i] = 0;
        }
        oprand_2.degree = oprand_1.degree;
    }else{
        for(unsigned i = oprand_1.degree + 1; i < oprand_2.degree + 1; i++){
            oprand_1.coefficient[i] = 0;
        }
        oprand_1.degree = oprand_2.degree;
    }

    Poly_t result;
    result.degree = oprand_1.degree;
    for(unsigned i = 0; i < result.degree + 1; i++){
        result.coefficient[i] = Euclidean_add_c(oprand_1.coefficient[i], oprand_2.coefficient[i]);
    }
    while((result.coefficient[result.degree] == 0)&&(result.degree > 0)){
        result.degree--;
    }
    return result;
}

unsigned Euclidean_mult_cc(unsigned oprand_1, unsigned oprand_2, const Field_t GF){
    unsigned result;
    if((oprand_1%GF.max_ele == 0)||(oprand_2%GF.max_ele == 0)){
        result = 0;
    }else{
        result = GF.gen[((GF.gen_inv[oprand_1%GF.max_ele] - 1) + (GF.gen_inv[oprand_2%GF.max_ele] - 1))%(GF.max_ele - 1) + 1];
    }
    return result;
}

Poly_t Euclidean_mult_pc(Poly_t oprand_1, unsigned oprand_2, const Field_t GF){
    Poly_t result;
    result.degree = oprand_1.degree;
    for(unsigned i = 0; i < result.degree + 1; i++){
        result.coefficient[i] = Euclidean_mult_cc(oprand_1.coefficient[i], oprand_2, GF);
    }
    if(oprand_2 == 0){
        result.degree = 0;
    }
    return result;
}

Poly_t Euclidean_mult_pz(Poly_t oprand, unsigned time, const Field_t GF){
    Poly_t result;
    result.degree = oprand.degree + time;
    for(unsigned i = 0; i < time; i++){
        result.coefficient[i] = 0;
    }
    for(unsigned i = time; i < result.degree + 1; i++){
        result.coefficient[i] = oprand.coefficient[i - time];
    }
    return result;
}

Poly_t Euclidean_mult_pp(Poly_t oprand_1, Poly_t oprand_2, const Field_t GF){
    Poly_t result;
    result.degree = oprand_1.degree + oprand_2.degree;
    for(unsigned i = 0; i < result.degree + 1; i++){
        result.coefficient[i] = 0;
    }

    for(unsigned i = 0; i < oprand_2.degree + 1; i++){
        result = Euclidean_add_p(result, Euclidean_mult_pz(Euclidean_mult_pc(oprand_1, oprand_2.coefficient[i], GF), i, GF), GF);
    }
    return result;
}

Poly_t Euclidean_div_pp(Poly_t divider, Poly_t divisor, const Field_t GF){
    Poly_t result;
    if(divisor.degree == 0){
        result.degree = divider.degree;
        for(unsigned i = 0; i < result.degree + 1; i++){
            if(divider.coefficient[i] != 0){
                result.coefficient[i] = GF.gen[(GF.gen_inv[divider.coefficient[i]%GF.max_ele] - GF.gen_inv[divisor.coefficient[0]%GF.max_ele] + (GF.max_ele - 1))%(GF.max_ele - 1) + 1];
            }
        }
    }else if(divider.degree < divisor.degree){
        result.degree = 0;
        result.coefficient[0] = 0;
    }else{
        result.degree = divider.degree - divisor.degree;
        for(unsigned i = 0; i < result.degree + 1; i++){
            result.coefficient[i] = 0;
        }
        while((divider.degree >= divisor.degree)&&(divider.degree > 0)){
            result.coefficient[divider.degree - divisor.degree] = GF.gen[(GF.gen_inv[divider.coefficient[divider.degree]%GF.max_ele] - GF.gen_inv[divisor.coefficient[divisor.degree]%GF.max_ele] + (GF.max_ele - 1))%(GF.max_ele - 1) + 1];
            divider = Euclidean_add_p(divider, Euclidean_mult_pz(Euclidean_mult_pc(divisor, result.coefficient[divider.degree - divisor.degree], GF), divider.degree - divisor.degree, GF), GF);
        }
    }
    return result;
}

Poly_t Euclidean_modp(Poly_t oprand, Poly_t modulo, const Field_t GF){
    return Euclidean_add_p(oprand, Euclidean_mult_pp(Euclidean_div_pp(oprand, modulo, GF), modulo, GF), GF);
}

Poly_t Euclidean_pow(Poly_t oprand, Poly_t modulo, unsigned power, const Field_t GF){
    Poly_t result;
    result.degree = 0;
    result.coefficient[0] = 1;
    for(unsigned i = 0; i < power; i++){
        result = Euclidean_modp(Euclidean_mult_pp(result, oprand, GF), modulo, GF);
    }
    return result;
}

Poly_t Euclidean_inv(Poly_t oprand, Poly_t modulo, const Field_t GF){
    Poly_t r_1, r_0, u_1, u_0, v_1, v_0;
    u_1.degree = 0;
    u_1.coefficient[0] = 1;
    u_0.degree = 0;
    u_0.coefficient[0] = 0;
    v_1.degree = 0;
    v_1.coefficient[0] = 0;
    v_0.degree = 0;
    v_0.coefficient[0] = 1;
    if(oprand.degree >= modulo.degree){
        r_1 = oprand;
        r_0 = modulo;
    }else{
        r_1 = modulo;
        r_0 = oprand;
    }
    while((r_0.degree != 0)||(r_0.coefficient[0] != 0)){
        Poly_t quotient = Euclidean_div_pp(r_1, r_0, GF);
        Poly_t remainder = Euclidean_add_p(r_1, Euclidean_mult_pp(r_0, quotient, GF), GF);
        Poly_t temp_u = Euclidean_add_p(u_1, Euclidean_mult_pp(u_0, quotient, GF), GF);
        Poly_t temp_v = Euclidean_add_p(v_1, Euclidean_mult_pp(v_0, quotient, GF), GF);
        r_1 = r_0;
        r_0 = remainder;
        u_1 = u_0;
        u_0 = temp_u;
        v_1 = v_0;
        v_0 = temp_v;
    }
    Poly_t result;
    if(oprand.degree >= modulo.degree){
        result = u_1;
    }else{
        result = v_1;
    }
    result = Euclidean_mult_pc(result, GF.gen[GF.max_ele - GF.gen_inv[r_1.coefficient[0]] + 1], GF);
    return result;
}

#include <string>

//P_X from degree7 down to 0.
//const unsigned int P_x=0x0000011b;

//This returns a*b in GF(2^8)
unsigned int galois_mul(unsigned int a, unsigned int b, unsigned int P_x);

//This gives you the inverse in GF(2^8);
unsigned int inverse(unsigned int input, unsigned int P_x);

//this returns input mod P(x)
unsigned int modp(unsigned int input, unsigned int P_x);

//this returns input^power mod P(x)
unsigned int galois_pow(unsigned input, int power, unsigned int P_x);

//this returns log(n)
int highest_degree(unsigned int a);
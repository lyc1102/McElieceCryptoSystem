#include "GF.h"
#include <iostream>
using namespace std;
int highest_degree(unsigned int a)
{
	int i=63;
	for ( ; i >= 0; --i)
	{
		if((a>>i&1)==1) return i;
	}

	return i;
}

static unsigned int galois_div(unsigned int dividor, unsigned int divisor)
{

	int highest_a=highest_degree(dividor), highest_b=highest_degree(divisor);
	unsigned int quotient=0;
	int diff=0, loop=0;
	while(highest_a>=highest_b)
	{
		loop++;
		if (loop>30)
		{
			break;
		}
		diff=highest_a-highest_b;
		quotient^=(0x00000001<<diff);
		dividor^=(divisor<<diff);
		highest_a=highest_degree(dividor);
	}

	return quotient;
}


static unsigned int mul(unsigned int a, unsigned int b)
{
	unsigned int output=0;
	for (int i = 0; i <=31; ++i)
	{
		if(((b>>i)&1)==1) output^=(a<<i);
	}
	return output;
}

unsigned int galois_mul(unsigned int a, unsigned int b, unsigned int P_x)
{
	unsigned int output = mul(a,b);
	return modp(output, P_x);
}

//This gives you the inverse in GF(2^8);
unsigned int inverse(unsigned int input, unsigned int P_x)
{
	unsigned int r0=P_x, r1=input, t0=0,t1=1;
	unsigned int q=0, temp=0;
	while(r1>0)
	{
		q=galois_div(r0,r1);
		temp=r0;
		r0=r1;
		r1=temp^mul(q,r1);

		temp=t0;
		t0=t1;
		t1=temp^mul(q,t1);
	}
	return modp(t0, P_x);
}

unsigned int modp(unsigned int input, unsigned int P_x)
{
	input^=mul(galois_div(input,P_x),P_x);
	return input;
}

unsigned int galois_pow(unsigned input, int power, unsigned int P_x)
{
	unsigned int output = 1;
	for (int i = 0; i < power; ++i)
	{
		output = galois_mul(output,input,P_x);
	}
	return output;
}
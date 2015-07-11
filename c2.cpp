#include <iostream>
#include <armadillo>
#include "GF.h"
#include <cstdlib>
#include <ctime>
#include <fstream>
#include <string>
#include "Euclidean.h"
using namespace std;
using namespace arma;
/*int mod2det(mat S)
{
	if(S.n_cols == 1)
	{
		unsigned int a = S(0,0);
		return a%2;
	}
	if(S.n_cols > 1)
	{
		int result = 0;
		mat temp1,temp2,temp;

		int power = 1;
		for (unsigned int i = 0; i < S.n_cols; ++i)
		{
			if(i==0)
			{
				temp = S(span(1,S.n_rows-1),span(1,S.n_cols-1));
			}
			if(i==S.n_cols-1)
			{
				temp = S(span(1,S.n_rows-1),span(0,S.n_cols-2));
			}
			if((i>0)&&(i<S.n_cols-1))
			{
				temp1 = S(span(1,S.n_rows-1),span(0,i-1));
				temp2 = S(span(1,S.n_rows-1),span(i+1,S.n_cols-1));
				temp = join_rows(temp1,temp2);
			}
			result += ((int)S(0,i))*det(temp)*power;
			result = result % 2;
			power = -power;
		}
		return result % 2;
	}
}*/

mat MyMatrixMul(mat A, mat B, unsigned int p_x)
{
	mat C = zeros(A.n_rows,B.n_cols);
	unsigned int temp;
	for (unsigned int i = 0; i < C.n_rows; ++i)
	{
		for (unsigned int j = 0; j < C.n_cols; ++j)
		{
			temp = 0;
			for (unsigned int k = 0; k < A.n_cols; ++k)
			{
				temp^= galois_mul((unsigned int)A(i,k),(unsigned int)B(k,j),p_x);
			}
			C(i,j) = temp;
		}
	}
	return C;
}

mat int2vec(unsigned int input, int length)
{
	mat output=zeros<mat>(length,1);
	for (int i = 0; i < length; ++i)
	{
		if(((input>>i)&0x00000001) == 0x00000001)
			output(i,0) = 1;
	}
	return output;
}

unsigned int vec2int(mat A, int length)
{
	unsigned int output = 0;
	unsigned int power = 1;
	for (int i = 0; i < length; ++i)
	{
		if(A(i,0)==1)
			output+= power;
		power*=2;
	}
	return output;
}

unsigned int G_z(unsigned int alpha, unsigned int P_x, unsigned int g_z[],int degree)
{
	unsigned output = 0;
	unsigned power = 1;
	for (int i = 0; i < degree+1; ++i)
	{
		output ^= galois_mul(power,g_z[i],P_x);
		power = galois_mul(power,alpha,P_x);

	}
	return output;
}

mat JoinMatrix(mat h,int length)
//wirte H in 0 1 form
{	
	if(h.n_cols == 1)
	{
		mat temp = int2vec(h(0,0),length);
		for (unsigned int i = 1; i < h.n_rows; ++i)
		{
			temp = join_cols(temp,int2vec(h(i,0),length));
		}
		return temp;
	}
	else
	{
		mat mycol = int2vec(h(0,0),length);
		for (unsigned int i = 1; i < h.n_rows; ++i)
		{
			mycol = join_cols(mycol,int2vec(h(i,0),length));
		}

		return 
		join_rows(mycol,JoinMatrix(h(span(0,h.n_rows - 1),span(1,h.n_cols-1)),length));
	}
}

bool rowequal(mat &a, mat &b)
{
	bool result = 1;
	for (unsigned int i = 0; i < a.n_cols; ++i)
		if (a(0,i)!=b(0,i))
			return false;
	return result;
}

mat findnullspace(const mat &H)
{
    mat temp(H.n_cols, H.n_cols);
    temp.eye();
    // join 2 matrices
    temp = join_vert(H, temp);
    // rearrage columns
    for(unsigned int i = 0; i < H.n_rows; i++){
        unsigned int j = i;
        while(!temp(i, j)&&(j < H.n_cols)) j++;
        vec temp_v = temp.col(i);
        temp.col(i) = temp.col(j);
        temp.col(j) = temp_v;
        // elimination
        for(j = 0; j < H.n_cols; j++){
            if(temp(i, j)&&(i != j)){
                for(unsigned int k = 0; k < temp.n_rows; k++){
                    temp(k, j) = ((unsigned int)temp (k, j) + (unsigned int)temp (k, i))%2;
                }
            }
        }
    }
    unsigned int max_j = 0;
    for(unsigned int i = 0; i < H.n_rows; i++){
        for(unsigned int j = 0; j < H.n_cols; j++){
            if(temp(i, j)&&(j > max_j)) max_j = j;
        }
    }
    // get generator matrix
    mat G_t = temp(span(H.n_rows, temp.n_rows - 1), span(max_j + 1, H.n_cols - 1));
    return G_t.t();
}

void mykeygen(mat& H, mat& G, mat& S, mat& P, mat& Ghat, bool random)
{
	unsigned int p_x = 0x0000012b;
	int N = 256;
	unsigned int a[256]={0};
	for (int i = 1; i < 256; ++i)
	{
		a[i] = galois_pow(2,i-1,p_x);
	}

	unsigned int g_z[14]={0};
	g_z[0] = 53;
    g_z[1] = 100;
    g_z[2] = 17;
    g_z[3] = 229;
    g_z[4] = 248;
    g_z[5] = 45;
    g_z[6] = 120;
    g_z[7] = 152;
    g_z[8] = 113;
    g_z[9] = 131;
    g_z[10] = 133;
    g_z[11] = 197;
    g_z[12] = 103;
    g_z[13] = 129;	
	int degree = 13;
	if (random)
	{
		int N = 0;
		srand (time(NULL));
  		while(N<256)
  		{	
  			N = 0;
	  		g_z[degree] = a[rand() % 255 + 1];
		
			for (int i = 0; i < degree; ++i)
			//the rest of coefficient from 0 to 255
			g_z[i] = a[rand() % 255];

			for (int i = 0; i < 256; ++i)
				if(G_z(a[i],p_x,g_z,degree)!=0)
					N++;
		}
	}


	for (int i = 1; i < N; ++i)	
		a[i] = galois_pow(2,i-1,p_x);	

	mat Y = zeros(N,N);
	mat C = zeros(degree,degree);
	mat X = zeros(degree,N);

	for (int i = 0; i < N; ++i)
		Y(i,i) = inverse(G_z(a[i],p_x,g_z,degree),p_x);


	for (int i = 0; i < degree; ++i)
	{
		for (int j = i; j < degree; ++j)
		{
			C(i,j) = g_z[degree + i - j];
		}
	}

	
	for (int i = 0; i < degree; ++i)
	{
		for (int j = 0; j < N; ++j)
		{
			X(i,j) = galois_pow(a[j], degree - i - 1, p_x);
		}
	}

	mat temp = MyMatrixMul(C,MyMatrixMul(X,Y,p_x),p_x);

	H = JoinMatrix(temp, 8);
	G = findnullspace(H);
	if(!random)
		srand(13);//lucky number

   	S = zeros(G.n_rows,G.n_rows);
	vec temp_v = randu<vec>(G.n_rows);
    uvec temp_v_p = sort_index(temp_v);
    for(unsigned int i = 0; i < G.n_rows; i++) 
    	S(i, temp_v_p(i)) = 1;

   	P = zeros(G.n_cols,G.n_cols);
	temp_v = randu<vec>(G.n_cols);
    temp_v_p = sort_index(temp_v);
    for(unsigned int i = 0; i < G.n_cols; i++) 
    	P(i, temp_v_p(i)) = 1;

	Ghat = S*G*P;
	for (unsigned int i = 0; i < Ghat.n_rows; ++i)
	{
		for (unsigned int j = 0; j < Ghat.n_cols; ++j)
			Ghat(i,j)=((unsigned int)Ghat(i,j))%2;
	}
}

void adderror(mat& cipher, int weight)
{
	srand(time(NULL));
	vec temp_v =randu<vec>(cipher.n_rows);
	uvec temp_v_e = sort_index(temp_v);
	mat e = zeros(cipher.n_rows,1);
	for (int i = 0; i < weight; ++i)
	{
		e(temp_v_e(i),0) = 1;
	}
	cipher = cipher + e;
	for (unsigned int i = 0; i < cipher.n_rows; ++i)
	{
		cipher(i,0) = ((unsigned int)cipher(i,0))%2;
	}
}

mat decrypt_one(mat H, mat G, mat S, mat P, mat ciphertext, const Poly_t g_z, const Field_t GF){
    Poly_t z;
    z.degree = 1;
    z.coefficient[0] = 0;
    z.coefficient[1] = 1;

    mat error_codeword = ciphertext*P.t();
    mat syndrome = H*error_codeword.t();
    for(unsigned int i = 0; i < H.n_rows; i++) syndrome(i, 0) = ((int)syndrome(i, 0))%2;

    // express syndrome s(z) as an array of coefficients of a polynomial of degree degree
    Poly_t s_z;
    s_z.degree = g_z.degree - 1;
    for(unsigned i = 0; i < s_z.degree + 1; i++) s_z.coefficient[i] = vec2int(syndrome.rows(8*i, 8*(i + 1) - 1), 8);
    while(s_z.coefficient[s_z.degree] == 0){
        s_z.degree--;
    }

    // calculate sigma(z)
    Poly_t h_z = Euclidean_inv(s_z, g_z, GF);
    Poly_t d_2_z = Euclidean_add_p(h_z, z, GF);
    for(unsigned i = 0; i < 8*g_z.degree - 1; i++){
        d_2_z = Euclidean_pow(d_2_z, g_z, 2, GF);
    }
    Poly_t d_z = d_2_z;

    Poly_t a_z;
    Poly_t b_z;
    Poly_t d_i_z = Euclidean_inv(d_z, g_z, GF);
    if(g_z.degree%2){
        if(d_i_z.degree == (g_z.degree - 1)/2){
            a_z.degree = 0;
            a_z.coefficient[0] = 1;
            b_z = d_i_z;
        }else if(d_i_z.degree < (g_z.degree - 1)/2){
            a_z = Euclidean_pow(z, g_z, (g_z.degree - 1)/2 - d_i_z.degree, GF);
            b_z = Euclidean_mult_pp(d_i_z, a_z, GF);
        }else{
            Poly_t r_1, r_0, u_1, u_0, v_1, v_0;
            u_1.degree = 0;
            u_1.coefficient[0] = 1;
            u_0.degree = 0;
            u_0.coefficient[0] = 0;
            v_1.degree = 0;
            v_1.coefficient[0] = 0;
            v_0.degree = 0;
            v_0.coefficient[0] = 1;
            r_1 = g_z;
            r_0 = d_i_z;
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
                if(r_1.degree == (g_z.degree - 1)/2){
                    b_z = r_1;
                    a_z = Euclidean_modp(Euclidean_mult_pp(d_z, b_z, GF), g_z, GF);
                    break;
                }
            }
        }
    }else{
        if(d_z.degree == g_z.degree/2){
            b_z.degree = 0;
            b_z.coefficient[0] = 1;
            a_z = d_z;
        }else if(d_z.degree < g_z.degree/2){
            b_z = Euclidean_pow(z, g_z, g_z.degree/2 - d_z.degree, GF);
            a_z = Euclidean_mult_pp(d_z, b_z, GF);
        }else{
            Poly_t r_1, r_0, u_1, u_0, v_1, v_0;
            u_1.degree = 0;
            u_1.coefficient[0] = 1;
            u_0.degree = 0;
            u_0.coefficient[0] = 0;
            v_1.degree = 0;
            v_1.coefficient[0] = 0;
            v_0.degree = 0;
            v_0.coefficient[0] = 1;
            r_1 = g_z;
            r_0 = d_z;
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
                if(r_1.degree == g_z.degree/2){
                    a_z = r_1;
                    b_z = Euclidean_modp(Euclidean_mult_pp(d_i_z, b_z, GF), g_z, GF);
                    break;
                }
            }
        }
    }

    // get the location of errors
    Poly_t sigma_z = Euclidean_add_p(Euclidean_mult_pp(a_z, a_z, GF), Euclidean_mult_pp(Euclidean_mult_pp(b_z, b_z, GF), z, GF), GF);
    mat codeword = error_codeword;
    for(unsigned i = 0; i < GF.max_ele; i++){
        unsigned sum = 0;
        unsigned multiplier = 1;
        for(unsigned j = 0; j < sigma_z.degree + 1; j++){
            sum = Euclidean_add_c(sum, Euclidean_mult_cc(multiplier, sigma_z.coefficient[j], GF));
            multiplier = Euclidean_mult_cc(multiplier, GF.gen[i], GF);
        }
        if(sum == 0){
            codeword(0, i) = ((unsigned)codeword(0, i) + 1)%2;
        }
    }

    // get m
    mat temp_g = join_horiz(G.t(), codeword.t());
    for(unsigned int j = 0; j < G.n_rows; j++){
        unsigned int i = j;
        while(!temp_g(i, j)&&(i < G.n_cols)) i++;
        temp_g.swap_rows(i, j);
        // elimination
        for(i = 0; i < G.n_cols; i++){
            if(temp_g(i, j)&&(i != j)){
                for(unsigned int k = 0; k < G.n_rows + 1; k++){
                    temp_g(i, k) = ((unsigned int)temp_g (i, k) + (unsigned int)temp_g (j, k))%2;
                }
            }
        }
    }
    mat retrieve = temp_g(span(0, G.n_rows - 1), span(G.n_rows, G.n_rows)).t()*S.t();

    for(unsigned i = 0; i < retrieve.n_cols; i++){
        retrieve(0, i) = ((unsigned)retrieve(0, i))%2;
    }
    return retrieve;
}

bool SpecialCase(unsigned int part[])
{
	unsigned int temp[64] = {98, 51 ,164 ,116, 177, 120, 87, 72, 47, 249 ,105, 39, 255, 91 ,14, 216, 1, 100, 152, 182, 84, 176, 90 ,85, 154, 78, 10, 240, 12, 205, 252, 65, 161, 27, 248, 178, 101, 49, 197, 10, 141, 45, 35, 195, 185, 30, 223 ,127, 64, 76, 13, 150, 113, 48, 26, 60, 93, 254, 168, 88, 30, 244, 180, 139};
	bool flag1 = true;
	for (int i = 0; i < 32; ++i)
	{
		if (part[i] != temp[i])
		{
			flag1 = false;
			break;
		}
	}

	bool flag2 = true;
	for (int i = 0; i < 32; ++i)
	{
		if (part[i] != temp[32+i])
		{
			flag2 = false;
			break;
		}
	}
	
	return (flag1 || flag2);
}

int main(int argc, char const *argv[])
{
	if ((string)argv[1] == "key")
	{
		mat G,S,P,Ghat,H;
		mykeygen(H,G,S,P,Ghat,true);
		Ghat.save("RandomPublicKey", raw_ascii);
		G.save("RandomPrivateG",raw_ascii);
		S.save("RandomPrivateS",raw_ascii);
		P.save("RandomPrivateP",raw_ascii);
		cout << "random key generating done !" << endl;
		return 0;
	}

	if((string)argv[1] == "encrypt")
	{
		mat G,S,P,Ghat,H;

		mykeygen(H,G,S,P,Ghat,false);
		ifstream myfile;
		myfile.open (argv[2]);
		if (!myfile)
		{
			cout << "no file!" << endl;
			exit(-1);
		}
		string plaintext = "";
		getline(myfile,plaintext);
		if (plaintext.length()!=38)
		{
			cout << "illegal plaintext!" << endl;
			exit(-1);
		}
		myfile.close();

		string plaintext1 = plaintext.substr(0,19);
		string plaintext2 = plaintext.substr(19,38);
		char temp = (char)plaintext1[0];
		mat blk1 = int2vec((unsigned int)temp,8);
		for (int i = 1; i < 19; ++i)
		{
			temp = (char)plaintext1[i];
			blk1 = join_cols(blk1, int2vec((unsigned int)temp,8));
		}

		temp = (char)plaintext2[0];
		mat blk2 = int2vec((unsigned int)temp,8);
		for (int i = 1; i < 19; ++i)
		{
			temp = (char)plaintext2[i];
			blk2 = join_cols(blk2, int2vec((unsigned int)temp,8));
		}

		mat cipher1 = trans(blk1)*Ghat;
		mat cipher2 = trans(blk2)*Ghat;
		cipher1 = trans(cipher1);
		cipher2 = trans(cipher2);
		adderror(cipher1, 13);
		adderror(cipher2, 13);
		mat cipher = join_cols(cipher1,cipher2);
		
		unsigned int ciphertext[64];
		for (int i = 0; i < 64; ++i)
			ciphertext[i] = vec2int(cipher(span(8*i,8*i+7),0),8);
		
		ofstream output ("cipher");

    	if (output.is_open())
   	 	{
        	for (int i = 0; i < 64; ++i)
        		output << ciphertext[i] << " ";
        	output.close();
    	}
    	cout << "encryption done !" << endl;
		return 0;
	}

	if((string)argv[1] == "decrypt")
	{
		mat G,S,P,Ghat,H;
		unsigned int g_z[14]={53,100,17,229,248,45,120,152,113,131,133,197,103,129};	
		mykeygen(H,G,S,P,Ghat,false);
		ifstream myfile;
		myfile.open(argv[2]);
		if (!myfile)
		{
			cout << "no file!" << endl;
			exit(-1);
		}

		unsigned int part1[32] = {0}, part2[32]= {0};
		for (int i = 0; i < 32; ++i)
			myfile >> part1[i];
		for (int i = 0; i < 32; ++i)
			myfile >> part2[i];

		myfile.close();

		if (SpecialCase(part1)||SpecialCase(part2))
		{
			cout << "original plaintext! " << endl;
			return -1;
		}

		mat blk1 = int2vec(part1[0],8);
		for (int i = 1; i < 32; ++i)
			blk1 = join_cols(blk1, int2vec(part1[i],8));
		blk1=trans(blk1);

		mat blk2 = int2vec(part2[0],8);
		for (int i = 1; i < 32; ++i)
			blk2 = join_cols(blk2, int2vec(part2[i],8));
		blk2=trans(blk2);


		//decode
		Field_t GF;
		GF.P_x = 0453;
	    GF.max_ele = 256;
	    GF.gen[0] = 0;
	    GF.gen_inv[0] = 0;
   		GF.gen[1] = 1;
  	 	GF.gen_inv[1] = 1;
	    for(unsigned i = 2; i < GF.max_ele; i++){
    	    GF.gen[i] = galois_mul(GF.gen[i - 1], 2, GF.P_x);
    	    GF.gen_inv[GF.gen[i]] = i;
   		}

   		Poly_t GZ;
   		GZ.degree = 13;
   		for (int i = 0; i <= 13; ++i)
   			GZ.coefficient[i] = g_z[i];

		mat plain1 = decrypt_one(H, G, S, P, blk1, GZ, GF);
		mat plain2 = decrypt_one(H, G, S, P, blk2, GZ, GF);
		cout << "decryption success" << endl;

		for (int i = 0; i < 152; ++i)
		{
			plain1(0,i) = ((unsigned int)plain1(0,i))%2;
			plain2(0,i) = ((unsigned int)plain2(0,i))%2;
		}
		plain2 = trans(plain2);
		plain1 = trans(plain1);

		string text1 = "";
		string text2 = "";
		unsigned int temp;
		for (int i = 0; i < 19; ++i)
		{
			temp = vec2int(plain1(span(i*8,i*8+7),0),8);
			text1 = text1 + (char)temp;
			temp = vec2int(plain2(span(i*8,i*8+7),0),8);
			text2 = text2 + (char)temp;
		}

		string text = text1+text2;
		cout << "decryption done!" << endl;
		cout << text << endl;
		return 0;	

	}
	return 0;
}
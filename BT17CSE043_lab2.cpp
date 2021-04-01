/*

Our scheme makes use of five algorithms: the Setup algorithm, the Keygen algorithm, the SPE_PP algorithm, the
Trapdoor algorithm and the Test algorithm. The scheme
works as follows:

Setup(1^lamda): Given a security parameter lamda, the algorithm
first selects two cyclic groups G1, GT with the same prime
order q and a bilinear pairing e : G1 × G1 ? GT . Then,
the algorithm chooses a generator P of G1 and two hash
functions h1 : G1 -> Z*q
q , h2 : {0, 1}* -> Z*q
q . Finally, the algorithm publishes the public parameter Para =
(G1, GT , q, P, e, h1, h2).

KeyGen(Para): DU and DS execute this algorithm to
generate the public/secret key pairs. P Ku = a P, P Ks = bP
are the public keys of them and the secret keys corresponding
to the public keys are sku = a, sks = b, where a, b belongs to Z*q
are randomly selected by DU and DS, respectively. Finally,
they send the public keys PKu, PKs to CA.

*/

#include <cstring>
#include <fstream>
#include <string>
#include <iostream>
#include <pbc/pbc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <assert.h>
#include <bits/stdc++.h>

# define MAX 100000000
using namespace std;

/*
	struct setup_output is a structure.
	The data type of each variable is explained here. Use of each variable is explained later.
*/
typedef struct setup_output 
{
	//GMP type for integers
 	mpz_t q;    // order of group (r)
 	
	 //pairings where elements belong; can initialize from sample pairing parameters bundled with PBC in the param subdirectory
    pairing_t pairing;
	
	//used to generate pairing parameters	
    pbc_param_t par;
    
    //elements of an algebraic structure
    element_t g1,g2,gt;	//elements of group G1, G1 and GT
    element_t P;    // Generator of group G1

}setup_result;

/*
	struct Keys is a structure.
	The data type of each variable is explained here. Use of each variable is explained later.
*/
typedef struct Keys 
{
    element_t PKu,PKs;	//PKu is the public key of data user of type element_t and PKs is the public key of data sender of type element_t
    mpz_t SKu, SKs;		//SKu is the public key of data user of type mpz and SKs is the public key of data sender of type mpz

}keys;

setup_result globle_setup;	//The global_setup is a variable of (type - setup_output) that holds values of the variables set by setup function to be used by other functions.

keys MyKeys;	//The MyKeys is a variable of (type - Keys) that holds values of the public and private keys to be used by other functions.


//function to convert a string to binary string 
string strToBinary(string s)
{
    int n = s.length();	//get length of string
  
    string result="";
    for (int i = 0; i <= n; i++)
    {
        // convert each char to ASCII value
        int val = int(s[i]);
  		// Convert ASCII value to binary
        string bin = "";
        while (val > 0)
        {
        	//if value is odd pust 1 in bin else push 0
            (val % 2)? bin.push_back('1') : bin.push_back('0');
            val /= 2;
        }
        reverse(bin.begin(), bin.end());
        result+=bin;	//append the result to bin
    
    }
    return result;	//return the converted binary string
}


//function to generate keys and store then in global variable MyKeys
void KeyGen()
{
    mpz_t a,b;	//a and b are two type mpz numbers
    //Initilizing a
    mpz_init(a);
    //Initilizing b
    mpz_init(b);
    // Setting them to random values
    mpz_random(a,MAX);
    mpz_random(b,MAX);
    
    //Printing the value of order of group
    //gmp_printf("q = %Zd \n", globle_setup.q);
    
    cout<<endl<<"==============================================================="<<endl;
    cout<<"Key Generation Algorithm"<<endl;
    cout<<"==============================================================="<<endl<<endl;
    //loop until we get a positive value for a and b
    while(mpz_cmp_ui(a, 0U) < 0)
    {
        mpz_random(a,MAX);
        mpz_mod(a, a, globle_setup.q);
    }
    while(mpz_cmp_ui(b, 0U) < 0)
    {
        mpz_random(b,MAX);
        mpz_mod(b, b, globle_setup.q);
    }
    
    //Taking mod so that they belong to Z*q and storing them in global_setup variable
    mpz_mod(a, a, globle_setup.q);
    mpz_mod(b, b, globle_setup.q);
    
    //Printing the values of a and b
    gmp_printf("(SKu) Data user secret key: %Zd \n", a);
    gmp_printf("(Sks) Data sender secret key: %Zd \n", b);
    
    //Initializing the values of MyKeys elements PKu and PKs 
    element_init_G1(MyKeys.PKu, globle_setup.pairing);
    element_init_G1(MyKeys.PKs, globle_setup.pairing);  
	
	//calculating and storing  public keys as PKu = aP and  PKs = bP where P is the generator og group G1   
    element_mul_mpz(MyKeys.PKu, globle_setup.P, a);
    element_mul_mpz(MyKeys.PKs, globle_setup.P, b);
    
    //Printing the values of Public key of data user and data sender
    element_printf("\n(PKu) Data user public key: %B", MyKeys.PKu);
    element_printf("\n(PKs) Data sender public key: %B\n", MyKeys.PKs);
    
   	//Initializing the values of MyKeys elements Sku and SKs 
    mpz_init(MyKeys.SKu);
    mpz_init(MyKeys.SKs);
    
     //Setting the values of secret keys of Data user and Data sender SKu = a and SKs = b
    mpz_set(MyKeys.SKu, a);
    mpz_set(MyKeys.SKs, b);
    
    cout<<endl;

}

//Hash function h1 : G1 -> Z*q
void hash1(element_t e,mpz_t h1_val)
{
	//Declaring variable data of type char *
    unsigned char data[100];
    //Declaring variable s of type element_t
    element_t s;
    //Initializing variable s of type element_t to pairing
    element_init_Zr(s, globle_setup.pairing);

	//storing char * form of e (element of group G1) sent using function element_to_bytes
    element_to_bytes(data, e);
    //Generate an element e deterministically from the len bytes stored in the buffer data  string ->element
    element_from_hash(s,data,20);
    //Converts e to a GMP integer z if such an operation makes sense
    element_to_mpz(h1_val,s);
 
}

//Hash function h2 : {0,1}* -> Z*q
void hash2(string str ,mpz_t h2_val)
{
    //Declaring variable s of type element_t
    element_t s;
     int len = str.length();
    unsigned char result[len];
    for(int i=0;i<len;i++)
    {
        result[i]=str[i];
    }
    //Initializing variable s of type element_t to pairing
    element_init_Zr(s, globle_setup.pairing);
   	//Generate an element e deterministically from the len bytes stored in the buffer data
    element_from_hash(s,result,20);
    //Converts e to a GMP integer z if such an operation makes sense
    element_to_mpz(h2_val, s); 

}

void setup(mpz_t security_parameter) 
{
	cout<<endl<<"==============================================================="<<endl;
    cout<<"Setup Algorithm"<<endl;
    cout<<"==============================================================="<<endl;
    pairing_t pairing;	//pairings where elements belong; can initialize from sample pairing parameters bundled with PBC in the param subdirectory
    pbc_param_t par;
    
    
    //=========================================type a curve starts here=================================================================
    mpz_t rb; 
    mpz_init(rb);
    
    /*
    	type a curve:(y^2 = x^3 + x)
		exp2, exp1, sign1, sign0, r:
  		r = 2^exp2 + sign1 * 2^exp1 + sign0 * 1 (Solinas prime)
		q, h:
  		r * h = q + 1
  		q is a prime, h is a multiple of 12 (thus q = -1 mod 12)
	*/
    int rbits=mpz_get_ui(security_parameter)+1;	//Here value of rbits is set to value one more than that of security_paramenter which is the bits of order of group
    mpz_set_ui(rb,rbits);
    int qbits=10;	//Value of q bits is set to 10
    pbc_param_init_a_gen(globle_setup.par,rbits,qbits);  // Initializing A type curve
    
	//pairing_init_pbc_param: Initialize a pairing with pairing parameters p
    pairing_init_pbc_param(globle_setup.pairing, globle_setup.par);
    pairing_init_pbc_param(pairing, globle_setup.par);
    cout<<endl<<"Curve paramenters: "<<endl<<endl;
    pbc_param_out_str(stdout, globle_setup.par);    // Printing the A type curve parameters
    
	
	//writing the values to a file
    FILE *stream;
    stream = fopen("a_param.txt", "w+");    
    pbc_param_out_str(stream, globle_setup.par);
    fclose(stream);
  
  	//reading values from file
    FILE *reading;
    char buff[1024];

    reading = fopen("a_param.txt", "r");
    fscanf(reading, "%s", buff);
    fgets(buff, 1024, (FILE*)reading);
    fgets(buff, 1024, (FILE*)reading);
    fgets(buff, 1024, (FILE*)reading);
    fgets(buff, 1024, (FILE*)reading);
    
	char * r= buff+2;
    mpz_init(globle_setup.q);
    mpz_init_set_str(globle_setup.q,r,10);
    
    //=========================================================type a curve ends here ================================================
    
    
    
    //========================================type a1 curve starts here=================================================================
    
    
    /*
    	type a1 curve: (y^2 = x^3 + x)
    
   	 	p, n, l:
		p + 1 = n * l
		p is prime, same as the q in a_param, n is the order of the group
	*/
	/*
    mpz_t k; 
	mpz_init(k);
	mpz_set(k, security_parameter);		// storing security parameter in k var
	mpz_t q; 
	mpz_init(q);
	random_prime_bits(q, k);			// Generating a random prime number of k bits
	gmp_printf("\nOrder of group G1, GT is q: %Zd\n", q);		// printing order of elliptic curve

	
	pairing_t pairing; 
	pbc_param_t par;
	
	pbc_param_init_a1_gen(par, q);		// Setting parameter "n",of a1 type elliptic curve 
	
	pairing_init_pbc_param(pairing, par);	//Initialize a pairing with pairing parameters par
	printf("\nParameters are (p + 1 = n * l where p is prime, same as the q in a_param, n is the order of the group): \n");
	pbc_param_out_str(stdout, par);			// Printing all parameters
	mpz_init(globle_setup.q);
    mpz_set(globle_setup.q,q);
    */
    
    //printing the order of group
    gmp_printf("\nOrder of group is: %Zd\n\n",globle_setup.q);
    
    //Declaring elements of group G1, G1 and GT
    element_t g1, g2, gt, p;
	
	//element is initialized it is associated with an algebraic structure
    element_init_G1(g1, pairing);
    element_init_G1(g2, pairing);
   	element_init_GT(gt, pairing);
    
    //Random elements are choosen to represent group as order of group is odd so every element of group is a generator
	element_random(g1);
    element_random(g2);

	//element is initialized it is associated with an algebraic structure
    element_init_G1(globle_setup.g1, pairing);
    element_init_G1(globle_setup.g2, pairing);

	//values asssigned to global variables
    element_set(globle_setup.g1, g1);
    element_set(globle_setup.g2, g2);
    
    //Values of g1, g2 are printed
    element_printf("Element of G1 group g1: %B\n", g1);
	element_printf("Element of G1 group g2: %B\n", g2);

	//Computes a pairing: out = e(in1, in2), where in1, in2, out must be in the groups G1, G2, GT.
	element_pairing(gt,g1,g2);

	//element global_setup.gt is initialized it is associated with an algebraic structure GT
    element_init_GT(globle_setup.gt,pairing);
    //element global_setup.gt is set to value of gt
    element_set(globle_setup.gt,gt);

	//values of pairing are printed
    element_printf("Applying bilinear pairing on g1 and g2, gt: %B\n", gt);
    
    //Generator is choosen as random element from group g1 and global_setup.P is assigned 
    element_init_G1(globle_setup.P,pairing);
    element_init_G1(p, pairing);
    element_random(p);
    element_set(globle_setup.P,p);
    
    //value of selected generator is printed
    element_printf("\nGenerator selected: %B\n", p);
    
    //Hashing
    
    //h1_val and h2_val are used to store results
    mpz_t h1_val,h2_val;
    //initizlizing
    mpz_init(h1_val);
    mpz_init(h2_val);
    
    //Hash 1
    hash1(p, h1_val);	//h1 : G1 -> Z*q
    
    //Hash 2
    string msg = "HelloWorld";	//msg whose hash value is to be calculated
	string bin = strToBinary(msg);	//message converted to binary

    hash2(bin,h2_val);	//h2 : {0, 1}* -> Z*q
    
    //Values of hash are printed 
    cout<<endl<<"Hash 1: ";
    element_printf("Element of group G1: %B -> ", p);
    gmp_printf("%Zd\n",h1_val);
    
    cout<<"Hash 2: (message) "<<msg<<endl<<bin<<" -> ";
    cout<<" -> ";
    gmp_printf("%Zd\n",h2_val);
    
}

int main () 
{
	// security paramater of type mpz
    mpz_t security_parameter;
    //Initializing security para,eter
    mpz_init(security_parameter);
    mpz_set_ui(security_parameter, 10);      // Setting lambda =6
    
    //Setup Algorithm
	setup(security_parameter);
    
    //Key Generation Algorithm
	KeyGen();
    
	return 0;
}
/*
//converting data to string and storing it in variable message
    string message="";
    for(int i=0;i<100;i++)
    {
        message+=data[i];
    }
    //Applying sha256 algorithm to calculate hash_value
    string hash_value = sha256(message); //string ->string
    int len = hash_value.length();
    unsigned char result[len];
    //Converting hash_value to char * array variable result
    for(int i=0;i<len;i++)
    {
        result[i]=hash_value[i];
    }*/

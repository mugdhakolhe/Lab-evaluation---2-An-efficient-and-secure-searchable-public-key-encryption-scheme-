#include <pbc/pbc.h>
#include <stdio.h>
#include <string.h>
#include <stdlib.h>
#include <gmp.h>
#include <time.h>
#include <assert.h>


//function to generate prime of n-bits 
void random_prime_bits(mpz_t result, mpz_t n) 
{
	mpz_t bits;
	mpz_init(bits);
	gmp_randstate_t state;
	gmp_randinit_default(state);		// Initialize state with default random number generator algorithm.
	gmp_randseed_ui(state, (rand()+1)*(rand()+1));		// Seed random number generator.
	if (mpz_cmp_ui(n,1) <= 0) 						//  If size is less than equal to 1 
	{
		printf("NO PRIME EXISTS\n");
	} 
	else 
	{
		mpz_t lower_limit;
		mpz_init(lower_limit);
		mpz_ui_pow_ui(lower_limit, 2, mpz_get_ui(n)-1);		// lower limit = 2^(n-1)
		while (1)
		{
			mpz_urandomb(bits, state ,mpz_get_ui(n));
			if (mpz_cmp(bits, lower_limit) > 0 && mpz_probab_prime_p(bits,mpz_get_ui(n))) 
			{
				mpz_set(result,bits);	// Break when the generated number is prime and greater than 2^(n-1)
				break;
			}
		}			
	}

}
char* stringToBinary(char* s) 		// Converts string to binary
{
    if(s == NULL) return 0;
    printf("%s\n", s); 
    size_t len = strlen(s);		
    char *binary = malloc(len*8 + 1); 
	binary[0] = '\0';
    for(size_t i = 0; i < len; ++i) 
	{
        char ch = s[i];
        for(int j = 7; j >= 0; --j)
		{
            if(ch & (1 << j)) 		//shift 1 left by j bits
			{
                strcat(binary,"1");
            } else 
			{
                strcat(binary,"0");
            }
        }
    }
    
    return binary;
}
void setup(mpz_t security_parameter) 
{
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
	//printf("\nParameters are (p + 1 = n * l where p is prime, same as the q in a_param, n is the order of the group): \n");
	//pbc_param_out_str(stdout, par);			// Printing all parameters
	element_t g1, g2, gt, temp;
	element_init_G1(g1, pairing);			//Initialize g1 to be an element of the group G1(input group)
	element_init_G1(g2, pairing);			//Initialize g2 to be an element of the group G1(input group)
	element_init_GT(gt, pairing);			//Initialize gt to be an element of the group GT(the output group)
	
	//element_random(g1);
	element_random(g2);
	
	mpz_t index;
	mpz_init(index);
	mpz_set_ui(index, 1);		
	element_t elements;
	element_init_G1(elements, pairing);
	//int j=20;
	element_t point;

	mpz_t roll_no;
	mpz_init(roll_no);
	mpz_set_ui(roll_no, 49);	// As roll no. was not present in the elliptic curve we take the number near to the roll no.
	mpz_t comp;
	mpz_init(comp);
	//mpz_set_ui(comp, 1);
	printf("\n\nElements of group are (by using generator): \n");	
	do
	{
		element_pow_mpz(elements,g2,index);
		gmp_printf("i = %Zd ", index);
		element_printf("i : %B\n",elements);
		element_init_same_as(point,element_x(elements)); //type
		element_set(point, element_x(elements)); //value
		// element_printf("i : %B\n",point);
		element_to_mpz(comp, point);	
		int ans = mpz_cmp(comp, roll_no);	//Comparing x-coordinate with the roll_no.
		
		if(ans==0)
		{
			//printf("Number near Roll number found!!\n\n");	// If found break the loop.
			element_set(g1, elements);
			//break;
		}
		mpz_add_ui(index, index, 1);

	}while(!element_is0(elements));		// Do the loop until we again encounter the first element
	
	printf("\nPairing of number near roll no. and random element from G1:");
	element_pairing(gt,g1,g2);
	element_printf("\nFirst element : %B\n", g1);
	element_printf("Second element: %B\n", g2);
	element_printf("Applying bilinear pairing on these 2 elements, gt: %B\n", gt);
	
	printf("\nFor verification: \n");
	mpz_t a, b, c;
	mpz_init(a);
	mpz_init(b);
	mpz_init(c);
	
	//gmp_printf("c: %Zd", c);
	mpz_set_ui(a, 2);
	mpz_set_ui(b, 3);
	mpz_mul(c, a, b);
	printf("For a=2 and b=3\n");
	
	element_t g1_new, g2_new, gt_new, g1_, g2_, gt_;
	element_init_G1(g1_new, pairing);			//Initialize g1_new to be an element of the group G1(input group)
	element_init_G1(g2_new, pairing);			//Initialize g2_new to be an element of the group G1(input group)
	element_init_GT(gt_new, pairing);
	element_init_G1(g1_, pairing);			
	element_init_G1(g2_, pairing);
	element_init_GT(gt_, pairing);
	
	element_random(g1_new);
	element_random(g2_new);
	
	element_printf("g1: %B\n", g1_new);
	element_printf("g2: %B\n", g2_new);
	
	element_pow_mpz(g1_, g1_new, a);
	element_pow_mpz(g2_, g2_new, b);
	element_pairing(gt_, g1_, g2_);
	element_printf("\nLHS: e(g1^a, g2^b) %B", gt_);
	
	element_pairing(gt_new, g1_new, g2_new);
	element_pow_mpz(gt_new, gt_new, c);
	element_printf("\nRHS: e(g1, g2)^ab %B\n", gt_new);
	
	
	//For hashing
	
	/*char *str="asdf";		
	char *s1;
	int jk;
	s1=stringToBinary(str);
	
	printf("/nHash code is:");
	size_t l=strlen(s1);
	for(jk=0; jk<l; jk++)
	{
		printf("%c", s1[jk]);
	}
	
	
	element_t h2_value;
	element_init_Zr(h2_value, pairing);
	element_from_hash(h2_value,s1, 44);

	element_printf("Hash value is: %B\n", h2_value);
	
	element_t e;
	element_init_G1(e, pairing);*/

	
	
}

int main () 
{
	mpz_t security_parameter;
	mpz_init(security_parameter);
	mpz_set_ui(security_parameter, 5);		// Setting lambda = 5
	setup(security_parameter);			
	return 0;
}






//mpz_t roll_no,temp1;
	//mpz_init(roll_no);

	//mpz_set_ui(roll_no, 44);
	//mpz_init(temp1);
	/*element_t *point;
	element_t zt;
	mpz_t roll_no,check;
	mpz_init(roll_no);
	mpz_init(check);
	mpz_set_ui(roll_no, 44);
	element_set_mpz(zt, roll_no);*/

	//element_t pt;
	//element_init_G1(pt, pairing);
	//element_set_mpz(pt, roll_no);
	//element_to_mpz(temp1, pt);
	//gmp_printf("%Zd\n",temp1);
//element_cmp(point, [44,0]);
		//point = element_x(elements);
		//element_to_mpz(check, point[1]);
		//int ans = element_cmp(point, zt);
		//int ans = mpz_cmp_ui(check, 44U);

	/*mpz_t z;
	mpz_init(z);
	element_to_mpz(z);
	gmp_printf("mpz %Zd\n", z);*/

	/*element_t ghash, ele;
	mpz_t x;
	int a;
	for(a=0; a<)
	element_pow_mpz(ele, g1, x)*/
	
	//element_t ghash;
	//element_init_G1(ghash, pairing);
	//element_random(ghash);
		//element_init_same_as(point,element_x(elements));
//signed long int i=43;
	//element_set_si(g1, i);
//element_init_G1([43,, pairing);
//mpz_t t;
	//mpz_init(t);
	//mpz_set(t, k);
	//mpz_add_ui(t,t, 1);
	
	//element_t *s;
	//s = element_x(g1);
	//element_printf("Element x of g1: %B\n", s);
	//element_printf("Element of G1 group g1: %B\n", g1);
	//element_printf("Element of G1 group g2: %B\n", g2);	
	//element_printf("Element of GT group gt: %B\n", gt);
	
	//element_pairing(gt,g1,g2);
	//element_printf("2. Element of G1 group g1: %B\n", g1);
	
	
	
	/*if (pairing_is_symmetric(pairing) )
	{
		printf("pairing is symmetric G1 AND G2 are same \n");
	}
	
	element_printf("3. Element of G1 group g1: %B\n", g1);*/

// Libraries
#include <stdio.h>
#include <gmp.h>  // ajouter -lgmp lors de la compilation, 
#include <stdlib.h>
#include <time.h>

// Struct
typedef struct
{
    mpz_t x, y;
}Point_aff;

typedef struct
{
    mpz_t x, y, z;
}Point_jac;

typedef struct
{
    Point_jac c0, c1_pm;
}Ciphertext;


// Constant declaration
#define A 0
#define B 3

// Functiuns declaration


// Function description
    // Operations on mpz_t type numbers
void inv_mult(mpz_t* res, const mpz_t nb, const mpz_t module)
{
    mpz_t exp;
    mpz_init(exp);
    mpz_sub_ui(exp, module, 2);
    mpz_powm(*res, nb, exp, module);

    mpz_clear(exp);   
}

void test_inv_mult()
{
    printf("\nTest inv_mult :\n");
    mpz_t res, nb, module;
    mpz_inits(res, nb, module, NULL);

    mpz_set_ui(nb, 3);
    mpz_set_ui(module, 7);

    inv_mult(&res, nb, module);
    
    gmp_printf("\t(1 / %Zd) %% %Zd = %Zd\n", nb, module, res);

    mpz_clears(res, nb, module, NULL);
}

void square_multiply(mpz_t* res, const mpz_t a, const mpz_t k, const mpz_t module)
{
    mpz_set(*res, a);
    
    for (long i = mpz_sizeinbase(k, 2) - 2; i >= 0; i--)
    {
        mpz_powm_ui(*res, *res, 2, module);
        if (mpz_tstbit(k, i))
        {
            mpz_mul(*res, *res, a);
            mpz_mod(*res, *res, module);
        }
    }
    printf("\n");
}

void test_square_multiply()
{
    printf("\nTest square_multiply :");
    mpz_t res, a, k, module;
    mpz_init(res);
    mpz_init_set_ui(a, 5);
    mpz_init_set_ui(k, 10);
    mpz_init_set_ui(module, 7);

    square_multiply(&res, a, k, module);
    gmp_printf("\t(%Zd ^ %Zd) %% %Zd = %Zd\n", a, k, module, res);

    mpz_clears(res, a, k, module, NULL);
}

    // Operations on points in affine coordinates
void p_aff_printf(const Point_aff p)
{
    gmp_printf("(%Zd, %Zd)\n", p.x, p.y);
}

void p_aff_init_em(Point_aff* p)
{
    mpz_inits(p->x, p->y, NULL);
}

void p_aff_init(Point_aff* p, const unsigned long int x, const unsigned long int y, const unsigned long int z)
{
    mpz_init_set_ui(p->x, x);
    mpz_init_set_ui(p->y, y);
}

void p_aff_clear(Point_aff* p)
{
    mpz_clears(p->x, p->y, NULL);
}

void p_aff_add(Point_aff* p_res, const Point_aff p1, const Point_aff p2, const mpz_t module)
{
    mpz_t inv, sous1, sous2, lamb;
    // On fait des define pour déclarer moins de variable en mémoire possible
    #define mul inv
    #define pow inv
    #define mod sous1
    mpz_inits(inv, sous1, sous2, lamb, NULL);

    mpz_sub(sous1, p2.x, p1.x);
    inv_mult(&inv, sous1, module);
    mpz_mod(mod, inv, module);
    mpz_sub(sous2, p2.y, p1.y);
    mpz_mul(mul, sous2, mod);
    mpz_mod(lamb, mul, module);

    mpz_powm_ui(pow, lamb, 2, module);
    mpz_sub(sous1, pow, p1.x);
    mpz_sub(sous2, sous1, p2.x);
    mpz_mod(p_res->x, sous2, module);

    mpz_sub(sous2, p1.x, p_res->x);
    mpz_mul(mul, lamb, sous2);
    mpz_mod(mod, mul, module);
    mpz_sub(sous2, mod, p1.y);
    mpz_mod(p_res->y, sous2, module);

    #undef mul
    #undef pow
    #undef mod
    mpz_clears(inv, sous1, sous2, lamb, NULL);
}

void test_p_aff_add()
{
    printf("\nTest p_aff_add :\n\t");

    Point_aff p1, p2, p3;
    mpz_t module;
    mpz_inits(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, module, NULL);

    mpz_set_ui(p1.x, 5);
    mpz_set_ui(p1.y, 3);
    mpz_set_ui(p2.x, 1);
    mpz_set_ui(p2.y, 5);
    mpz_set_ui(module, 7);

    p_aff_add(&p3, p1, p2, module);
    
    p_aff_printf(p1);
    printf("\t  +\n\t");
    p_aff_printf(p2);
    printf("\t  =\n\t");
    p_aff_printf(p3);

    mpz_clears(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, module, NULL);
}

void p_aff_double(Point_aff* p_res, const Point_aff p, const mpz_t module)
{
    mpz_t pow, mul, r_add, inv, lamb;
    #define mod1 pow
    #define mod2 r_add
    #define sub inv
    mpz_inits(pow, mul, r_add, inv, lamb, NULL);

    mpz_powm_ui(pow, p.x, 2, module);
    mpz_mul_ui(mul, pow, 3);
    mpz_mod(mod1, mul, module);
    mpz_add_ui(r_add, mod1, A);
    mpz_mul_ui(mul, p.y, 2);
    mpz_mod(mod1, mul, module);
    inv_mult(&inv, mod1, module);
    mpz_mul(mul, r_add, inv);
    mpz_mod(lamb, mul, module);

    mpz_powm_ui(pow, lamb, 2, module);
    mpz_mul_ui(mul, p.x, 2);
    mpz_mod(mod2, mul, module);
    mpz_sub(sub, pow, mod2);
    mpz_mod(p_res->x, sub, module);

    mpz_sub(sub, p.x, p_res->x);
    mpz_mul(mul, lamb, sub);
    mpz_mod(mod1, mul, module);
    mpz_sub(sub, mod1, p.y);
    mpz_mod(p_res->y, sub, module);

    #undef mod1
    #undef mod2
    #undef sub
    mpz_clears(pow, mul, r_add, inv, lamb, NULL);

}

void test_p_aff_double()
{
    printf("\nTest p_aff_double :\n");

    Point_aff p1, p2;
    mpz_t module;
    mpz_inits(p1.x, p1.y, p2.x, p2.y, module, NULL);

    mpz_set_ui(p1.x, 1);
    mpz_set_ui(p1.y, 5);
    mpz_set_ui(module, 7);

    p_aff_double(&p2, p1, module);

    printf("\t 2 * ");
    p_aff_printf(p1);
    printf("\t     =\n\t   ");
    p_aff_printf(p2);

    mpz_clears(p1.x, p1.y, p2.x, p2.y, module, NULL);
}

    // Operations on points in jacobian coordinates
void p_jac_printf(const Point_jac p)
{
    gmp_printf("(%Zd, %Zd, %Zd)\n", p.x, p.y, p.z);
}

void p_jac_cp(Point_jac* dest, const Point_jac src)
{
    mpz_set(dest->x, src.x);
    mpz_set(dest->y, src.y);
    mpz_set(dest->z, src.z);    
}

void p_jac_init_em(Point_jac* p)
{
    mpz_inits(p->x, p->y, p->z, NULL);
}

void p_jac_init(Point_jac* p, const unsigned long int x, const unsigned long int y, const unsigned long int z)
{
    mpz_init_set_ui(p->x, x);
    mpz_init_set_ui(p->y, y);
    mpz_init_set_ui(p->z, z);
}

void p_jac_clear(Point_jac* p)
{
    mpz_clears(p->x, p->y, p->z, NULL);
}

void p_jac_to_p_aff(Point_aff* p_aff, const Point_jac p_jac, const mpz_t module)
{
    mpz_t pow, inv;
    mpz_inits(pow, inv, NULL);

    inv_mult(&inv, p_jac.z, module);
    mpz_powm_ui(pow, inv, 2, module);
    mpz_mul(inv, p_jac.x, pow);
    mpz_mod(p_aff->x, inv, module);

    inv_mult(&inv, p_jac.z, module);
    mpz_powm_ui(pow, inv, 3, module);
    mpz_mul(inv, p_jac.y, pow);
    mpz_mod(p_aff->y, inv, module);

    mpz_clears(pow, inv, NULL);
}

void test_p_jac_to_p_aff()
{
    printf("\nTest p_jac_to_p_aff :\n\t");

    Point_aff p_aff;
    Point_jac p_jac;
    mpz_t module;

    mpz_inits(p_aff.x, p_aff.y, NULL);
    mpz_init_set_ui(p_jac.x, 5);
    mpz_init_set_ui(p_jac.y, 3);
    mpz_init_set_ui(p_jac.z, 1);
    mpz_init_set_ui(module, 7);

    p_jac_to_p_aff(&p_aff, p_jac, module);
    p_aff_printf(p_aff);

    mpz_clears(p_aff.x, p_aff.y, p_jac.x, p_jac.y, p_jac.z, module, NULL);
}

void p_jac_add(Point_jac* p_res, const Point_jac p1, const Point_jac p2, const mpz_t module)
{
    mpz_t pow, mul, mod, u1, u2, s1, s2, h, v;
    #define j u1
    #define i u2
    #define r h
    #define add mul
    #define pow1 mul
    #define sub1 mul
    #define sub mod
    #define mod1 pow
    mpz_inits(pow, mul, mod, u1, u2, s1, s2, h, v, NULL);

    mpz_powm_ui(pow, p2.z, 2, module);
    mpz_mul(mul, p1.x, pow);
    mpz_mod(u1, mul, module);

    mpz_powm_ui(pow, p1.z, 2, module);
    mpz_mul(mul, p2.x, pow);
    mpz_mod(u2, mul, module);

    mpz_mul(mul, p1.y, p2.z);
    mpz_mod(mod, mul, module);
    mpz_powm_ui(pow, p2.z, 2, module);
    mpz_mul(mul, mod, pow);
    mpz_mod(s1, mul, module);

    mpz_powm_ui(pow, p1.z, 3, module);
    mpz_mul(mul, p2.y, pow);
    mpz_mod(s2, mul, module);

    mpz_sub(h, u2, u1);

    mpz_mul_ui(mul, h, 2);
    mpz_mod(mod, mul, module);
    mpz_powm_ui(i, mod, 2, module);

    mpz_mul(mul, u1, i);
    mpz_mod(v, mul, module);

    mpz_mul(mul, h, i);
    mpz_mod(j, mul, module);

    mpz_add(add, p1.z, p2.z);
    mpz_powm_ui(pow, add, 2, module);
    mpz_powm_ui(pow1, p1.z, 2, module);
    mpz_sub(sub, pow, pow1);
    mpz_powm_ui(pow, p2.z, 2, module);
    mpz_sub(mod, sub, pow);
    mpz_mul(mul, mod, h);
    mpz_mod(p_res->z, mul, module);

    mpz_sub(sub, s2, s1);
    mpz_mul_ui(mul, sub, 2);
    mpz_mod(r, mul, module);

    mpz_powm_ui(pow, r, 2, module);
    mpz_sub(sub, pow, j);
    mpz_mul_ui(mul, v, 2);
    mpz_mod(mod1, mul, module);
    mpz_sub(mul, sub, mod1);
    mpz_mod(p_res->x, mul, module);

    mpz_sub(sub, v, p_res->x);
    mpz_mul(mul, r, sub);
    mpz_mod(mod, mul, module);
    mpz_mul_ui(mul, s1, 2);
    mpz_mod(mod1, mul, module);
    mpz_mul(mul, mod1, j);
    mpz_mod(mod1, mul, module);
    mpz_sub(sub1, mod, mod1);
    mpz_mod(p_res->y, sub1, module);

    #undef j
    #undef i
    #undef r
    #undef add
    #undef pow1
    #undef sub1
    #undef sub
    #undef mod1
    mpz_clears(pow, mul, mod, u1, u2, s1, s2, h, v, NULL);
}

void test_p_jac_add()
{
    printf("\nTest p_jac_add :\n\t");

    Point_jac p1, p2, p_res;
    mpz_t module;

    p_jac_init(&p1, 5, 3, 1);
    p_jac_init(&p2, 1, 5, 1);
    p_jac_init_em(&p_res);
    mpz_init_set_ui(module, 7);

    p_jac_add(&p_res, p1, p2, module);

    p_jac_printf(p_res);

    p_jac_clear(&p1);
    p_jac_clear(&p2);
    p_jac_clear(&p_res);
    mpz_clear(module);
}

void p_jac_opp(Point_jac* p_res, const Point_jac p, const mpz_t module)
{
    mpz_t opp;
    mpz_init(opp);

    mpz_set(p_res->x, p.x);
    mpz_sub(opp, module, p.y);
    mpz_set(p_res->y, opp);
    mpz_set(p_res->z, p.z);

    mpz_clear(opp);
}
void p_jac_dbl(Point_jac* p_res, const Point_jac p, const mpz_t module)
{
    mpz_t pow, mul, mod, a, b, c, e, f;
    #define d b
    #define add mod
    #define sub mul
    #define sub1 pow
    #define mod1 pow
    mpz_inits(pow, mul, mod, a, b, c, e, f, NULL);


    mpz_powm_ui(a, p.x, 2, module);
    
    mpz_powm_ui(b, p.y, 2, module);

    mpz_powm_ui(c, b, 2, module);

    mpz_add(add, p.x, b);
    mpz_powm_ui(pow, add, 2, module);
    mpz_sub(sub, pow, a);
    mpz_sub(sub1, sub, c);
    mpz_mul_ui(mul, sub1, 2);
    mpz_mod(d, mul, module);

    mpz_mul_ui(mul, a, 3);
    mpz_mod(e, mul, module);

    mpz_powm_ui(f, e, 2, module);

    mpz_mul_ui(mul, d, 2);
    mpz_mod(mod, mul, module);
    mpz_sub(sub, f, mod);
    mpz_mod(p_res->x, sub, module);

    mpz_sub(sub1, d, p_res->x);
    mpz_mul(mul, e, sub1);
    mpz_mod(mod, mul, module);
    mpz_mul_ui(mul, c, 8);
    mpz_mod(mod1, mul, module);
    mpz_sub(sub, mod, mod1);
    mpz_mod(p_res->y, sub, module);

    mpz_mul_ui(mul, p.y, 2);
    mpz_mod(mod, mul, module);
    mpz_mul(mul, mod, p.z);
    mpz_mod(p_res->z, mul, module);

    #undef d
    #undef add
    #undef sub
    #undef sub1
    #undef mod1
    mpz_clears(pow, mul, mod, a, b, c, e, f, NULL);
}

void test_p_jac_dbl()
{
    printf("\nTest p_jac_dbl :\n\t");

    Point_jac p, p_res;
    mpz_t module;

    p_jac_init_em(&p_res);
    p_jac_init(&p, 5, 3, 1);
    mpz_init_set_ui(module, 7);

    p_jac_dbl(&p_res, p, module);

    p_jac_printf(p_res);

    p_jac_clear(&p_res);
    p_jac_clear(&p);
    mpz_clear(module);
}

void p_jac_mult_scal(Point_jac* p_0, const Point_jac p, const mpz_t k, const mpz_t module)
{
    // We implement it with Montgomery ladder
    Point_jac p_res_calc;
    p_jac_init_em(&p_res_calc);

    p_jac_cp(p_0, p);

    Point_jac p_1;
    p_jac_init_em(&p_1);

    p_jac_dbl(&p_1, *p_0, module);
    int n = mpz_sizeinbase (k, 2);
    
    for (int i = n - 2; i >= 0; i--)
    {
        if (mpz_tstbit(k, i))
        {
            p_jac_add(&p_res_calc, *p_0, p_1, module);
            p_jac_cp(p_0, p_res_calc);
            p_jac_dbl(&p_res_calc, p_1, module);
            p_jac_cp(&p_1, p_res_calc);
        }
        else
        {
            p_jac_add(&p_res_calc, p_1, *p_0, module);
            p_jac_cp(&p_1, p_res_calc);
            p_jac_dbl(&p_res_calc, *p_0, module);
            p_jac_cp(p_0, p_res_calc);           
        }
    }

    p_jac_clear(&p_res_calc);
    p_jac_clear(&p_1);
}

void test_p_jac_mult_scal()
{
    printf("\nTest p_jac_mult_scal :\n\t");

    Point_jac p_res, p;
    Point_aff p_res_aff;
    mpz_t k, module;

    p_jac_init_em(&p_res);
    p_jac_init(&p, 1, 5, 1);
    p_aff_init_em(&p_res_aff);
    mpz_init_set_ui(k, 3);
    mpz_init_set_ui(module, 7);

    p_jac_mult_scal(&p_res, p, k, module);
    p_jac_to_p_aff(&p_res_aff, p_res, module);

    p_jac_printf(p);
    gmp_printf("\t  * %Zd = \n\t", k);
    p_jac_printf(p_res);
    printf("\t    =\n\t  ");
    p_aff_printf(p_res_aff);

    p_jac_clear(&p_res);
    p_jac_clear(&p);
    p_aff_clear(&p_res_aff);
    mpz_clears(k, module, NULL);
}

void alice_bob()
{
    time_t t;
    const int ordre = 13;
    mpz_t a, b, module;
    Point_aff p_A1_aff, p_B1_aff;
    Point_jac p_j, p_A, p_B, p_A1, p_B1;

    srand((unsigned) time(&t));

    mpz_init_set_ui(a, (unsigned int) ((rand() % ordre) + 1));
    mpz_init_set_ui(b, (unsigned int) ((rand() % ordre) + 1));
    mpz_init_set_ui(module, 7);

    p_aff_init_em(&p_A1_aff);
    p_aff_init_em(&p_B1_aff);

    p_jac_init(&p_j, 5, 3, 1);
    p_jac_init_em(&p_A);
    p_jac_init_em(&p_B);
    p_jac_init_em(&p_A1);
    p_jac_init_em(&p_B1);

    p_jac_mult_scal(&p_A, p_j, a, module);
    p_jac_mult_scal(&p_B, p_j, b, module);

    p_jac_mult_scal(&p_A1, p_A, b, module);
    p_jac_mult_scal(&p_B1, p_B, a, module);

    p_jac_to_p_aff(&p_A1_aff, p_A1, module);
    p_jac_to_p_aff(&p_B1_aff, p_B1, module);


    gmp_printf("\nAlice et Bob \n\ta : %Zd, b : %Zd\n\t", a, b);
    p_aff_printf(p_A1_aff);
    printf("\t");
    p_aff_printf(p_B1_aff);

    mpz_clears(a, b, module, NULL);
    p_aff_clear(&p_A1_aff);
    p_aff_clear(&p_B1_aff);
    p_jac_clear(&p_j);
    p_jac_clear(&p_A);
    p_jac_clear(&p_B);
    p_jac_clear(&p_A1);
    p_jac_clear(&p_B1);
}

void cpt_aff(const Ciphertext cpt)
{
    p_jac_printf(cpt.c0);
    p_jac_printf(cpt.c1_pm);
}

// EC Elgamal
void f_mpz_to_p_aff(Point_aff* p, const mpz_t nb, const mpz_t module)
{
    mpz_t sub, r, add, k;

    mpz_inits(sub, r, add, k, NULL);

    mpz_sub_ui(sub, module, 3);
    mpz_mod_ui (r, sub, 4);

    if (mpz_cmp_ui(r, 0)) // Check if our module has the form " module = 4k+3" 
    {
        mpz_cdiv_q_ui(k, sub, 4);
        
        mpz_set(p->x, nb);

        mpz_add_ui(add, k, 1);

        //Trouver y^2  = x^3 + 3, puis on peut trouver y avec 

        mpz_powm(p->y, add, module);
    }
    else
        printf("Module is not supported. Need : module = 4k+3.\n");

    mpz_clears(sub, r, add, k, NULL);
}

void ec_elgamal()
{
    // Initialisation
    Point_jac p, y;
    Point_aff y_pub;
    mpz_t ordre, ordre_m1, module, x, rdm;
    gmp_randstate_t state;

    p_jac_init(&p, 6, 4, 1);
    p_jac_init_em(&y);
    p_aff_init_em(&y_pub);
    mpz_init_set_ui(ordre, 13);
    mpz_init_set_ui(module, 7);
    mpz_inits(x, ordre_m1, rdm, NULL);

    mpz_sub_ui(ordre_m1, ordre, 1);

    // Generate public key
    // 1.1.Generate random x
    unsigned long seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    mpz_urandomm(rdm, state, ordre_m1);
    mpz_add_ui(x, rdm, 1);

    // 1.2. Calculate our public key
    p_jac_mult_scal(&y, p, x, module);
    p_jac_to_p_aff(&y_pub, y, module);

    printf("Public key is : ");
    p_aff_printf(y_pub);
    gmp_printf("x = %Zd\n", x);

    // Encryption
    Ciphertext cpt;
    Point_jac c1, pm;
    Point_aff pm_aff;
    mpz_t k;

    mpz_init(k);

    p_jac_init_em(&cpt.c0);
    p_jac_init_em(&cpt.c1_pm);
    p_jac_init_em(&c1);
    p_aff_init_em(&pm_aff);

    // 2.1. Generate random k
    mpz_urandomm(rdm, state, ordre_m1);
    mpz_add_ui(k, rdm, 1);


    p_jac_mult_scal(&cpt.c0, p, k, module);
    p_jac_mult_scal(&c1, y, k, module);

    p_jac_init(&pm, 2, 5, 1);

    p_jac_add(&cpt.c1_pm, c1, pm, module);

    printf("Encryption :\n\tPlain :");
    p_jac_printf(pm);
    p_jac_to_p_aff(&pm_aff, pm, module);
    p_aff_printf(pm_aff);
    printf("\tCiphertext :");
    cpt_aff(cpt);

    // Decryption
    Point_jac minus_c1;
    p_jac_init_em(&minus_c1);

    p_jac_mult_scal(&c1, cpt.c0, x, module);
    p_jac_opp(&minus_c1, c1, module);
    p_jac_add(&pm, cpt.c1_pm, minus_c1, module);

    printf("Decryption :\n\tPlain :");
    p_jac_printf(pm);
    p_jac_to_p_aff(&pm_aff, pm, module);
    p_aff_printf(pm_aff);

    //Faire les clear
    p_jac_clear(&p);
    p_jac_clear(&y);
    p_aff_clear(&y_pub);
    p_jac_clear(&cpt.c0);
    p_jac_clear(&cpt.c1_pm);
    p_jac_clear(&c1);
    p_aff_clear(&pm_aff);
    p_jac_clear(&pm);
    p_jac_clear(&minus_c1);

    mpz_clears(x, k, module, ordre, ordre_m1, rdm, NULL);
}

int main(void)
{
/*     test_inv_mult();
    test_square_multiply();

    test_p_aff_add();
    test_p_aff_double();

    test_p_jac_to_p_aff();
    test_p_jac_add();
    test_p_jac_dbl();
    test_p_jac_mult_scal();

    alice_bob(); */

    //ec_elgamal();

    Point_aff p;
    mpz_t nb, module;

    p_aff_init_em(&p);
    mpz_init_set_ui(nb, 2);
    mpz_init_set_ui(module, 7);

    f_mpz_to_p_aff(&p, nb, module);

    p_aff_printf(p);

    return 0;
}
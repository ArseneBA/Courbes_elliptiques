// Libraries
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "affine.h"
#include "jacobien.h"


// Tests functions description
    // Operations on mpz_t type numbers
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

// Test Diffie Hellman with Alice and bob
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


// EC Elgamal
typedef struct
{
    Point_jac c0, c1_pm;
}Ciphertext;

void cpt_aff(const Ciphertext cpt)
{
    p_jac_printf(cpt.c0);
    p_jac_printf(cpt.c1_pm);
}

void f_mpz_to_p_aff(Point_aff* p, const mpz_t nb, const mpz_t module)
{
    mpz_t sub, r, add, k, pow;

    mpz_inits(sub, r, add, k, pow, NULL);

    mpz_sub_ui(sub, module, 3);
    mpz_mod_ui (r, sub, 4);

    if (mpz_cmp_ui(r, 0) == 0) // Check if our module has the form " module = 4k+3" 
    {
        mpz_cdiv_q_ui(k, sub, 4);

        mpz_set(p->x, nb);

        // Find y^2  = x^3 + 3 
        mpz_powm_ui(pow, nb, 3, module);
        mpz_add_ui(sub, pow, 3);
        mpz_mod(pow, sub, module);  // y^2

        // y = (y^2)^(k+1) [module]
        mpz_add_ui(add, k, 1);  // k+1 
        mpz_powm(p->y, sub, add, module);
    }
    else
        printf("Module is not supported. Need : module = 4k+3.\n");

    mpz_clears(sub, r, add, k, pow, NULL);
}

void ec_elgamal()
{
    printf("\n\nElGamal with Elliptic curve :\n\n");
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
    //      Generate random x
    unsigned long seed = time(NULL);
    gmp_randinit_default(state);
    gmp_randseed_ui(state, seed);

    mpz_urandomm(rdm, state, ordre_m1);
    mpz_add_ui(x, rdm, 1);

    //      Calculate our public key
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

    //      Generate random k
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

    // Clears
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
    test_inv_mult();
    test_square_multiply();

    test_p_aff_add();
    test_p_aff_double();

    test_p_jac_to_p_aff();
    test_p_jac_add();
    test_p_jac_dbl();
    test_p_jac_mult_scal();

    alice_bob();

    ec_elgamal();

    return 0;
}

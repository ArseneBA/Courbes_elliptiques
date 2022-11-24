// Libraries
#include <stdio.h>
#include <gmp.h>  // ajouter -lgmp lors de la compilation, 

// Struct
typedef struct
{
    mpz_t x, y;
}Point_aff;

typedef struct
{
    mpz_t x, y, z;
}Point_jac;

// Constant declaration
#define A 0
#define B 3

// Functiuns declaration


// Function description
void aff_point_affine(Point_aff p)
{
    gmp_printf("(%Zd, %Zd)\n", p.x, p.y);
}

void inv_mult(mpz_t* res, mpz_t nb, const mpz_t module)
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
}

void add_aff(Point_aff* p_res, const Point_aff p1, const Point_aff p2, const mpz_t module)
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

void test_add_aff()
{
    printf("\nTest add_aff :\n\t");

    Point_aff p1, p2, p3;
    mpz_t module;
    mpz_inits(p1.x, p1.y, p2.x, p2.y, p3.x, p3.y, module, NULL);

    mpz_set_ui(p1.x, 5);
    mpz_set_ui(p1.y, 3);
    mpz_set_ui(p2.x, 1);
    mpz_set_ui(p2.y, 5);
    mpz_set_ui(module, 7);

    add_aff(&p3, p1, p2, module);
    
    aff_point_affine(p1);
    printf("\t  +\n\t");
    aff_point_affine(p2);
    printf("\t  =\n\t");
    aff_point_affine(p3);
}

void double_aff(Point_aff* p_res, const Point_aff p, const mpz_t module)
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

void test_double_aff()
{
    printf("\nTest double_aff :\n");

    Point_aff p1, p2;
    mpz_t module;
    mpz_inits(p1.x, p1.y, p2.x, p2.y, module, NULL);

    mpz_set_ui(p1.x, 1);
    mpz_set_ui(p1.y, 5);
    mpz_set_ui(module, 7);

    double_aff(&p2, p1, module);

    printf("\t 2 * ");
    aff_point_affine(p1);
    printf("\t     =\n\t   ");
    aff_point_affine(p2);
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
}

int main(void)
{
    test_inv_mult();
    test_add_aff();
    test_double_aff();
    test_square_multiply();

    return 0;
}
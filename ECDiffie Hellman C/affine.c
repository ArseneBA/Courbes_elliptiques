#include "affine.h"

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

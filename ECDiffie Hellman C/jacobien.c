#include "jacobien.h"

// Function description
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

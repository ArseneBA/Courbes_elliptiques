#ifndef JACOBIEN_H
#define JACOBIEN_H

// Libraries
#include <stdio.h>
#include <gmp.h>

#include "affine.h"


// Stuct
typedef struct
{
    mpz_t x, y, z;
}Point_jac;


// Function declaration
    // Operations on points in jacobian coordinates
void p_jac_printf(const Point_jac p);
void p_jac_cp(Point_jac* dest, const Point_jac src);
void p_jac_init_em(Point_jac* p);
void p_jac_init(Point_jac* p, const unsigned long int x, const unsigned long int y, const unsigned long int z);
void p_jac_clear(Point_jac* p);
void p_jac_to_p_aff(Point_aff* p_aff, const Point_jac p_jac, const mpz_t module);
void p_jac_add(Point_jac* p_res, const Point_jac p1, const Point_jac p2, const mpz_t module);
void p_jac_opp(Point_jac* p_res, const Point_jac p, const mpz_t module);
void p_jac_dbl(Point_jac* p_res, const Point_jac p, const mpz_t module);
void p_jac_mult_scal(Point_jac* p_0, const Point_jac p, const mpz_t k, const mpz_t module);


// Constant declaration
#define A 0
#define B 3

#endif
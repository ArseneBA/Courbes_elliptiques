#ifndef AFFINE_H
#define AFFINE_H

// Libraries
#include <stdio.h>
#include <gmp.h>


// Struct
typedef struct
{
    mpz_t x, y;
}Point_aff;


// Function declaration
    // Operations on mpz_t type numbers
void inv_mult(mpz_t* res, const mpz_t nb, const mpz_t module);
void square_multiply(mpz_t* res, const mpz_t a, const mpz_t k, const mpz_t module);

    // Operations on points in affine coordinates
void p_aff_printf(const Point_aff p);
void p_aff_init_em(Point_aff* p);
void p_aff_init(Point_aff* p, const unsigned long int x, const unsigned long int y, const unsigned long int z);
void p_aff_clear(Point_aff* p);
void p_aff_add(Point_aff* p_res, const Point_aff p1, const Point_aff p2, const mpz_t module);
void p_aff_double(Point_aff* p_res, const Point_aff p, const mpz_t module);


// Constant declaration
#define A 0
#define B 3

#endif
/*
 * */

#ifndef __FACTOR_H_
#define __FACTOR_H_

#include "gmp_helper.h"

// number up to which to attempt trial division:
const unsigned long TRIAL_LIMIT = 1000000;

void factor(mpz_t&, mpz_t&);

bool trial_division(mpz_t&, mpz_t&);
bool fermat_factorization(mpz_t&, mpz_t&);
bool pollard_rho_factorization(mpz_t&, mpz_t&);
bool pollard_p_1_factorization(mpz_t&, mpz_t&);

void gcd(mpz_t&, mpz_t&, mpz_t&);

#endif

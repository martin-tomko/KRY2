/*
 * */

#ifndef __FACTOR_H_
#define __FACTOR_H_

#include "gmp_helper.h"
#include "elliptic_curves.h"
#include "sieve.h"

// number up to which to attempt trial division:
const unsigned long TRIAL_LIMIT = 1000000;

void factor(mpz_t&, mpz_t&);

bool trial_division(mpz_t&, mpz_t&);
bool fermat_factorization(mpz_t&, mpz_t&);
bool pollard_rho_factorization(mpz_t&, mpz_t&);
bool pollard_p_1_factorization(mpz_t&, mpz_t&);


bool lenstra_ec_factorization(mpz_t&, mpz_t&);
bool pick_curve(ECPoint&, mpz_t&, mpz_t&);
  // Cummulative least common multiple - of integers up to K
void lcm_cum(mpz_t&, unsigned long, mpz_t&);
bool is_m_to_r(mpz_t&, mpz_t&, mpz_t&);

#endif

/*
 * */

#ifndef __FACTOR_H_
#define __FACTOR_H_

#include "gmp_helper.h"
#include "elliptic_curves.h"
#include "sieve.h"

// number up to which to attempt trial division:
const unsigned long TRIAL_LIMIT = 1000000;

// Lenstra EC parameters:
const unsigned long B0 = 100000;
const double Q_TRIAL = 1.02;
const unsigned PHASE2_MULTIPLIER = 10;

const unsigned E_SELECTION[] = {1, 2, 3, 6, 12, 18, 24, 30, 60, 0};

void factor(mpz_t&, mpz_t&);

bool trial_division(mpz_t&, mpz_t&);
bool fermat_factorization(mpz_t&, mpz_t&);
bool pollard_rho_factorization(mpz_t&, mpz_t&);
bool pollard_p_1_factorization(mpz_t&, mpz_t&);
bool brent_factorization(mpz_t&, mpz_t&);

inline bool is2pow(unsigned long i) { return (i & (i-1)) == 0; }

inline void apply_f(mpz_t& Y, mpz_t& C, mpz_t& N, mpz_t& tmp)
  /* Y := Y^2 + C mod N*/
  { mpz_set(tmp, C); mpz_addmul(tmp, Y, Y); mpz_mod(Y, tmp, N); }
  //{ mpz_mul(Y, Y, Y); mpz_add(Y, Y, C); mpz_mod(Y, Y, N); }



bool lenstra_ec_factorization(mpz_t&, mpz_t&);
bool pick_curve(ECPoint&, mpz_t&, mpz_t&);
  // Cummulative least common multiple - of integers up to K
void lcm_cum(mpz_t&, unsigned long, mpz_t&);
bool is_m_to_r(mpz_t&, mpz_t&, unsigned long);
bool ecm_try_mult(mpz_t&, ECPoint&, unsigned long);

bool ecm_phase_one(mpz_t&, ECPoint&, unsigned long);
bool ecm_phase_two(mpz_t&, ECPoint&, unsigned long);
unsigned choose_e(unsigned long);


#endif

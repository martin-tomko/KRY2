/*
 * */

#ifndef __PRIMES_H_
#define __PRIMES_H_

#include "gmp_helper.h"

const unsigned SECURITY_PARAMETER = 100;

  // probabilistic test, employs Miller-Rabin
bool is_probably_prime(mpz_t&);

#endif

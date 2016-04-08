/*
 * */

#include "factor.h"
#include "kry.h"
#include <iostream>

void factor(mpz_t& P, mpz_t& N) {
  /* Find a prime factor of N: */
  if (mpz_cmp_ui(N, 2) < 0) { // if N < 2
    mpz_set(P, N);
    return;
  }

  // Attempt to find a factor using trial division:
  if (trial_division(P, N)) { 
    std::cerr << "Found using Trial division" << std::endl;
    return;
  }

  // Attempt to find a factor using Pollard p-1 factorization:
  if (pollard_p_1_factorization(P, N)) {
    std::cerr << "Found using Pollard p-1" << std::endl;
    return;
  }

  // Attempt to find a factor using Fermat factorization:
  if (fermat_factorization(P, N)) {
    std::cerr << "Found using Fermat" << std::endl;
    return;
  }

  // Attempt to find a factor using Pollard rho factorization:
  if (pollard_rho_factorization(P, N)) {
    std::cerr << "Found using Pollard rho" << std::endl;
    return;
  }

  // No factor was found, N is a prime (or the above does not implement a full
  // search), return N:
  mpz_set(P, N);
}

bool trial_division(mpz_t& P, mpz_t& N) {
  // Set the limit of trial division:
  unsigned long limit = TRIAL_LIMIT;
  mpz_sqrt(P, N);
  if (mpz_cmp_ui(P, limit) < 0) {
    limit = mpz_get_ui(P);
  }
  // limit = min(sqrt(N), TRIAL_LIMIT)
  

  // Trial division:

  /* As every prime (except 2 and 3) is of the form (6k-1) or (6k+1), where
   * k is a natural number, the algorithm is simplified as follows:
   * first, division by 2 and 3 is tested, and then for every multiple k of 6
   * (up to a limit) division by k-1 and k+1 is tested.
   * This leads to a 3x decrease in division tests needed compared to a naive
   * trial division. */
  if (mpz_divisible_ui_p(N, 2)) {
    mpz_set_ui(P, 2);
    return true;
  }
  if (mpz_divisible_ui_p(N, 3)) {
    mpz_set_ui(P, 3);
    return true;
  }
  for (unsigned long k = 6; k < TRIAL_LIMIT; k += 6) {
    if (mpz_divisible_ui_p(N, k-1)) {
      mpz_set_ui(P, k-1);
      return true;
    }
    if (mpz_divisible_ui_p(N, k+1)) {
      mpz_set_ui(P, k+1);
      return true;
    }
  }

  return false; // Trial division found no factors.
}

bool fermat_factorization(mpz_t& P, mpz_t& N) {
  bool found = false;
  GMPNum gmp[5];
  
  mpz_t &X  = gmp[0],
        &XX = gmp[1],
        &Y  = gmp[2],
        &YY = gmp[3],
        &R  = gmp[4];

  mpz_sqrtrem(X, R, N);
  if (mpz_cmp_ui(R, 0) != 0) {
    mpz_add_ui(X, X, 1);
  }
  // X = ceil(sqrt(N))

  for(; mpz_cmp(X, N) <= 0; mpz_add_ui(X, X, 1)) {
//    if (mpz_divisible_ui_p(X, 1000)) {
//      print_hex_mpz(X); std::cout << std::endl;
//    }
    mpz_mul(XX, X, X);
    mpz_sub(YY, XX, N);          // YY = XX - N
    mpz_sqrtrem(Y, R, YY);       // compute sqrt and reminder of YY
    if (mpz_cmp_ui(R, 0) == 0) { // YY is square, a factor was found:
      mpz_sub(P, X, Y);          // X - Y
      found = true;
      break;
    }
    mpz_add_ui(X, X, 1);         // X = X + 1
  }

  return found;
}

bool pollard_rho_factorization(mpz_t& P, mpz_t& N) {
  /* Pollard's rho algorithm using the polynomial g(x) = x^2 + 1 */

  GMPNum gmp[2];
  mpz_t &X = gmp[0],
        &Y = gmp[1];

  mpz_set_ui(X, 2);
  mpz_set_ui(Y, 2);
  mpz_set_ui(P, 1);
  while (mpz_cmp_ui(P, 1) == 0) {
    // X = g(X):
    mpz_mul(X, X, X);
    mpz_add_ui(X, X, 1);

    // Y = g(g(Y)):
    mpz_mul(Y, Y, Y);
    mpz_add_ui(Y, Y, 1);
    mpz_mul(Y, Y, Y);
    mpz_add_ui(Y, Y, 1);

    // P = gcd(|X-Y|, n):
    mpz_sub(P, X, Y);
    mpz_abs(P, P);
    gcd(P, P, N);
  }

  return (mpz_cmp(P, N) != 0); // return P == N ? false : true;
}

bool pollard_p_1_factorization(mpz_t& P, mpz_t& N) {
  Randomizer rnd;

  GMPNum gmp[2];
  mpz_t &fct = gmp[0],
        &k   = gmp[1];

  // TODO: upper bound
  // TODO: optimize
  while(true) {
    //mpz_set_ui(k, TRIAL_LIMIT);
    mpz_set_ui(k, 2);

    // choose A:
    mpz_sub_ui(P, N, 2);     // P = N-2
    mpz_urandomm(fct, rnd, P);    // 0 <= fct < N-2
    mpz_add_ui(fct, fct, 2); // 2 <= fct < N
    mpz_set(P, fct);

    while(mpz_cmp(P, N) != 0) {
      mpz_powm(fct, fct, k, N);
      mpz_sub_ui(P, fct, 1);
      gcd(P, P, N);
      if ((mpz_cmp_ui(P, 1) != 0) && (mpz_cmp(P, N) != 0)) {
        return true;
      }
      mpz_add_ui(k, k, 1);
    }
  }

  return false;
}

void gcd(mpz_t& R, mpz_t& M, mpz_t& N) {
  /* GCD using Euclid's algorithm */

  static GMPNum gmp[2];
  static mpz_t &A = gmp[0],
               &B = gmp[1];

  // Create copies of M and N, so as not to overwrite them:
  mpz_set(A, M);
  mpz_set(B, N);

  // Compute:
  while (mpz_cmp_ui(B, 0) != 0) {
    mpz_set(R, B);
    mpz_mod(B, A, B);
    mpz_set(A, R);
  }

  // If the loop finished, R contains gcd(A, B); no further assignment needed.
}

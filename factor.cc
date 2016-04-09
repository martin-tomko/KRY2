/*
 * */

#include "factor.h"
#include "kry.h"
#include <iostream>
#include <cmath>
//#include "sieve.h"

#include "elliptic_curves.h"

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
  if (lenstra_ec_factorization(P, N)) {
    std::cerr << "Found using Lenstra EC factorization" << std::endl;
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

  /*
  Randomizer rnd;
  size_t bound1 = 10000000; // = 10000000;
//       bound2 = bound1 * bound1;
  const bool_vec& is_prime = Sieve::get_vector(bound1);

  size_t pexp;
//  double logb;
  GMPNum gmp[2];
  mpz_t &m = gmp[0],
        &a = gmp[1];

  // choose a:
  mpz_sub_ui(P, N, 2);     // P = N-2
  mpz_urandomm(a, rnd, P); // 0 <= fct < N-2
  mpz_add_ui(a, a, 2);     // 2 <= fct < N

  while (true) {
    Sieve::stretch(bound1);

    // Compute m = lcm(1, ..., bound1):
    mpz_set_ui(m, 1);
    for (size_t p = 2; p < bound1; p++) { // for each prime:
      if (!is_prime[p]) { // skip composite numbers
        continue;
      }

      pexp = p;
      while (p < bound1) {
        pexp *= p;
      }
      pexp /= p;

      mpz_mul_ui(m, m, pexp);
    }

    


    gcd(P, a, N);
    if (mpz_cmp_ui(P, 1) != 0) { // if gcd(a,n) != 1, we found a factor.
      return true;
    }

    mpz_powm(a, a, m, N); // a = a^m mod N
    mpz_sub_ui(a, a, 1);  // a = a - 1
    gcd(P, a, N);

    if (mpz_cmp_ui(P, 1) == 0) { // gcd == 1, increase bound and continue
      bound1 *= 2;
    } else if (mpz_cmp(P, N) == 0) { // gcd == N, change a
      // choose a:
      mpz_sub_ui(P, N, 2);     // P = N-2
      mpz_urandomm(a, rnd, P); // 0 <= fct < N-2
      mpz_add_ui(a, a, 2);     // 2 <= fct < N
    } else { // otherwise, a factor has been found.
      return true;
    }

  }
  */

//  /*
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
/*
    mpz_sub_ui(P, N, 2);     // P = N-2
    mpz_urandomm(fct, rnd, P);    // 0 <= fct < N-2
    mpz_add_ui(fct, fct, 2); // 2 <= fct < N
*/
    rnd.get_random(fct, 2, N);
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
//  */
}

bool lenstra_ec_factorization(mpz_t& P, mpz_t& N) {
  /* It is assumed that N was already analyzed with a simpler method
   * and is not divisible by 2 or 3. */

  // 1. TODO: check whether m^r = n for some r?
  unsigned long K = 1000000;
  unsigned attempts = 0;

  GMPNum gmp[1];
  mpz_t &k = gmp[0];
  lcm_cum(k, K, N);
  mpz_mod(k, k, N);

  std::cout << "Lenstra is happening." << std::endl;
  while (true) {
    ECPoint E;
    // Pick curve:
    if (!pick_curve(E, P, N)) {
      std::cout << "We are done, somehow" << std::endl;

      // A factor was found and written to P.
      return true;
    }

  
    // Print point:
    std::cout << "Observe the point:" <<  std::endl;
    print_hex_mpz(E.get_x());
    std::cout << std::endl;
    print_hex_mpz(E.get_y());
    std::cout << std::endl;
    print_hex_mpz(E.get_a());
    std::cout << std::endl;


    if (! E.pow(P, k)) { // Addition failed.
      if (mpz_cmp(P, N) == 0) {
        K /= 2;
        break;
      } else { // A factor was found:
        return true;
      }
    }

    if (E.is_inf()) {
      continue;
    }

    attempts++;
    if ((attempts & 255) == 255) {
      std::cerr << "Increasing K from " << K << std::endl;
      K *= 2;
      lcm_cum(k, K, N);
      mpz_mod(k, k, N);
    }
  }

  return false;
}

bool pick_curve(ECPoint& E, mpz_t& remainder, mpz_t& N) {
  static Randomizer rnd;

  GMPNum gmp[3];
  mpz_t &x = gmp[0],
        &y = gmp[1],
        &a = gmp[2];

  while(true) {
    rnd.get_random(x, 2, N);
    rnd.get_random(y, 2, N);
    rnd.get_random(a, 2, N);

    E.generate(a, x, y, N);
    if (E.nonsingular(remainder)) { // E is nonsingular, return.
      return true;
    }

    if (mpz_cmp(remainder, N) != 0) { // A factor was found.
      return false;
    }

    // Otherwise, we must try again.
  }
}

void lcm_cum(mpz_t& lcm, unsigned long K, mpz_t& N) {
  Sieve& is_prime = Sieve::get_instance(K);

  unsigned m;

  mpz_set_ui(lcm, 1);
  for (unsigned p = 2; p < K; p++) {
    if (is_prime[p]) {
      m = p;
      while (m < K) {
        m *= p;
      }
  //    m /= p;

      mpz_mul_ui(lcm, lcm, m);
      if ((p & 255) == 255) {
        mpz_mod(lcm, lcm, N);
      }
    }
  }
}

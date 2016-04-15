/*
 * */

#include "primes.h"

bool is_probably_prime(mpz_t& N) {
  Randomizer rnd;
  GMPNum gmp[4];
  mpz_t &N1 = gmp[0],
        &R  = gmp[1],
        &A  = gmp[2],
        &Y  = gmp[3];

  unsigned t = SECURITY_PARAMETER,
           s;

  if (mpz_cmp_ui(N, 1) <= 0) {             // If N <= 1:
    // Return false, no primes below 1 exist.
    return false;
  } else if (mpz_divisible_2exp_p(N, 1)) { // If 2|N:

    // Return true if N == 2, false otherwise.
    return (mpz_cmp_ui(N, 2) == 0)
      ? true
      : false;
    /* verbose for clarity */
  }

  // From here on we are sure N is an odd integer and N >= 3.
  // Perform the Miller-Rabin test:

  // Write N-1 as 2^s * R, R odd:
  mpz_sub_ui(N1, N, 1);      // N = N-1
  mpz_tdiv_q_2exp(R, N1, 1); // R = (N-1) / 2
  s = 1;
  while (mpz_divisible_2exp_p(R, 1)) { // While 2|R:
    mpz_tdiv_q_2exp(R, R, 1);
    s++;
  }
  // N-1 == 2^s * R, R is odd

  for (unsigned i = 0; i < t; i++) { // Perform t random tests:
    rnd.get_random(A, 2, N1); // 2 <= A <= N-2

    mpz_powm(Y, A, R, N);     // Y = A^R mod N
    /* Our love, our love, in a ball of yarn */

    if ((mpz_cmp_ui(Y, 1) != 0) && (mpz_cmp(Y, N1) != 0)) {
      // If Y != 1 and Y != N-1:

      for (unsigned j = 1; (j < s) && (mpz_cmp(Y, N1) != 0); j++) {
        // For j in [1, s), while Y != N-1:
        mpz_powm_ui(Y, Y, 2, N);
        if (mpz_cmp_ui(Y, 1) == 0) { // If Y == 1:
          // N is definitely composite:
          return false; 
        }
      }

      if (mpz_cmp(Y, N1) != 0) { // If Y != N-1:
        // N is definitely composite:
        return false;
      }
    }
  }

  // Passed all attempts to disprove primality, N is probably prime:
  return true;
  /* return very_likely; */
}

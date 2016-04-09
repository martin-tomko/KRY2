/*
 * */

#include "gmp_helper.h"
#include <fstream>
#include <iostream>


void Randomizer::seed() {
  // TODO: get seed somehow
  //   1. /dev/random, if not avaliable, then:
  //   2. system time, if not available, then:
  //   3. whatever

  unsigned long seed = 48321841;
  // If all else fails, the generator will be initialized with this number.
  
  std::ifstream rnd_stream;
  try {
    rnd_stream.open("/dev/random", std::ios::in | std::ios::binary);
    rnd_stream.read(reinterpret_cast<char*>(&seed), sizeof(seed));
    rnd_stream.close();
  } catch (...) {
    if (rnd_stream.is_open()) {
      rnd_stream.close();
    }
  }
  
  gmp_randseed_ui(state, seed); 
}

void Randomizer::get_random(mpz_t& n, unsigned long min_val, mpz_t& max_val) {
  mpz_sub_ui(n, max_val, min_val);
  mpz_urandomm(n, state, n);
  mpz_add_ui(n, n, min_val);
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

bool invert(mpz_t& result, mpz_t& a, mpz_t& n) {
  static GMPNum gmp[5];
  static mpz_t &t   = gmp[0],
               &tn  = gmp[1],
               &r   = gmp[2],
               &rn  = gmp[3],
               &quo = gmp[4];



  // Initialize:
  mpz_set_ui(t, 0);
  mpz_set_ui(tn, 1);
  mpz_set(r, n);
  mpz_set(rn, a);

  while (mpz_cmp_ui(rn, 0) != 0) {
    // quo = r / rn:
    mpz_tdiv_q(quo, r, rn);

    // (r, rn) = (rn, r - quo*rn):
    mpz_swap(r, rn);
    mpz_submul(rn, quo, r);

    // (t, tn) = (tn, t - quo*tn):
    mpz_swap(t, tn);
    mpz_submul(tn, quo, t);
  }

  if (mpz_cmp_ui(r, 1) > 0) { // r > 1, a is not invertible
    return false;
  }

  if (mpz_cmp_ui(t, 0) < 0) { // t < 0, needs to be modulated:
    mpz_add(result, t, n);
  } else {
    mpz_set(result, t);
  }
  return true;
}

/*
 * */

#ifndef __GMP_HELPER_H_
#define __GMP_HELPER_H_

#include <cstddef>
#include <gmp.h>
#include <gmpxx.h>


// wrapper for static initialization / deinitialization
class GMPNum {
 public:
  GMPNum() { mpz_init(value); }
  ~GMPNum() { mpz_clear(value); }
  operator mpz_t&() { return value; } // automatic conversion to mpz_t&
 private:
  mpz_t value;
};

class Randomizer {
 public:
  Randomizer() { gmp_randinit_mt(state); seed(); } // Mersenne Twister
  ~Randomizer() { gmp_randclear(state); }
  void seed();
  operator gmp_randstate_t&() { return state; }
 private:
  gmp_randstate_t state;
  
};

#endif

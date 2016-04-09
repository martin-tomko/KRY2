/*
 * */

#ifndef __ELLIPTIC_CURVES_H
#define __ELLIPTIC_CURVES_H

#include "gmp_helper.h"

// TODO: point "at infinity"! not "infinite" point
//       (but the rest are "finite" points

class ECPoint {
 public:
  ECPoint(bool is_inf) : valid {is_inf}, inf {is_inf} {}
  ECPoint() : valid {false}, inf {false} {}
  ECPoint(ECPoint& other)             // copy ctr
    { duplicate(other); }
  ECPoint& operator =(ECPoint& other) // assignment operator
    { duplicate(other); return *this; }

  void duplicate(ECPoint&);

  bool operator ==(ECPoint&); 

  // initialize point with A, X, Y, N
  void generate(mpz_t&, mpz_t&, mpz_t&, mpz_t&);

  bool add(mpz_t&, ECPoint&);
  bool pow(mpz_t&, mpz_t&);
  bool pow(mpz_t&, unsigned long);

  bool is_valid() { return valid; }
  bool is_inf() { return inf; }
  bool nonsingular(mpz_t&);

  void set_inf() { inf = true; valid = true; }

//  mpz_t& set_x(mpz_t& new_x) { mpz_set(x, new_x); }
//  mpz_t& set_y(mpz_t& new_y) { mpz_set(y, new_y); }

  mpz_t& get_x() { return x; }
  mpz_t& get_y() { return y; }
  mpz_t& get_a() { return a; }
  mpz_t& get_b() { return b; }
  mpz_t& get_n() { return n; }
 private:
  bool valid;
  bool inf;
  GMPNum a, b, n, x, y;
};


void test_curves();

#endif

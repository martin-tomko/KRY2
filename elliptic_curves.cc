/*
 * */

#include "elliptic_curves.h"
#include <algorithm>

#include <iostream>
#include "kry.h"
using namespace std;
/*
ECPoint::ECPoint(const ECPoint& other) {
  if (other.is_valid()) {
    if (other.is_inf()) {
      inf = true;
    } else {
      // set EC parameters:
      mpz_set(a, other.get_a());
      mpz_set(b, other.get_b());
      mpz_set(n, other.get_n());

      // set coordinates:
      mpz_set(x, other.get_x());
      mpz_set(y, other.get_y());

      inf = false;
    }
    valid = true;
  } else {
    valid = false;
  }
}

ECPoint& ECPoint::operator =(const ECPoint& other) {
*/

void ECPoint::duplicate(ECPoint& other) {
  if (other.is_valid()) {
    if (other.is_inf()) {
      inf = true;
    } else {
      // set EC parameters:
      mpz_set(a, other.get_a());
      mpz_set(b, other.get_b());
      mpz_set(n, other.get_n());

      // set coordinates:
      mpz_set(x, other.get_x());
      mpz_set(y, other.get_y());

      inf = false;
    }
    valid = true;
  } else {
    valid = false;
  }
}

bool ECPoint::operator ==(ECPoint& other) { // compare P and Q:
  if (valid && other.is_valid()) { // both are valid:
    if (inf) { // P is infinite:
      return other.is_inf();
    } else if (!other.is_inf()) { // P is finite, if Q is also finite:
      return (mpz_cmp(x, other.get_x()) == 0)
          && (mpz_cmp(y, other.get_y()) == 0);
    } else { // P is finite, Q infinite:
      return false;
    }
  } else { // at least one is not valid, comparison fails:
    return false;
  }
}

void ECPoint::generate(mpz_t& A, mpz_t& X, mpz_t& Y, mpz_t& N) {
  mpz_set(a, A);
  mpz_set(x, X);
  mpz_set(y, Y);
  mpz_set(n, N);

  // compute b = y^2 - x^3 - ax:
  mpz_mul(b, x, x);    // b == x^2
  mpz_mul(b, b, x);    // b == x^3
  mpz_submul(b, y, y); // b == x^3 - y^2
  mpz_neg(b, b);       // b == y^2 - x^3
  mpz_submul(b, a, x); // b == y^2 - x^3 - ax
  mpz_mod(b, b, n);

  valid = true;
  inf = false;
}

bool ECPoint::add(mpz_t& remainder, ECPoint& P2_ref) {
  /* Add point P2 to this, report success in return value.
   * If unsuccessful, report the computed remainder.
   *
   * Assumes both points to be valid!
   * */

  static GMPNum gmp[3];
  ECPoint P2 {P2_ref};

  mpz_t &x2 = P2.get_x(),
        &y2 = P2.get_y();
  static mpz_t &slope = gmp[0],
               &tmp   = gmp[1],
               &tmp2  = gmp[2];



  if (is_inf()) { // if P1 is infinite, P1 + P2 = P2:
    duplicate(P2);
    return true;
  } else if (P2.is_inf()) { // if P2 is infinite, P1 + P2 = P1:
    return true;
  } else if (mpz_cmp(x, x2) == 0) { // P1.x == P2.x
    mpz_neg(tmp, y2);
    if (mpz_cmp(y, tmp) == 0) { // P1.y == -P2.y, therefore P1 == -P2:
      inf = true;
      return true;
    }
  }


  // Compute the slope:
  if ((mpz_cmp(x, x2) == 0) && mpz_cmp(y, y2) == 0) { // P1 == P2:
    // Compute slope = (3*x1^2 + a) / (2*y1):

    mpz_mul_ui(slope, y, 2);           // 2*y1
    mpz_mod(slope, slope, n);          // mod n
    if (!mpz_invert(tmp, slope, n)) {  // tmp = (2*y1)^(-1) mod n
      gcd(remainder, slope, n); // if the inversion failed, return the gcd.
      valid = false;
      return false;
    }
    
    mpz_mul(slope, x, x);              // x1^2
    mpz_mul_ui(slope, slope, 3);       // 3*x1^2
    mpz_add(slope, slope, a);          // 3*x1^2 + a
    mpz_mul(slope, slope, tmp);        // slope = (3*x1^2 + a) / (2*y1)

  } else {                                            // P1 != P2:
    // Compute slope = (3*x1^2 + a) / (2*y1):
    // Compute slope = (y2 - y1) / (x2 - x1):
    //
    mpz_sub(slope, x2, x);             // x2-x1
    mpz_mod(slope, slope, n);          // mod n
    if (!mpz_invert(tmp, slope, n)) {  // tmp = (x2-x1)^(-1) mod n
      gcd(remainder, slope, n); // if the inversion failed, return the gcd.
      valid = false;
      return false;
    }

    mpz_sub(slope, y2, y);             // y2-y1
    mpz_mul(slope, slope, tmp);        // (y2-y1) / (x2-x1)
  }

  // Slope is computed, compute the new coordinates:
  //   x = slope^2 - x1 - x2,
  //   y = slope*(x1 - x) - y1:

  mpz_set(tmp2, x);           // tmp2 = x1
  mpz_mul(tmp, slope, slope); // tmp = slope^2
  mpz_sub(tmp, tmp, x);       // tmp = slope^2 - x1
  mpz_sub(x, tmp, x2);        // x = slope^2 - x1 - x2

  mpz_sub(tmp, tmp2, x);      // x1 - x
  mpz_mul(tmp, slope, tmp);   // slope*(x1 - x)
  mpz_sub(y, tmp, y);         // y = slope*(x1 - x) - y1

  // Normalize:
  mpz_mod(x, x, n);
  mpz_mod(y, y, n);

 /*
  cout << "\n(";
  print_dec_mpz(x);
  cout << ", ";
  print_dec_mpz(y);
  cout << ")\n\n";
 */

  return true;
}

bool ECPoint::pow(mpz_t& remainder, mpz_t& exponent) {
  static GMPNum gmp[1];
  static mpz_t &e = gmp[0];

  ECPoint Y {true};

  mpz_set(e, exponent);

  while (mpz_cmp_ui(e, 1) > 0) {       // while e > 1:
    if (!mpz_divisible_2exp_p(e, 1)) { // if e%2 == 1:
      if (!Y.add(remainder, *this)) {  // Y += P
        valid = false; // if the addition failed, return the gcd
        return false;
      }
    }
    mpz_tdiv_q_2exp(e, e, 1);          // e /= 2
    if (!add(remainder, *this)) {      // P += P
      valid = false; // if the addition failed, return the gcd
      return false;
    }
  }

  if (!add(remainder, Y)) {            // P += Y
    valid = false; // if the addition failed, return the gcd
    return false;
  }

  return true;
}

bool ECPoint::pow(mpz_t& remainder, unsigned long e) {
  ECPoint Y {true}; // infinite

  while (e > 1) {                      // while e > 1:
    if (e&1) {                         // if e%2 == 1:
      if (!Y.add(remainder, *this)) {  // Y += P
        valid = false; // if the addition failed, return the gcd
        return false;
      }
    }
    e <<= 1;                           // e /= 2
    if (!add(remainder, *this)) {      // P += P
      valid = false; // if the addition failed, return the gcd
      return false;
    }
  }

  if (!add(remainder, Y)) {            // P += Y
    valid = false; // if the addition failed, return the gcd
    return false;
  }

  return true;
}



bool ECPoint::nonsingular(mpz_t& divisor) {
  static GMPNum gmp[2];
  static mpz_t &a3 = gmp[0],
               &b2 = gmp[1];

  mpz_mul(a3, a, a);
  mpz_mul(a3, a3, a);
  mpz_mul(b2, b, b);

  mpz_mul_ui(a3, a3, 4);
  mpz_mul_ui(b2, b2, 27);
  mpz_add(a3, a3, b2);

  gcd(divisor, a3, n);
  return mpz_cmp_ui(divisor, 1) == 0;
}














void test_curves() {
  static const char* N_STR =
    "0xffffffff00000001000000000000000000000000ffffffffffffffffffffffff";

  static const char* A_STR =
    "-0x3";

  static const char* B_STR = 
    "0x5ac635d8aa3a93e7b3ebbd55769886bc651d06b0cc53b0f63bce3c3e27d2604b";

  static const char* PX_STR =
    "0x6b17d1f2e12c4247f8bce6e563a440f277037d812deb33a0f4a13945d898c296";
  static const char* PY_STR =
    "0x4fe342e2fe1a7f9b8ee7eb4a7c0f9e162bce33576b315ececbb6406837bf51f5";

  static const char* KX_STR =
    "0x52910a011565810be90d03a299cb55851bab33236b7459b21db82b9f5c1874fe";
  static const char* KY_STR =
    "0xe3d03339f660528d511c2b1865bcdfd105490ffc4c597233dd2b2504ca42a562";

  GMPNum gmp[7];
  mpz_t &n  = gmp[0],
        &a  = gmp[1],
        &b  = gmp[2],
        &px = gmp[3],
        &py = gmp[4],
        &kx = gmp[5],
        &ky = gmp[6];


  mpz_set_str(n, N_STR, 0);
  mpz_set_str(a, A_STR, 0);
  mpz_set_str(b, B_STR, 0);
  mpz_set_str(px, PX_STR, 0);
  mpz_set_str(py, PY_STR, 0);
  mpz_set_str(kx, KX_STR, 0);
  mpz_set_str(ky, KY_STR, 0);

  mpz_mod(a, a, n);

//  print_hex_mpz(a);
//  cout << endl;

  ECPoint P, K, Q;

  P.generate(a, px, py, n);
  K.generate(a, kx, ky, n);

  // assignment test:
  Q = P;
  P = K;
  K = Q;

  Q = K;
  K = P;
  P = Q;

/*
  Q.set_inf();
  Q.add(a, P);
  print_hex_mpz(Q.get_x());
  cout << " ";
  print_hex_mpz(Q.get_y());
  cout << endl;
*/

/*
  if ((mpz_cmp(P.get_b(), b) != 0) || (mpz_cmp(K.get_b(), b) != 0)) {
    cerr << "b wasn't computed correctly." << endl;
    print_hex_mpz(b); cout << endl;
    print_hex_mpz(P.get_b()); cout << endl;
    print_hex_mpz(K.get_b()); cout << endl;
  }

  Q = P;
  if (!Q.add(a, K)) {
    cerr << "addition P+K failed for some reason" << endl;
  }
  if (mpz_cmp(Q.get_b(), b) != 0) {
    cerr << "b wasn't computed correctly." << endl;
  }
  cout << endl;
  print_dec_mpz(Q.get_x());
  cout << endl;
  print_dec_mpz(Q.get_y());
  cout << endl;

  Q = P;
  if (!Q.add(a, P)) {
    cerr << "addition P+P failed for some reason" << endl;
  }
  if (mpz_cmp(Q.get_b(), b) != 0) {
    cerr << "b wasn't computed correctly." << endl;
  }
  cout << endl;
  print_dec_mpz(Q.get_x());
  cout << endl;
  print_dec_mpz(Q.get_y());
  cout << endl;
*/

//  P.add(a, K);

///*
  mpz_set_ui(b, 69);
  P.pow(a, b);
  

  if (P == K) {
    cerr << "They are equal, great!" << endl;
  } else {
    cerr << "Damn, they aren't equal." << endl;
  }

  print_hex_mpz(P.get_x());
  cout << " ";
  print_hex_mpz(P.get_y());
  cout << endl;
//*/
/*
  print_hex_mpz(K.get_x());
  cout << " ";
  print_hex_mpz(K.get_y());
  cout << endl;
*/
}

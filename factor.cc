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

  // Attempt to find a factor using Brent factorization:
  if (brent_factorization(P, N)) {
    std::cerr << "Found using Brent" << std::endl;
    return;
  }

  // Attempt to find a factor using Lenstra EC factorization:
  if (lenstra_ec_factorization(P, N)) {
    std::cerr << "Found using Lenstra EC factorization" << std::endl;
    return;
  }

  // Attempt to find a factor using Brent factorization:
  if (brent_factorization(P, N)) {
    std::cerr << "Found using Brent" << std::endl;
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



bool brent_factorization(mpz_t& P, mpz_t& N) {
  Randomizer rnd;

  GMPNum gmp[11];
  mpz_t &Y   = gmp[0],
        &C   = gmp[1],
        &M   = gmp[2],
        &R   = gmp[3],
        &Q   = gmp[4],
        &K   = gmp[5],
        &X   = gmp[6],
        &YY  = gmp[7],
        &I   = gmp[8],
        &LIM = gmp[9],
        &tmp = gmp[10];

  if (mpz_divisible_2exp_p(N, 1)) { // if 2|N:
    mpz_set_ui(P, 2);
    return true;
  }



  rnd.get_random(Y, 1, N);
  rnd.get_random(C, 1, N);
  rnd.get_random(M, 1, N);

  std::cerr << "N: ";
  mpz_out_str(stderr, 10, N);
  std::cerr << std::endl; 

  std::cerr << "Y: ";
  mpz_out_str(stderr, 10, Y);
  std::cerr << std::endl; 

  std::cerr << "C: ";
  mpz_out_str(stderr, 10, C);
  std::cerr << std::endl; 

  std::cerr << "M: ";
  mpz_out_str(stderr, 10, M);
  std::cerr << std::endl; 


/*
  mpz_set_ui(Y, 2);
  mpz_set_ui(C, 1);
//  mpz_set_ui(M, 2);
  rnd.get_random(M, 1, N);
*/

  mpz_set_ui(P, 1);
  mpz_set_ui(R, 1);
  mpz_set_ui(Q, 1);

  while(mpz_cmp_ui(P, 1) == 0) { // while P == 1:
    mpz_set(X, Y);                 // X = Y
                                   // for I in [0, R):
    mpz_set_ui(I, 0);
    while (mpz_cmp(I, R) < 0) {
                                     // Y = Y^2 + C mod N:
      apply_f(Y, C, N, tmp);

      mpz_add_ui(I, I, 1);
    }
    mpz_set_ui(K, 0);              // K = 0
                                   // while K < R and G == 1:
    while ((mpz_cmp(K, R) < 0) && (mpz_cmp_ui(P, 1) == 0)) {
      mpz_set(YY, Y);                // YY = Y
                                     // for I in [0, min(M, R-K)):
      mpz_set_ui(I, 0);
      // LIM = min(M, R-K):
      mpz_sub(LIM, R, K);
      if (mpz_cmp(M, LIM) < 0) { mpz_set(LIM, M); }

      while (mpz_cmp(I, LIM) < 0) {
                                       // Y = Y^2 + C mod N:
        apply_f(Y, C, N, tmp);
                                       // Q *= |X-Y| mod N:
        mpz_sub(tmp, X, Y);
        mpz_abs(tmp, tmp);
        mpz_mul(Q, Q, tmp);
        mpz_mod(Q, Q, N);

        mpz_add_ui(I, I, 1);
      }
      gcd(P, Q, N);                  // P = gcd(Q, N)
      mpz_add(K, K, M);              // K += M
    }
    mpz_mul_2exp(R, R, 1);         // R *= 2
  }

  std::cerr << "R: ";
  mpz_out_str(stderr, 10, R);
  std::cerr << std::endl; 

  unsigned cnt = 0;
  if (mpz_cmp(P, N) == 0) {      // if P == N:
    std::cerr << "Running extra loop:\n";
    do {                           // repeat:
                                     // YY = YY^2 + C mod N:
      apply_f(YY, C, N, tmp);

                                     // tmp = |X-YY|:
      mpz_sub(tmp, X, YY);
      mpz_abs(tmp, tmp);
                                     // P = gcd(|X-YY|, N)
      gcd(P, tmp, N);
      cnt++;
    } while (mpz_cmp_ui(P, 1) <= 0); // until P > 1
    std::cerr << "Successfully finished extra loop after " << cnt
              << " iterations.\n";
  }

  return true;
}


bool brent_factorization_old(mpz_t& P, mpz_t& N) {

  GMPNum gmp[3];
  mpz_t &xi = gmp[0],
        &xm = gmp[1],
        &df = gmp[2]; 

  mpz_set_ui(xi, 2);
  mpz_set(xm, xi);
  unsigned long i = 1;

  if (mpz_divisible_2exp_p(N, 1)) { // if 2|N:
    mpz_set_ui(P, 2);
    return true;
  }

  while(true) {
    // TODO: co ak je prvocislo?
    mpz_mul(xi, xi, xi);
    mpz_add_ui(xi, xi, 1);
    mpz_mod(xi, xi, N);
    mpz_sub(df, xi, xm);

    gcd(P, df, N);
    if ((mpz_cmp_ui(P, 1) != 0) && (mpz_cmp(P, N) != 0)) {
      return true;
    }

    if (is2pow(i)) {
      mpz_set(xm, xi);
    }
    i++;
  }

  return false;
}



bool lenstra_ec_factorization(mpz_t& P, mpz_t& N) {
  /* It is assumed that N was already analyzed with a simpler method
   * and is not divisible by 2 or 3. */

  unsigned long B1 = B0;
  unsigned trial_cnt = 0;

  //GMPNum gmp[1];
  //mpz_t &k = gmp[0];

  // 0. Initial checks:

  // It is assumed that N is not divisible by 2 or 3.
  // Is n a perfect power of some integer?
  if (is_m_to_r(P, N, TRIAL_LIMIT)) {
    return true;
  }

  //std::cout << "Lenstra is happening." << std::endl;
  while (true) {

    // 1. Pick a curve:

    ECPoint E;
    if (!pick_curve(E, P, N)) {
      // A factor was found and written to P.
      return true;
    }

    // 2. Phase one - attempt to multiply E by small primes:
    if (ecm_phase_one(P, E, B1)) {
      return true;
    }


    // 3. Phase two - 
    if (ecm_phase_two(P, E, B1)) {
      return true;
    }

    trial_cnt += 1;
    B1 *= Q_TRIAL;                // increase B1

    //if ((trial_cnt & 255) == 255) {
      std::cerr << "Trial #" << trial_cnt << " finished, "
                << "new B1 = " << B1 << std::endl;
   // }
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
      m /= p;

      mpz_mul_ui(lcm, lcm, m);
      if ((p & 255) == 255) {
        mpz_mod(lcm, lcm, N);
      }
    }
  }
}

bool is_m_to_r(mpz_t& m, mpz_t& n, unsigned long limit) {
  /* Check whether n = m^r for integers m, r such that m >= limit.
   *  */
  GMPNum gmp;
  mpz_t& mr = gmp;

  /* GMP has a function for this, which is however not very useful, as it
   * doesn't provide the root:
   *
   *   return mpz_perfect_power_p(n);
   * */

  unsigned long r = 2;
  mpz_root(m, n, r);
  while (mpz_cmp_ui(m, limit) >= 0) {  // while m >= limit:
    mpz_pow_ui(mr, m, r);           // mr = m^r
    if (mpz_cmp(mr, n) == 0) {      // mr == n  -->  n == m^r:
      return true;
    }
    r++;
    mpz_root(m, n, r);
  }
  return false;
}

bool ecm_try_mult(mpz_t& remainder, ECPoint& E, unsigned long K) {
  Sieve& is_prime = Sieve::get_instance(K);

  unsigned m;

  for (unsigned p = 2; p < K; p++) {
    if (is_prime[p]) {
      m = p;
      while (m < K) {
        m *= p;
      }
      m /= p;

      if (! E.pow(remainder, m)) {
        return true;
      }
    }
  }
  return false;
}

bool ecm_phase_one(mpz_t& divisor, ECPoint& E, unsigned long B1) {
  Sieve& is_prime = Sieve::get_instance(B1);

  unsigned long m;

  for (unsigned long p = 2; p < B1; p++) {
    if (!is_prime[p]) {
      continue;
    }

    m = p;
    while (m < B1) { m *= p; }
    m /= p;

    if (!E.pow(divisor, m)) {
      return true;
    }

    if (E.is_inf()) {
      break;
    }
  }

  return false;
}

bool ecm_phase_two(mpz_t& divisor, ECPoint& E, unsigned long B1) {
  Randomizer rnd;

  unsigned long B2 = B1 * PHASE2_MULTIPLIER;
  unsigned e = choose_e(B1); // 1 <= e <= 60

  unsigned t = 5;
  while (static_cast<unsigned>((1 << t) - 1) <= (50*e)) {
    // 2^t - 1 <= 50*e:
    t++;
  }
  t--;

  unsigned r = (1 << t) - 1;
  unsigned rr = B2/r;

  std::vector<GMPNum> xi { r };

  GMPNum gmp[7];
  mpz_t &lim  = gmp[0],
        &tmp  = gmp[1],
        &prod = gmp[2],
        &X    = gmp[3],
        &tmp2 = gmp[4];

  ECPoint Q;

  mpz_set_ui(lim, (1<<30) / (e+2));
  unsigned long ui, vi, uj, vj;

  rnd.get_random(tmp, 1, lim);
  ui = mpz_get_ui(tmp);
  rnd.get_random(tmp, 1, lim);
  vi = mpz_get_ui(tmp);
  rnd.get_random(tmp, 1, lim);
  uj = mpz_get_ui(tmp);
  rnd.get_random(tmp, 1, lim);
  vj = mpz_get_ui(tmp);

  // compute polynomial:
  for (unsigned i = 0; i < r; i++) {
    mpz_set_ui(tmp, ui*i + vi);
    mpz_powm_ui(tmp, tmp, e, E.get_n());
    Q = E;
    if (!Q.pow(divisor, tmp)) {
      return true;
    }
    mpz_set(xi[i], Q.get_x());
  }

  // evaluate product:
  mpz_set_ui(prod, 1);
  for (unsigned j = 0; j < rr; j++) {
    mpz_set_ui(tmp, uj*j + vj);
    mpz_powm_ui(tmp, tmp, e, E.get_n());
    Q = E;
    if (!Q.pow(divisor, tmp)) {
      return true;
    }

    //evaluate polynomial for j:
    mpz_set(X, Q.get_x());
    mpz_set_ui(tmp, 1);
    for (unsigned i = 0; i < r; i++) {
      mpz_sub(tmp2, X, xi[i]);
      mpz_mul(tmp, tmp, tmp2);
    }
    mpz_mul(prod, prod, tmp);
  }

  gcd(divisor, prod, E.get_n());

  if ((mpz_cmp_ui(divisor, 1) != 0)
   && (mpz_cmp(divisor, E.get_n()) != 0)) {
    return true;
  }

  return false;
}

unsigned choose_e(unsigned long B1) {
  unsigned e = 1;
  unsigned i = 0;
  B1 /= 1250;

  while((e = E_SELECTION[++i]) != 0) {
    if ((e * e) > B1) {
      return E_SELECTION[i-1];
    }
  }
  return (e > 0 ? e : 1);
}

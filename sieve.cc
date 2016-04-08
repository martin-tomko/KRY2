#include "sieve.h"
#include <cmath>
#include <iostream>
#include <algorithm>

using namespace std;

void Sieve::resize(size_t new_size) {
  if (new_size > size) {
    size = new_size;
    prime.resize(size, true);
    compute();
  }
  // otherwise do nothing, new_size is a minimum requested size.
}


/* Compute Eratosthenes' sieve for the required size.
 *
 * A very basic implementation,
 * works fast enough for reasonable sizes (around 10 000 000),
 * but will get very slow very fast for larger values */
void Sieve::compute() {
  // The maximal value of a prime for which sieving is performed:
  size_t limit = std::floor (std::sqrt (size - 1));

  if (first_uncomputed <= 2) {
    // This is most likely the first call, initialize the array:
    prime[0] = prime[1] = false;
    first_uncomputed = 2;

    /* The rest of the vector should be initialized to true by the constructor
     * or resize(). */
  } else {
    // Some computation has already been performed, apply known primes to
    // new values:

    for (size_t p = 2; p < std::min(first_uncomputed, limit+1); p++) {
      if (!prime[p]) {
        // p is composite, do nothing.
        continue;
      }

      // The first multiple of p not yet computed:
      for (size_t k = (((first_uncomputed-1) / p) + 1) * p; k < size; k += p) {
        // mark all multiples of p as composite:
        prime[k] = false;

        /* The first multiple of p greater than last_value is l*p, where
         * l = ((first_uncomputed - 1) / p) + 1 */
      }
    }
  }


  for (size_t p = first_uncomputed; p <= limit; p++) {
    if (!prime[p]) { // p was already shown to be composite:
      continue;
    }

    for (size_t k = p*p; k < size; k += p) {
      // Mark all multiples of p as composite:
      prime[k] = false;

      /* p*p is the lowest unmarked multiple */
    } 
  }

  first_uncomputed = size;
}

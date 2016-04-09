/*
 * */

#ifndef __SIEVE_H
#define __SIEVE_H

#include <cstddef>
#include <vector>

using std::size_t;

const size_t DEFAULT_SIEVE_SIZE = 10000000;


class Sieve {
 public:
  ~Sieve() { prime.clear(); }

  static Sieve& get_instance()
    { static Sieve s {DEFAULT_SIEVE_SIZE}; return s; }

  static Sieve& get_instance(size_t min_size)
    { Sieve& s = get_instance(); s.resize(min_size); return s; }

  static void stretch(size_t min_size)
    { get_instance().resize(min_size); }

  bool operator[](size_t index) { return prime[index]; }
 private:
  Sieve() = delete;
  Sieve(size_t s) : prime (s, true), size {s}, first_uncomputed {0}
    { compute(); }

  void resize(size_t);     //
  void compute();          // Perform Eratosthenes' Sieve

  std::vector<bool> prime; // Vector containing boolean primality information
    /* Somehow wasteful, but std::bitset cannot be resized at runtime. */

  size_t size;             // Size of the vector
  size_t first_uncomputed; // The first index without a meaningful value
};

#endif

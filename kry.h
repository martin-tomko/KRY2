/*
 * */

#ifndef __KRY_H_
#define __KRY_H_

#include "gmp_helper.h"
#include <iostream>

extern const char* HELP_MSG;

enum class Mode {GENERATE, ENCRYPT, DECRYPT, FACTOR, HELP};
extern const char* MODE_NAME[];

struct Params {
  Mode  mode;
  GMPNum* val;
  size_t val_cnt;

  Params() : val {NULL}, val_cnt {0} { }
  ~Params() { delete[] val; }
  void set(size_t cnt)
    { if (val) delete[] val; val_cnt = cnt; val = new GMPNum[cnt]; }
};

void get_params(int&, char**&, Params&);
inline void print_help() { std::cerr << HELP_MSG; }
void print_params(Params&);
inline void print_hex_mpz(mpz_t& n)
  { std::cout << "0x"; mpz_out_str(stdout, 16, n); }
inline void print_dec_mpz(mpz_t& n)
  { mpz_out_str(stdout, 10, n); }

#endif

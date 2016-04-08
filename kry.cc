/*
 * */

#include "kry.h"
#include <iostream>
#include <cstring>
#include <stdexcept>
#include <string>
#include <sstream>

#include "factor.h"

#include "sieve.h"

using namespace std;

const char* HELP_MSG =
"usage:\n"
"  ./kry -g B     -- generate a B-bit RSA key\n"
"  ./kry -e E N M -- encrypt M using RSA with (E,N) as a key\n"
"  ./kry -d D N C -- decrypt C using RSA with (D,N) as a key\n"
"  ./kry -b N     -- find a prime factor of N\n";

const char* MODE_NAME[] = {
  "generate", "encrypt", "decrypt", "factor", "help"
};

int main(int argc, char** argv) {
  Params par;

  try {
    get_params(argc, argv, par);
//    print_params(par);
  } catch (runtime_error& e) {
    cerr << "Exception: " << e.what() << endl;
    print_help();
    return 1;
  }


  try {
    switch(par.mode) {
      case Mode::GENERATE: // generate RSA key:
        throw runtime_error("GENERATE not yet implemented!");
        break;

      case Mode::ENCRYPT: { // Encrypt message:
        GMPNum gmp_C;
        mpz_t& C = gmp_C;           // ciphertext
        mpz_t& E = par.val[0];      // public exponent
        mpz_t& N = par.val[1];      // RSA module
        mpz_t& M = par.val[2];      // plaintext

        mpz_powm(C, M, E, N);       // c = m^e mod n

        print_hex_mpz(C);           // print result to stdout
        cout << endl;
      } break;

      case Mode::DECRYPT: { // Decrypt message:
        GMPNum gmp_M;
        mpz_t& M = gmp_M;           // plaintext
        mpz_t& D = par.val[0];      // public exponent
        mpz_t& N = par.val[1];      // RSA module
        mpz_t& C = par.val[2];      // ciphertext

        /* Does not decrypt using a secure version of powm and is therefore
         * susceptible to side-channel attacks.
         * Newer versions of GMP contain the function "mpz_powm_sec" which
         * would be better suited to this application. */
        mpz_powm(M, C, D, N);       // m = c^d mod n

        print_hex_mpz(M);           // print result to stdout
        cout << endl;
      } break;

      case Mode::FACTOR: { // find a factor of a number:
        // throw runtime_error("FACTOR not yet implemented!");
        GMPNum gmp_P;
        mpz_t& P = gmp_P;           // prime factor
        mpz_t& N = par.val[0];      // number to be factored

        factor(P, N);               // factor the number

        print_hex_mpz(P);           // print the factor to stdout
        cout << endl;
      } break;

      case Mode::HELP: // print help msg:
        print_help();
        return 1;
        break;

      default:
        throw runtime_error("Invalid program mode (default branch reached).");
        break;
    }
  } catch (runtime_error& e) {
    cerr << "Exception: " << e.what() << endl;
    return 1;
  }


  return 0;
}

void get_params(int& argc, char**& argv, Params& par) {
  int expected_argc = -1;

  if (argc <= 1) {
    throw runtime_error("Invalid arguments.");
  } else if (strcmp(argv[1], "-g") == 0) { // generate:
    par.mode = Mode::GENERATE;
    expected_argc = 3;
  } else if (strcmp(argv[1], "-e") == 0) { // encrypt:
    par.mode = Mode::ENCRYPT;
    expected_argc = 5;
  } else if (strcmp(argv[1], "-d") == 0) { // decrypt:
    par.mode = Mode::DECRYPT;
    expected_argc = 5;
  } else if (strcmp(argv[1], "-b") == 0) { // factor:
    par.mode = Mode::FACTOR;
    expected_argc = 3;
  } else {
    throw runtime_error("Invalid arguments.");
  }

  if (argc != expected_argc) {
    throw runtime_error("Invalid number of parameters.");
  }

  // argument count is valid, read the numbers:
  // TODO: prerobit na cyklus (spolu s upravou struct Params)
  par.set(expected_argc - 2);
  for (size_t i = 0; i < par.val_cnt; i++) {
    if (mpz_set_str(par.val[i], argv[i+2], 0) != 0) {
      ostringstream oss("Parameter ");
      oss << i << " not a valid number.";
      throw runtime_error(oss.str());
    }
  }
}

void print_params(Params& par) {
  cerr << "Mode: " << MODE_NAME[static_cast<size_t>(par.mode)] << endl;
  cerr << "Value count: " << par.val_cnt << endl;
  for (size_t i = 0; i < par.val_cnt; i++) {
    cerr << "val[" << i << "] = ";
    mpz_out_str(stderr, 10, par.val[i]);
    cerr << endl;
  }
}

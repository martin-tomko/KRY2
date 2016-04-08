/*
 * */

#include "gmp_helper.h"

void Randomizer::seed() {
  // TODO: get seed somehow
  //   1. /dev/random, if not avaliable, then:
  //   2. system time, if not available, then:
  //   3. whatever
  unsigned long seed = 48321841; // TODO: improve
  gmp_randseed_ui(state, seed); 
}

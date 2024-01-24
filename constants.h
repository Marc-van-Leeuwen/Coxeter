/*
  This is constants.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef CONSTANTS_H  /* guard against multiple inclusions */
#define CONSTANTS_H

#include "globals.h"
#include <limits.h>

namespace constants {
  using namespace globals;
};

#define BITS(x) (CHAR_BIT*sizeof(x))  /* size of x in bits */

namespace constants {

  const Ulong CHARFLAGS = ((Ulong)1 << CHAR_BIT)-1; // |0xFF| as |Ulong|

  extern Ulong *lmask; // bit posistion |==i| masks, for |i| up to |BITS(Ulong)|
  extern Ulong *leqmask; // same for bit position |<=i| masks
  extern unsigned *firstbit; // |2^CHAR_BIT| "lowest set bit" positions
  extern unsigned *lastbit;

  unsigned firstBit(Ulong f);
  void initConstants(); // set up pointers above, pointing into static arrays
  unsigned lastBit(Ulong f);
};

#endif

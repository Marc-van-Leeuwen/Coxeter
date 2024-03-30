/*
  This is constants.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef CONSTANTS_H  /* guard against multiple inclusions */
#define CONSTANTS_H

#include "globals.h"
#include <limits.h> // for |CHAR_BITS|, the number of bits per |char|
#include <limits>

#define BITS(x) (CHAR_BIT*sizeof(x))  /* size of x in bits */

namespace constants {

  namespace aux {
    template<unsigned long n> struct BaseShift // value gives 2log of n
    { static constexpr unsigned long value = 1 + BaseShift<n/2>::value; };
    template<> struct BaseShift<1ul>
    { static constexpr unsigned long value = 0; };
  }

  static constexpr unsigned longBits = // typically 64
    std::numeric_limits<Ulong>::digits;
  static const unsigned posBits = // mask (often widened), typically |0x3F==63|
    longBits - 1; // part of bit address for position within ull word
  static const unsigned long baseBits = // typically |0xFFFFFFFFFFFFFFC0|
    ~static_cast<Ulong>(posBits); // mask for non bit-positions in a bit address
  static constexpr unsigned long baseShift = // typically 6
    aux::BaseShift<longBits>::value;

  const Ulong CHARFLAGS = ((Ulong)1 << CHAR_BIT)-1; // |0xFF| as |Ulong|

  extern Ulong *eq_mask; // bit position |==i| masks, for |i<BITS(Ulong)|
  extern Ulong *lt_mask; // bit position |<i| masks, for |i<=BITS(Ulong)|
  extern Ulong *leq_mask; // bit position |<=i| masks, for |i<BITS(Ulong)|
  extern unsigned *firstbit; // |2^CHAR_BIT| "lowest set bit" positions
  extern unsigned *lastbit; // |2^CHAR_BIT| "highest set bit" positions

  void initConstants(); // set up pointers above, pointing into static arrays
  unsigned first_bit(Ulong f); // position first set bit, assumed to exist
  unsigned last_bit(Ulong f); // position last set bit, assumed to exist
  unsigned firstBit(Ulong f); // first set bit, or |BITS(Ulong)| if none
  unsigned lastBit(Ulong f); // last set bit, or |BITS(Ulong)| if none

};

#endif

/*
  This is constants.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "constants.h"

#include <cassert>

namespace constants {
  Ulong *eq_mask;
  Ulong *lt_mask;
  Ulong *leq_mask;
  unsigned *firstbit;
  unsigned *lastbit;
};

/****************************************************************************

  This module provides some basic constants that will be used elsewhere
  in the program. The idea is that the call to initConstants() should be
  the very first one in the program; therefore this module should not
  depend on anything else.

 ****************************************************************************/

namespace constants {

void initConstants()

{
  static Ulong d_eq_mask[BITS(Ulong)];
  static Ulong d_lt_mask[BITS(Ulong)+1];

  eq_mask = &d_eq_mask[0]; // point pointer to that static array
  lt_mask = &d_lt_mask[0];
  leq_mask = &d_lt_mask[1];

  lt_mask[0] = 0L;

  for (Ulong j = 0; j < BITS(Ulong); j++)
    {
      eq_mask[j] = 1ul << j;
      leq_mask[j] = lt_mask[j] | eq_mask[j];
    }

  static unsigned d_firstbit[1<<CHAR_BIT];
  static unsigned d_lastbit[1<<CHAR_BIT];
  firstbit = d_firstbit; // point pointer to that static array
  lastbit = d_lastbit; // point pointer to that static array

  d_firstbit[0] = CHAR_BIT; // "out of range" value: no such set bit
  d_firstbit[0] = 0; // basis for recursive propagation below

  d_lastbit[0] = CHAR_BIT; // "out of range" value: no such set bit
  d_lastbit[1] = 0; // basis for recursive propagation below

  for (unsigned j = 1; j < (1 << (CHAR_BIT-1)); ++j)
  {
    d_firstbit[2*j]   = d_firstbit[j]+1; // even number is its half, shifted
    d_firstbit[2*j+1] = 0; // first set bit for any odd number is 1
    d_lastbit[2*j]   = d_lastbit[j]+1; // even number is its half, shifted
    d_lastbit[2*j+1] = d_lastbit[j]+1; // odd number has same last bit
  }

}


// get first 'set' bit position in cases where it is known to exist
unsigned first_bit(Ulong f) // bit position of the first set bit in |f|
{
  assert(f!=0);
  unsigned shift = 0;
  if ((f&0xFFFFFFFF)==0)
    shift=32,f>>=32;
  if ((f&0xFFFF)==0)
    shift+=16,f>>=16;
  if ((f&0xFF)==0)
    shift+=8,f>>=8;
  return shift + firstbit[f&0xFF];
}

unsigned firstBit(Ulong f) // bit position of the first set bit in |f|
{
  if (f == 0)
    return BITS(Ulong); // "out of range" value: no such set bit

  return first_bit(f);
}


unsigned lastBit(Ulong f) // bit position of the first set bit in |f|
{
  if (f == 0)
    return BITS(Ulong); // "out of range" value: no such set bit

  unsigned shifted = 0;
  while ((f&~CHARFLAGS)!=0)
  {
    shifted += CHAR_BIT;
    f >>= CHAR_BIT;
  }

  return shifted + lastbit[f];
}

};

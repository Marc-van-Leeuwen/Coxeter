/*
  This is constants.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "constants.h"

namespace constants {
  Ulong *lmask;
  Ulong *leqmask;
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
  static Ulong d_lmask[BITS(Ulong)];
  static Ulong d_leqmask[BITS(Ulong)];

  lmask = d_lmask; // point pointer to that static array
  leqmask = d_leqmask;

  leqmask[0] = 1L;
  lmask[0] = 1L;

  for (Ulong j = 1; j < BITS(Ulong); j++)
    {
      lmask[j] = lmask[j-1] << 1;
      leqmask[j] = leqmask[j-1] | lmask[j];
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

unsigned firstBit(Ulong f) // bit position of the first set bit in |f|
{
  if (f == 0)
    return BITS(Ulong); // "out of range" value: no such set bit

  unsigned shifted = 0;
  while ((f&CHARFLAGS)==0)
  {
    shifted += CHAR_BIT;
    f >>= CHAR_BIT;
  }
  return shifted + firstbit[f&CHARFLAGS];
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

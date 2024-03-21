/*
  This is globals.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef GLOBALS_H  /* guard against multiple inclusions */
#define GLOBALS_H

#include <stddef.h>
#include <stdlib.h>
#include <cstdint> // for |uint32_t| and |uint64_t|

// these definitions are in the global namspace
// so it suffices to #include this file to have them
using Ulong = std::uint64_t; // used to be |unsigned long|, in practice is same

using Lflags = std::uint64_t; // can hold two bits per Coxeter genertor
using GenSet = std::uint32_t; // can hold one bit per Coxeter generator

#endif

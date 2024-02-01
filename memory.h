/*
  This is memory.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef MEMORY_H  /* guard against multiple inclusions */
#define MEMORY_H

#include <stdio.h> // for |FILE|
#include "globals.h"
#include "constants.h"

/******** type declarations *************************************************/

namespace memory {
  union Align;
  class Arena;
};

/******** function declarations *********************************************/

void* operator new (size_t size, memory::Arena& a);
void* operator new[] (size_t size, memory::Arena& a);

namespace memory {
  Arena& arena(); // get a reference to our unique static root |Arena| object
  void pause(); // do nothing; presumably for setting break point in debugging
};

/******** Type definitions **************************************************/

union memory::Align {
  Ulong d_ulong;
  void *d_voidptr;
};

class memory::Arena
{
  struct MemBlock { MemBlock *next; }; // basis, actual memory block follows

  // the following data are separated by size |2^i| of blocks; |x[i]| is:
  MemBlock* d_list[BITS(Ulong)]; // linked list of free size |2^i| blocks
  Ulong d_used[BITS(Ulong)]; // number of size |2^i| blocks in actual use
  Ulong d_allocated[BITS(Ulong)]; // total number of size |2^i| blocks allocated

  unsigned d_bsBits; // minimal power of 2 for which we buy memory chunks
  unsigned d_count; // total number of words that were ever |calloc|ed
  void newBlock(unsigned b); // ensure a bloick of size |2^b| is available

 public:
/* constructors and destructors */
  Arena(Ulong bsBits); // create empty |Arena|; |bsBits| only sets its |d_bsBits|
  ~Arena();

/* modifiers */
  void *alloc(size_t n);
  void *realloc(void *ptr, size_t old_size, size_t new_size);
  void free(void *ptr, size_t n);

/* accessors */
  static Ulong allocSize(Ulong n, Ulong m);
  static Ulong byteSize(Ulong n, Ulong m);
  void print(FILE* file) const;
};

/******** Inline implementations *****************************************/

inline void* operator new(size_t size, memory::Arena& a)
  {return a.alloc(size);}
inline void* operator new[](size_t size, memory::Arena& a)
  {return a.alloc(size);}

#endif

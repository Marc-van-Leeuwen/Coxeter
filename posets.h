/*
  This is posets.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef POSETS_H  /* guard against multiple inclusions */
#define POSETS_H

#include "globals.h"

/******** type declarations *************************************************/

namespace posets {
  typedef Ulong PosetElt;
  class Poset;
};

/******** type definitions **************************************************/

#include "bits.h"
#include "list.h"
#include "memory.h"
#include "wgraph.h"

namespace posets {

class Poset {
  containers::vector<bitmap::BitMap> d_closure;
 public:
/* constructors and destructors */
  void operator delete(void* ptr, size_t size)
    {return memory::arena().free(ptr,sizeof(Poset));}
  Poset();
  Poset(const Ulong &n);
  Poset(const wgraph::OrientedGraph& G);
/* manipulators */
/* accessors */
  void findMaximals(const bits::BitMap& D, bits::Set& a) const;
  containers::vector<Ulong> maxima_within(bitmap::BitMap D) const;
  bool isTriangular() const;
  Ulong size() const;
  void hasseDiagram(wgraph::OrientedGraph& H) const;
/* input/output */
};

};

/******** inline implementations ********************************************/

namespace posets {

inline Ulong Poset::size() const {return d_closure.size();}

};

#endif

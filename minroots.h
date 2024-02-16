/*
  This is minroots.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef MINROOTS_H  /* guarantee single inclusion */
#define MINROOTS_H

#include <limits.h>
#include "globals.h"
#include "graph.h"
#include "list.h"
#include "memory.h"

/******** type declarations *************************************************/

namespace minroots {
  typedef unsigned MinNbr;
  typedef char DotProduct;
  class MinTable;
};

/* constants */

namespace minroots {
  const MinNbr MINNBR_MAX = UINT_MAX-4;  /* top values are reserved */
  const MinNbr MINROOT_MAX = MINNBR_MAX; /* should not exceed MINNBR_MAX */
  const MinNbr undef_minnbr = MINNBR_MAX + 1;
  const MinNbr not_minimal = MINNBR_MAX + 2;
  const MinNbr not_positive = MINNBR_MAX + 3;
};

/******** function declarations *********************************************/

#include "bits.h"
#include "coxtypes.h"
#include "dotval.h"
#include "io.h"

namespace minroots {
  std::string& append(std::string& str, const dotval::DotVal& a);
  bits::Lflags descent(MinTable& T, MinNbr r);
  coxtypes::Length depth(MinTable& T, MinNbr r);
  void print(FILE *file, MinTable& T);
  coxtypes::CoxWord& reduced(MinTable& T, MinNbr r);
  bits::Lflags support(MinTable& T, MinNbr r);
};

/******* type definitions ****************************************************/


class minroots::MinTable {
 protected:
  coxtypes::Rank d_rank;
  MinNbr d_size;
  list::List<MinNbr*> d_min;
  list::List<DotProduct*> d_dot;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
  {return memory::arena().free(ptr,sizeof(MinTable));}
  void* operator new(size_t, void* ptr) {return ptr;}
  void operator delete(void* ptr, void* placement) {};
  MinTable() {};
  MinTable(const graph::CoxGraph& G);
  ~MinTable();
/* manipulators */
  void fill(const graph::CoxGraph& G);
/* accessors */
  bits::Lflags descent(const coxtypes::CoxWord& g) const;
  dotval::DotVal dot(MinNbr r, coxtypes::Generator s) const
    { return dotval::DotVal(d_dot[r][s]); }
  int insert
    (coxtypes::CoxWord& g, const coxtypes::Generator& s,
     const bits::Permutation& order) const;
  const coxtypes::CoxWord& inverse(coxtypes::CoxWord& g) const;
  bool inOrder(const coxtypes::CoxWord& g, const coxtypes::CoxWord& h) const;
  bool inOrder
    (list::List<coxtypes::Length>& a,
     const coxtypes::CoxWord& g, const coxtypes::CoxWord& h)
    const;
  bool isDescent(const coxtypes::CoxWord& g, const coxtypes::Generator& s) const;
  bits::Lflags ldescent(const coxtypes::CoxWord& g) const;
  const coxtypes::CoxWord& normalForm
    (coxtypes::CoxWord& g, const bits::Permutation& order) const;
  MinNbr min(MinNbr r, coxtypes::Generator s) const { return d_min[r][s]; }
  int prod(coxtypes::CoxWord& g, const coxtypes::Generator& s) const;
  int prod
    (coxtypes::CoxWord& g, coxtypes::CoxLetter *const h, const Ulong& n) const;
  int prod(coxtypes::CoxWord& g, const coxtypes::CoxWord& h) const;
  coxtypes::Rank rank() const { return d_rank; }
  bits::Lflags rdescent(const coxtypes::CoxWord& g) const;
  const coxtypes::CoxWord& reduced
    (coxtypes::CoxWord& g, coxtypes::CoxWord& h) const;
  MinNbr size() const { return d_size; }
  const coxtypes::CoxWord& power(coxtypes::CoxWord& a, const Ulong& m) const;
}; // |class minroots::MinTable|

#endif

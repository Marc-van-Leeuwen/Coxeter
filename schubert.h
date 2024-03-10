/*
  This is schubert.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef SCHUBERT_H  /* guard against multiple inclusions */
#define SCHUBERT_H

#include "globals.h"
#include "coxtypes.h"
#include "graph.h"
#include "containers.h"
#include "bitmap.h"
#include "sl_list.h"
#include "io.h"
#include "list.h"
#include "stack.h"
#include "bits.h"
#include "interface.h"

/******** type declarations *************************************************/

namespace schubert {
  class ClosureIterator;
  class SchubertContext;
  class StandardSchubertContext;

  struct NFCompare;

  using CoxNbrList = containers::vector<coxtypes::CoxNbr>;
  typedef list::List<coxtypes::BettiNbr> Homology;
};

/******** function declarations *********************************************/


namespace schubert {
  void betti(Homology& h, coxtypes::CoxNbr y, const SchubertContext& p);
  void extractInvolutions(const SchubertContext& p, bits::BitMap& b);
  void extractMaximals
    (const SchubertContext& p, list::List<coxtypes::CoxNbr>& c);
  void extractMaximals
    (const SchubertContext& p, list::List<coxtypes::CoxNbr>& c,
     list::List<Ulong>& a);
  void select_maxima_for
    (const SchubertContext& p, bits::BitMap& b, const Lflags& f);
  void select_maxima_for
    (const SchubertContext& p, bitmap::BitMap& b, Lflags f);
  Ulong min(const bits::Set& c, NFCompare& nfc);
  Ulong minDescent(const GenSet& f, const bits::Permutation& order);
  void minimize
    (const SchubertContext& p, bits::BitMap& b, const Lflags& f);
  void print(FILE* file, const SchubertContext& p);
  void printBitMap(FILE* file, const bits::BitMap& pi, const SchubertContext& p,
		   const interface::Interface& I);
  void printList
   (FILE* file, const list::List<coxtypes::CoxNbr>& v, const SchubertContext& p,
    const interface::Interface& I);
  void printPartition
    (FILE* file, const bits::Partition& pi,
     const SchubertContext& p, const interface::Interface& I);
  void printPartition
    (FILE* file, const bits::Partition& pi, const bits::BitMap& b,
     const SchubertContext& p, const interface::Interface& I);
  void readBitMap(list::List<coxtypes::CoxNbr>& c, const bits::BitMap& b);
  void read_bitmap
    (CoxNbrList& c, const bits::BitMap& b);
  bool shortLexOrder(const SchubertContext& p, coxtypes::CoxNbr x,
		     coxtypes::CoxNbr y, const bits::Permutation& order);
  void setBitMap(bits::BitMap& b, const list::List<coxtypes::CoxNbr>& c);
  Ulong sum(const Homology& h);
};

/******** type definitions *************************************************/

namespace schubert {

struct NFCompare {
  const SchubertContext& p;
  const bits::Permutation& order;
  NFCompare(const SchubertContext& q,
	    const bits::Permutation& generator_ordering)
    :p(q),order(generator_ordering) {};
  ~NFCompare() {};
  bool operator()(coxtypes::CoxNbr x, coxtypes::CoxNbr y)
  {return shortLexOrder(p,x,y,order);}
};

class AbstractSchubertContext {
  friend class ClosureIterator;
 public:
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(AbstractSchubertContext));}
  virtual ~AbstractSchubertContext() {};
/* accessors */
  virtual coxtypes::CoxWord& append
    (coxtypes::CoxWord& g, coxtypes::CoxNbr x) const = 0;
  virtual Lflags ascent(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::CoxNbr contextNumber(const coxtypes::CoxWord& g) const = 0;
  virtual Lflags descent(coxtypes::CoxNbr x) const = 0;
  virtual const bits::BitMap& downset(coxtypes::Generator s) const = 0;
  virtual void extendSubSet(bits::SubSet& q, coxtypes::Generator s) const = 0;
  virtual void extractClosure(bits::BitMap& b, coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::Generator firstDescent(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::Generator firstLDescent(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::Generator firstRDescent(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::Generator firstDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const = 0;
  virtual coxtypes::Generator firstLDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const = 0;
  virtual coxtypes::Generator firstRDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const = 0;
  virtual const CoxNbrList& hasse(coxtypes::CoxNbr x) const = 0;
  virtual bool inOrder(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const = 0;
  virtual bool isDescent(coxtypes::CoxNbr x, coxtypes::Generator s) const = 0;
  virtual GenSet lascent(coxtypes::CoxNbr x) const = 0;
  virtual GenSet ldescent(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::Length length(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::CoxNbr lshift
    (coxtypes::CoxNbr x, coxtypes::Generator s) const = 0;
  virtual coxtypes::CoxNbr maximize
    (coxtypes::CoxNbr x, const Lflags& f) const = 0;
  virtual coxtypes::Length maxlength() const = 0;
  virtual coxtypes::CoxNbr minimize
    (coxtypes::CoxNbr x, const Lflags& f) const = 0;
  virtual coxtypes::CoxWord& normalForm
    (coxtypes::CoxWord& g, coxtypes::CoxNbr x,
     const bits::Permutation& order) const = 0;
  virtual Ulong nStarOps() const = 0;
  virtual const bitmap::BitMap& parity(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::Rank rank() const = 0;
  virtual GenSet rascent(coxtypes::CoxNbr x) const = 0;
  virtual GenSet rdescent(coxtypes::CoxNbr x) const = 0;
  virtual coxtypes::CoxNbr rshift
    (coxtypes::CoxNbr x, coxtypes::Generator s) const = 0;
  virtual GenSet S() const = 0;
  virtual coxtypes::CoxNbr shift
    (coxtypes::CoxNbr x, coxtypes::Generator s) const = 0;
  virtual coxtypes::CoxNbr size() const = 0;
  virtual coxtypes::CoxNbr star
    (coxtypes::CoxNbr x, const Ulong& r) const = 0;
  virtual Lflags twoDescent(coxtypes::CoxNbr x) const = 0;
  virtual const type::Type& type() const = 0;
/* modifiers */
  virtual coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g) = 0;
  virtual void permute(const bits::Permutation& a) = 0;
  virtual void revertSize(const Ulong& n) = 0;
  virtual void setSize(const Ulong& n) = 0;
/* input-output */
  virtual std::string& append(std::string&, coxtypes::CoxNbr x) const = 0;
  virtual std::string& append
    (std::string&, coxtypes::CoxNbr x, const interface::Interface& I) const = 0;
  virtual void print(FILE* file, coxtypes::CoxNbr x) const = 0;
  virtual void print
    (FILE* file, coxtypes::CoxNbr x, const interface::Interface& I) const = 0;
};

class SchubertContext
{
 private:
/* private class declaration */
  class ContextExtension;
  const graph::CoxGraph& d_graph;
  coxtypes::Rank d_rank;
  coxtypes::Length d_maxlength;
  coxtypes::CoxNbr d_size;
  containers::vector<coxtypes::Length> d_length;
  containers::vector<CoxNbrList> d_hasse;
  containers::vector<Lflags> d_descent;
  containers::matrix<coxtypes::CoxNbr> d_shift;
  list::List<coxtypes::CoxNbr*> d_star; // indexed by |CoxNbr|, then |Ulong|
  containers::vector<bitmap::BitMap> d_downset; // length |2*d_rank|
  bitmap::BitMap d_parity[2]; // array of TWO parity bitmaps
  containers::stack<ContextExtension> d_history;
/* private member functions */
  void fillCoatoms(const Ulong& first, coxtypes::Generator s);
  void fillDihedralShifts(coxtypes::CoxNbr x, coxtypes::Generator s);
  void fillShifts(coxtypes::CoxNbr first, coxtypes::Generator s);
  void fillStar(coxtypes::CoxNbr first);
  void fullExtension(bits::SubSet& q, coxtypes::Generator s);
  void subSetExtension(bits::SubSet& q, coxtypes::Generator s) const;
 public:
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(SchubertContext));}
/* constructors and destructors */
  SchubertContext(const graph::CoxGraph& G);
  ~SchubertContext();

  // inlined accessor methods
  coxtypes::Rank rank() const { return d_rank; }
  coxtypes::CoxNbr size() const { return d_size; }
  coxtypes::Length length(coxtypes::CoxNbr x) const { return d_length[x]; }
  const bitmap::BitMap& parity(coxtypes::CoxNbr x) const
    { return d_parity[d_length[x]%2]; }
  coxtypes::Length maxlength() const { return d_maxlength; }

  Lflags descent(coxtypes::CoxNbr x) const
     {return d_descent[x];}
  Lflags ascent(coxtypes::CoxNbr x) const
    { return ~d_descent[x]&constants::lt_mask[2*d_rank]; }
  GenSet rdescent(coxtypes::CoxNbr x) const
    {return d_descent[x] & constants::lt_mask[d_rank];}
  GenSet rascent(coxtypes::CoxNbr x) const
    {return ~rdescent(x)&constants::lt_mask[d_rank];}
  GenSet ldescent(coxtypes::CoxNbr x) const
    { return d_descent[x] >> d_rank; } // left descents as (neutral) generators
  GenSet lascent(coxtypes::CoxNbr x) const
    { return ~ldescent(x)&constants::lt_mask[d_rank]; }
  coxtypes::Generator firstDescent(coxtypes::CoxNbr x) const
    { return constants::firstBit(descent(x));}
  coxtypes::Generator firstLDescent(coxtypes::CoxNbr x) const
    { return constants::firstBit(ldescent(x)); }
  coxtypes::Generator firstRDescent(coxtypes::CoxNbr x) const
    { return constants::firstBit(rdescent(x));}
  coxtypes::Generator firstDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const
    { return firstRDescent(x,order); }
  coxtypes::Generator firstLDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const
    { return minDescent(ldescent(x),order); }
  coxtypes::Generator firstRDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const
    { return minDescent(rdescent(x),order); }

  bool isDescent(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { return (d_descent[x]&constants::eq_mask[s])!=0; } // whether Right descent

  coxtypes::CoxNbr shift(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { return d_shift.entry(x,s); } // left or right shift
  coxtypes::CoxNbr rshift(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { assert(s<d_rank); return d_shift.entry(x,s); }
  coxtypes::CoxNbr lshift(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { assert(s<d_rank); return d_shift.entry(x,d_rank+s); }

  const bitmap::BitMap& down_set(coxtypes::Generator s) const
    { return d_downset[s]; }
  bits::BitMap downset(coxtypes::Generator s) const
  { return bits::BitMap(down_set(s)); } // convert
  const CoxNbrList& hasse(coxtypes::CoxNbr x) const { return d_hasse[x]; }
  const type::Type& type() const { return d_graph.type(); }

  Lflags S() const // mask for simple generators, to then restrict to |GenSet|
    { return constants::lt_mask[d_rank]; }

/* accessors */
  bitmap::BitMap closure(coxtypes::CoxNbr x) const;
  containers::sl_list<coxtypes::Generator> word(coxtypes::CoxNbr x) const;
  void spread_subset
    (bitmap::BitMap& q, CoxNbrList& elements, coxtypes::Generator s) const;

  coxtypes::CoxWord& append(coxtypes::CoxWord& g, coxtypes::CoxNbr x) const;
  coxtypes::CoxNbr contextNumber(const coxtypes::CoxWord& g) const;
  void extendSubSet(bits::SubSet& q, coxtypes::Generator s) const;
  void extractClosure(bits::BitMap& b, coxtypes::CoxNbr x) const;
  Lflags twoDescent(coxtypes::CoxNbr x) const;
  bool inOrder(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const;
  bool isSuperExtremal(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const;
  coxtypes::CoxNbr maximize(coxtypes::CoxNbr x, const Lflags& f) const;
  coxtypes::CoxNbr minimize(coxtypes::CoxNbr x, const Lflags& f) const;
  coxtypes::CoxWord& normalForm
    (coxtypes::CoxWord& g, coxtypes::CoxNbr x, const bits::Permutation& order)
    const;
  Ulong nStarOps() const { return d_graph.finite_edges().size(); }
  coxtypes::CoxNbr star(coxtypes::CoxNbr x, const Ulong& r) const
    { return d_star[x][r]; }
/* manipulators */
  coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);
  void increase_size(const Ulong& n);
/* i/o */
  std::string& append(std::string&, coxtypes::CoxNbr x) const;
  std::string& append(std::string&, coxtypes::CoxNbr x,
		      const interface::Interface& I) const;
  void print(FILE* file, coxtypes::CoxNbr x) const;
  void print(FILE* file, coxtypes::CoxNbr x, const interface::Interface& I)
    const;

}; // |class SchubertContext|


class ClosureIterator {
 private:
  struct node
  { bitmap::BitMap closure; // the Bruhat interval $[e,current]$
    coxtypes::CoxNbr closure_size; // the size of that interval (for speed)
    coxtypes::CoxNbr current; // the current Coxeter element
    GenSet asc; // its not yet explored ascent set
  };
  const SchubertContext& d_schubert;
  containers::stack<node> state;
  CoxNbrList elements; // closure, for fast traversal
  bitmap::BitMap d_visited; // everything seen so far
 public:
/* constructors and destructors */
  ClosureIterator(const SchubertContext& p);
  ~ClosureIterator() {};
/* iterator operators */
  operator bool() const { return not state.empty(); }
  void operator++();
// accessors (only valid if not past the end)
  coxtypes::CoxNbr current() const { return state.top().current; }
  const bitmap::BitMap& closure() const { return state.top().closure; }
}; // |class ClosureIterator|

}; // |namespace schubert|

#endif

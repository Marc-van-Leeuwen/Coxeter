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
#include "bits.h"
#include "interface.h"

/******** type declarations *************************************************/

namespace schubert {
  class ClosureIterator;
  class SchubertContext;
  class StandardSchubertContext;

  struct NFCompare;

  using CoxNbrList = containers::vector<coxtypes::CoxNbr>;
  using Homology   = containers::vector<coxtypes::BettiNbr>;
};

/******** function declarations *********************************************/


namespace schubert {
  Homology betti(coxtypes::CoxNbr y, const SchubertContext& p);
  bool is_involution(const SchubertContext& p,coxtypes::CoxNbr x);
  containers::sl_list<Ulong> indices_of_maxima
  (const SchubertContext& p, containers::vector<coxtypes::CoxNbr>& c);
  void select_maxima_for
    (const SchubertContext& p, Lflags f, bitmap::BitMap& b);
  coxtypes::Generator first_flagged(GenSet f, const bits::Permutation& order);
  bool shortlex_leq(const SchubertContext& p, const bits::Permutation& order,
		    coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  Ulong sum(const Homology& h);
  bitmap::BitMap closure(const SchubertContext& p,coxtypes::CoxNbr x);

  void print(FILE* file, SchubertContext& p);
  void printPartition
    (FILE* file, const bits::Partition& pi,
     const SchubertContext& p, const interface::Interface& I);

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
  { return shortlex_leq(p,order,x,y); }
};


class SchubertContext
{
  const graph::CoxGraph& d_graph;
  containers::vector<coxtypes::Length> d_length;
  containers::vector<CoxNbrList> d_hasse;
  containers::vector<Lflags> d_descent;
  containers::matrix<coxtypes::CoxNbr> d_shift;
  containers::matrix<coxtypes::CoxNbr> d_star; // indexed by |CoxNbr,edge_nr|
  containers::vector<bitmap::BitMap> d_downset; // length |2*d_rank|
  bitmap::BitMap d_parity[2]; // array of TWO (complementary) parity bitmaps
  coxtypes::Length d_maxlength; // maximal length

public:
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(SchubertContext));}
/* constructors and destructors */
  SchubertContext(const graph::CoxGraph& G);

  // inlined accessor methods
  coxtypes::Rank rank() const { return d_graph.rank(); }
  coxtypes::CoxNbr size() const { return d_length.size(); }
  coxtypes::Length length(coxtypes::CoxNbr x) const { return d_length[x]; }
  const bitmap::BitMap& parity(coxtypes::CoxNbr x) const
    { return d_parity[d_length[x]%2]; }
  coxtypes::Length max_length() const { return d_maxlength; }

  bool in_context(coxtypes::CoxNbr x) const { return x<size(); }

  Lflags descent(coxtypes::CoxNbr x) const
    { return d_descent[x]; }
  GenSet rdescent(coxtypes::CoxNbr x) const
    { return d_descent[x] & constants::lt_mask[rank()]; }
  GenSet ldescent(coxtypes::CoxNbr x) const
    { return d_descent[x] >> rank(); } // left descents as (neutral) generators
  template<char side> Lflags descent_set(coxtypes::CoxNbr x) const
    { return side=='l' ? ldescent(x) : side=='r' ? rdescent(x) : descent(x); }

  Lflags ascent(coxtypes::CoxNbr x) const
    { return ~d_descent[x]&constants::lt_mask[2*rank()]; }
  GenSet rascent(coxtypes::CoxNbr x) const
    {return ~rdescent(x)&constants::lt_mask[rank()];}
  GenSet lascent(coxtypes::CoxNbr x) const
    { return ~ldescent(x)&constants::lt_mask[rank()]; }
  template<char side> Lflags ascent_set(coxtypes::CoxNbr x) const
    { return side=='l' ? lascent(x) : side=='r' ? rascent(x) : ascent(x); }

  template<char side>
    coxtypes::Generator first_descent(coxtypes::CoxNbr x) const
    { return constants::first_bit(descent_set<side>(x)); }
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
    { return first_flagged(ldescent(x),order); }
  coxtypes::Generator firstRDescent
    (coxtypes::CoxNbr x, const bits::Permutation& order) const
    { return first_flagged(rdescent(x),order); }

  template<char side>
    coxtypes::Generator is_descent
    (coxtypes::CoxNbr x, coxtypes::Generator s) const
    { return (descent_set<side>(x)&constants::eq_mask[s])!=0; }
  bool isDescent(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { return (d_descent[x]&constants::eq_mask[s])!=0; } // whether |s| is descent
  bool is_ascent(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { return (d_descent[x]&constants::eq_mask[s])==0; } // whether |s| is ascent

  coxtypes::CoxNbr shift(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { return d_shift.entry(x,s); } // left or right shift
  coxtypes::CoxNbr rshift(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { assert(s<rank()); return d_shift.entry(x,s); }
  coxtypes::CoxNbr lshift(coxtypes::CoxNbr x, coxtypes::Generator s) const
    { assert(s<rank()); return d_shift.entry(x,rank()+s); }

  const bitmap::BitMap& down_set(coxtypes::Generator s) const
    { return d_downset[s]; }
  const CoxNbrList& hasse(coxtypes::CoxNbr x) const { return d_hasse[x]; }
  const type::Type& type() const { return d_graph.type(); }

  Lflags S() const // mask for simple generators, to then restrict to |GenSet|
    { return constants::lt_mask[rank()]; }

/* accessors */
  bitmap::BitMap closure(coxtypes::CoxNbr x) const;
  void spread_subset
    (bitmap::BitMap& q, CoxNbrList& elements, coxtypes::Generator s) const;

  coxtypes::Cox_word word(coxtypes::CoxNbr x) const;
  coxtypes::CoxNbr context_number(const coxtypes::CoxWord& g) const;
  coxtypes::CoxNbr context_number(const coxtypes::Cox_word& g) const;
  void extendSubSet(bits::SubSet& q, coxtypes::Generator s) const;
  Lflags twoDescent(coxtypes::CoxNbr x) const;
  bool Bruhat_leq(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const;
  bool isSuperExtremal(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const;
  coxtypes::CoxNbr maximize(coxtypes::CoxNbr x, const Lflags& f) const;
  coxtypes::CoxNbr minimize(coxtypes::CoxNbr x, const Lflags& f) const;
  coxtypes::CoxWord& normalForm
    (coxtypes::CoxWord& g, coxtypes::CoxNbr x, const bits::Permutation& order)
    const;
  Ulong nStarOps() const { return d_graph.finite_edges().size(); }

  // manipulators
  coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g);
  coxtypes::CoxNbr extend_context(const coxtypes::Cox_word& g);
  coxtypes::CoxNbr star(coxtypes::CoxNbr x, Ulong r)
  { assert(x<size());
    if (d_star.nr_rows()<size())
      fill_star_table();
    return d_star.entry(x,r);
  }
  template<char side> // either 'l' or 'r'
    const coxtypes::CoxNbr* star_base(coxtypes::CoxNbr x)
  { assert(x<size());
    if (d_star.nr_rows()<size())
      fill_star_table();
    return d_star.row(x) + (side=='l' ? nStarOps() : 0);
  }

  void revertSize(Ulong n);
  void permute(const bits::Permutation& a);
/* i/o */
  std::string& append(std::string&, coxtypes::CoxNbr x) const;
  std::string& append(std::string&, coxtypes::CoxNbr x,
		      const interface::Interface& I) const;
  void print(FILE* file, coxtypes::CoxNbr x) const;
  void print(FILE* file, coxtypes::CoxNbr x, const interface::Interface& I)
    const;

private:
  void extend_context
    (bitmap::BitMap& q, CoxNbrList& elements, coxtypes::Generator s);
  void increase_size(Ulong n);
  void fill_Hasse(Ulong first, coxtypes::Generator s);
  void fill_shifts_and_descents(coxtypes::CoxNbr first, coxtypes::Generator s);

  void fill_star_table();
  coxtypes::CoxNbr* right_star_row(coxtypes::CoxNbr x)
    { return d_star.row(x); }
  coxtypes::CoxNbr* left_star_row(coxtypes::CoxNbr x)
    { return d_star.row(x)+nStarOps(); }
  template<bool left> // get decreasing "star" image, or undefined
    coxtypes::CoxNbr falling_star(coxtypes::CoxNbr x, GenSet st) const;

}; // |class SchubertContext|


class ClosureIterator {
 private:
  struct node
  { bitmap::BitMap closure; // the Bruhat interval $[e,current]$
    coxtypes::CoxNbr closure_size; // the size of that interval (for speed)
    coxtypes::CoxNbr current; // the current Coxeter element
    GenSet asc; // its not yet explored ascent set
    node(Ulong capacity) : closure(capacity),closure_size(),current(),asc() {}
  };
  const SchubertContext& d_schubert;
  containers::vector<node> state;
  Ulong sp; // state pointer (equals index of top |node|)
  CoxNbrList elements; // closure, for fast traversal
  bitmap::BitMap visited; // everything seen so far
 public:
/* constructors and destructors */
  ClosureIterator(const SchubertContext& p);
  ~ClosureIterator() {};
/* iterator operators */
  operator bool() const { return not state.empty(); }
  void operator++();
// accessors (only valid if not past the end)
  coxtypes::CoxNbr current() const { return state[sp].current; }
  const bitmap::BitMap& closure() const { return state[sp].closure; }
}; // |class ClosureIterator|

}; // |namespace schubert|

#endif

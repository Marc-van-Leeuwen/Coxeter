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

  typedef list::List<coxtypes::CoxNbr> CoatomList;
  typedef list::List<coxtypes::BettiNbr> Homology;
};

/******** function declarations *********************************************/


namespace schubert {
  void betti(Homology& h, const coxtypes::CoxNbr& y, const SchubertContext& p);
  void extractInvolutions(const SchubertContext& p, bits::BitMap& b);
  void extractMaximals(const SchubertContext& p, list::List<coxtypes::CoxNbr>& c);
  void extractMaximals(const SchubertContext& p, list::List<coxtypes::CoxNbr>& c,
		       list::List<Ulong>& a);
  void select_maxima_for
    (const SchubertContext& p, bits::BitMap& b, const bits::Lflags& f);
  Ulong min(const bits::Set& c, NFCompare& nfc);
  Ulong minDescent(const bits::Lflags& f, const bits::Permutation& order);
  void minimize(const SchubertContext& p, bits::BitMap& b, const bits::Lflags& f);
  void print(FILE* file, const SchubertContext& p);
  void printBitMap(FILE* file, const bits::BitMap& pi, const SchubertContext& p,
		   const interface::Interface& I);
  void printList(FILE* file, const list::List<coxtypes::CoxNbr>& v, const SchubertContext& p,
		 const interface::Interface& I);
  void printPartition(FILE* file, const bits::Partition& pi,
		      const SchubertContext& p, const interface::Interface& I);
  void printPartition(FILE* file, const bits::Partition& pi, const bits::BitMap& b,
		      const SchubertContext& p, const interface::Interface& I);
  void readBitMap(list::List<coxtypes::CoxNbr>& c, const bits::BitMap& b);
  void read_bitmap(containers::vector<coxtypes::CoxNbr>& c, const bits::BitMap& b);
  bool shortLexOrder(const SchubertContext& p, const coxtypes::CoxNbr& x,
		     const coxtypes::CoxNbr& y, const bits::Permutation& order);
  void setBitMap(bits::BitMap& b, const list::List<coxtypes::CoxNbr>& c);
  Ulong sum(const Homology& h);
};

/******** type definitions *************************************************/

namespace schubert {

struct NFCompare {
  const SchubertContext& p;
  const bits::Permutation& order;
  NFCompare(const SchubertContext& q, const bits::Permutation& generator_ordering)
    :p(q),order(generator_ordering) {};
  ~NFCompare() {};
  bool operator()(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)
  {return shortLexOrder(p,x,y,order);}
};

class SchubertContext {
  friend class ClosureIterator;
 public:
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(SchubertContext));}
  virtual ~SchubertContext() {};
/* accessors */
  virtual coxtypes::CoxWord& append(coxtypes::CoxWord& g, const coxtypes::CoxNbr& x) const = 0;
  virtual bits::Lflags ascent(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::CoxNbr contextNumber(const coxtypes::CoxWord& g) const = 0;
  virtual bits::Lflags descent(const coxtypes::CoxNbr& x) const = 0;
  virtual const bits::BitMap& downset(const coxtypes::Generator& s) const = 0;
  virtual void extendSubSet(bits::SubSet& q, const coxtypes::Generator& s) const = 0;
  virtual void extractClosure(bits::BitMap& b, const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::Generator firstDescent(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::Generator firstLDescent(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::Generator firstRDescent(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::Generator firstDescent(const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const = 0;
  virtual coxtypes::Generator firstLDescent(const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const = 0;
  virtual coxtypes::Generator firstRDescent(const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const = 0;
  virtual const CoatomList& hasse(const coxtypes::CoxNbr& x) const = 0;
  virtual bool inOrder(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const = 0;
  virtual bool isDescent(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const = 0;
  virtual bits::Lflags lascent(const coxtypes::CoxNbr& x) const = 0;
  virtual bits::Lflags ldescent(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::Length length(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::CoxNbr lshift(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const = 0;
  virtual coxtypes::CoxNbr maximize(const coxtypes::CoxNbr& x, const bits::Lflags& f) const = 0;
  virtual coxtypes::Length maxlength() const = 0;
  virtual coxtypes::CoxNbr minimize(const coxtypes::CoxNbr& x, const bits::Lflags& f) const = 0;
  virtual coxtypes::CoxWord& normalForm(coxtypes::CoxWord& g, const coxtypes::CoxNbr& x,
			      const bits::Permutation& order) const = 0;
  virtual Ulong nStarOps() const = 0;
  virtual const bits::BitMap& parity(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::Rank rank() const = 0;
  virtual bits::Lflags rascent(const coxtypes::CoxNbr& x) const = 0;
  virtual bits::Lflags rdescent(const coxtypes::CoxNbr& x) const = 0;
  virtual coxtypes::CoxNbr rshift(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const = 0;
  virtual bits::Lflags S() const = 0;
  virtual coxtypes::CoxNbr shift(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const = 0;
  virtual coxtypes::CoxNbr size() const = 0;
  virtual coxtypes::CoxNbr star(coxtypes::CoxNbr x, const Ulong& r) const = 0;
  virtual bits::Lflags twoDescent(const coxtypes::CoxNbr& x) const = 0;
  virtual const type::Type& type() const = 0;
/* modifiers */
  virtual coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g) = 0;
  virtual void permute(const bits::Permutation& a) = 0;
  virtual void revertSize(const Ulong& n) = 0;
  virtual void setSize(const Ulong& n) = 0;
/* input-output */
  virtual std::string& append(std::string&, const coxtypes::CoxNbr& x) const = 0;
  virtual std::string& append(std::string&, const coxtypes::CoxNbr& x, const interface::Interface& I)
    const = 0;
  virtual void print(FILE* file, const coxtypes::CoxNbr& x) const = 0;
  virtual void print(FILE* file, const coxtypes::CoxNbr& x, const interface::Interface& I)
    const = 0;
};

class StandardSchubertContext:public SchubertContext {
 private:
/* private class declaration */
  class ContextExtension {
  private:
    StandardSchubertContext& d_schubert;
    Ulong d_size;
    coxtypes::CoxNbr* d_shift;
    coxtypes::CoxNbr* d_star;
  public:
    void* operator new(size_t size) {return memory::arena().alloc(size);}
    void operator delete(void* ptr)
      {return memory::arena().free(ptr,sizeof(ContextExtension));}
    ContextExtension(StandardSchubertContext& p, const Ulong& c);
    ~ContextExtension();
    Ulong size() {return d_size;}
  };
  const graph::CoxGraph& d_graph;
  coxtypes::Rank d_rank;
  coxtypes::Length d_maxlength;
  coxtypes::CoxNbr d_size;
  list::List<coxtypes::Length> d_length;
  list::List<CoatomList> d_hasse;
  list::List<bits::Lflags> d_descent;
  list::List<coxtypes::CoxNbr*> d_shift;
  list::List<coxtypes::CoxNbr*> d_star;
  bits::BitMap* d_downset;
  bits::BitMap* d_parity;
  bits::SubSet d_subset;
  stack::Stack<ContextExtension*> d_history;
/* private member functions */
  void fillCoatoms(const Ulong& first, const coxtypes::Generator& s);
  void fillDihedralShifts(const coxtypes::CoxNbr& x, const coxtypes::Generator& s);
  void fillShifts(const coxtypes::CoxNbr& first, const coxtypes::Generator& s);
  void fillStar(const coxtypes::CoxNbr& first);
  void fullExtension(bits::SubSet& q, const coxtypes::Generator& s);
  void subSetExtension(bits::SubSet& q, const coxtypes::Generator& s) const;
 public:
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(StandardSchubertContext));}
/* friend declaration */
  friend ContextExtension::ContextExtension(StandardSchubertContext&,
					    const Ulong& c);
  friend ContextExtension::~ContextExtension();
/* constructors and destructors */
  StandardSchubertContext(const graph::CoxGraph& G);
  ~StandardSchubertContext();
/* accessors */
  coxtypes::CoxWord& append(coxtypes::CoxWord& g, const coxtypes::CoxNbr& x) const;
  bits::Lflags ascent(const coxtypes::CoxNbr& x) const;                          /* inlined */
  coxtypes::CoxNbr contextNumber(const coxtypes::CoxWord& g) const;
  bits::Lflags descent(const coxtypes::CoxNbr& x) const;                         /* inlined */
  const bits::BitMap& downset(const coxtypes::Generator& s) const;               /* inlined */
  void extendSubSet(bits::SubSet& q, const coxtypes::Generator& s) const;
  void extractClosure(bits::BitMap& b, const coxtypes::CoxNbr& x) const;
  coxtypes::Generator firstDescent(const coxtypes::CoxNbr& x) const;                 /* inlined */
  coxtypes::Generator firstLDescent(const coxtypes::CoxNbr& x) const;                /* inlined */
  coxtypes::Generator firstRDescent(const coxtypes::CoxNbr& x) const;                /* inlined */
  coxtypes::Generator firstDescent(const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const;                                                       /* inlined */
  coxtypes::Generator firstLDescent(const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const;                                                       /* inlined */
  coxtypes::Generator firstRDescent(const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const;                                                       /* inlined */
  const CoatomList& hasse(const coxtypes::CoxNbr& x) const;                /* inlined */
  bool inOrder(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const;
  bool isDescent(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const;     /* inlined */
  bool isSuperExtremal(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y) const;
  bits::Lflags lascent(const coxtypes::CoxNbr& x) const;                         /* inlined */
  bits::Lflags ldescent(const coxtypes::CoxNbr& x) const;                        /* inlined */
  coxtypes::Length length(const coxtypes::CoxNbr& x) const;                          /* inlined */
  coxtypes::CoxNbr lshift(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const;      /* inlined */
  coxtypes::CoxNbr maximize(const coxtypes::CoxNbr& x, const bits::Lflags& f) const;
  coxtypes::Length maxlength() const;                                      /* inlined */
  coxtypes::CoxNbr minimize(const coxtypes::CoxNbr& x, const bits::Lflags& f) const;
  coxtypes::CoxWord& normalForm(coxtypes::CoxWord& g, const coxtypes::CoxNbr& x, const bits::Permutation& order)
    const;
  Ulong nStarOps() const;                                      /* inlined */
  const bits::BitMap& parity(const coxtypes::CoxNbr& x) const;                   /* inlined */
  coxtypes::Rank rank() const;                                             /* inlined */
  bits::Lflags rascent(const coxtypes::CoxNbr& x) const;                         /* inlined */
  bits::Lflags rdescent(const coxtypes::CoxNbr& x) const;                        /* inlined */
  coxtypes::CoxNbr rshift(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const;      /* inlined */
  bits::Lflags S() const;                                              /* inlined */
  coxtypes::CoxNbr shift(const coxtypes::CoxNbr& x, const coxtypes::Generator& s) const;       /* inlined */
  coxtypes::CoxNbr size() const;                                           /* inlined */
  coxtypes::CoxNbr star(coxtypes::CoxNbr x, const Ulong& r) const;                 /* inlined */
  bits::Lflags twoDescent(const coxtypes::CoxNbr& x) const;
  const type::Type& type() const;                                      /* inlined */
/* manipulators */
  coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);
  void setSize(const Ulong& n);
/* i/o */
  std::string& append(std::string&, const coxtypes::CoxNbr& x) const;
  std::string& append(std::string&, const coxtypes::CoxNbr& x,
		      const interface::Interface& I) const;
  void print(FILE* file, const coxtypes::CoxNbr& x) const;
  void print(FILE* file, const coxtypes::CoxNbr& x, const interface::Interface& I) const;
};

class ClosureIterator {
 private:
  const SchubertContext& d_schubert;
  bits::SubSet d_subSet;
  coxtypes::CoxWord d_g;
  list::List<Ulong> d_subSize;
  bits::BitMap d_visited;
  coxtypes::CoxNbr d_current;
  bool d_valid;
/* private functions */
  void update(const coxtypes::CoxNbr& x, const coxtypes::Generator& s);
 public:
/* constructors and destructors */
  ClosureIterator(const SchubertContext& p);
  ~ClosureIterator() {};
/* iterator operators */
  operator bool() const;                                         /* inlined */
  void operator++();
  const bits::SubSet& operator()() const;                              /* inlined */
/* accessors */
  const coxtypes::CoxNbr& current() const;                                 /* inlined */
};

};

/******** inline definitions **********************************************/

namespace schubert {
  inline bits::Lflags StandardSchubertContext::ascent(const coxtypes::CoxNbr& x) const
    {return ~d_descent[x]&constants::leqmask[2*d_rank-1];}
  inline bits::Lflags StandardSchubertContext::descent(const coxtypes::CoxNbr& x) const
    {return d_descent[x];}
  inline const bits::BitMap& StandardSchubertContext::downset(const coxtypes::Generator& s)
    const {return d_downset[s];}
  inline coxtypes::Generator StandardSchubertContext::firstDescent(const coxtypes::CoxNbr& x) const
    {return constants::firstBit(descent(x));}
  inline coxtypes::Generator StandardSchubertContext::firstLDescent(const coxtypes::CoxNbr& x)
    const {return constants::firstBit(ldescent(x));}
  inline coxtypes::Generator StandardSchubertContext::firstRDescent(const coxtypes::CoxNbr& x)
    const {return constants::firstBit(rdescent(x));}
  inline coxtypes::Generator StandardSchubertContext::firstDescent(const coxtypes::CoxNbr& x,
    const bits::Permutation& order) const {return firstRDescent(x,order);}
  inline coxtypes::Generator StandardSchubertContext::firstLDescent(const coxtypes::CoxNbr& x,
    const bits::Permutation& order) const {return minDescent(ldescent(x),order);}
  inline coxtypes::Generator StandardSchubertContext::firstRDescent(const coxtypes::CoxNbr& x,
    const bits::Permutation& order) const {return minDescent(rdescent(x),order);}
  inline const CoatomList& StandardSchubertContext::hasse(const coxtypes::CoxNbr& x)
    const {return d_hasse[x];}
  inline bool StandardSchubertContext::isDescent(const coxtypes::CoxNbr& x,
						 const coxtypes::Generator& s)
    const {return d_descent[x]&constants::lmask[s];}
  inline bits::Lflags StandardSchubertContext::lascent(const coxtypes::CoxNbr& x) const
    {return ~ldescent(x)&constants::leqmask[d_rank-1];}
  inline bits::Lflags StandardSchubertContext::ldescent(const coxtypes::CoxNbr& x) const
    {return d_descent[x] >> d_rank;}
  inline coxtypes::Length StandardSchubertContext::length(const coxtypes::CoxNbr& x) const
    {return d_length[x];}
  inline coxtypes::CoxNbr StandardSchubertContext::lshift(const coxtypes::CoxNbr& x,
					    const coxtypes::Generator& s)
    const {return d_shift[x][d_rank+s];}
  inline coxtypes::Length StandardSchubertContext::maxlength() const
    {return d_maxlength;}
  inline Ulong StandardSchubertContext::nStarOps() const
    {return d_graph.starOps().size();}
  inline const bits::BitMap& StandardSchubertContext::parity(const coxtypes::CoxNbr& x) const
    {return d_parity[d_length[x]%2];}
  inline coxtypes::Rank StandardSchubertContext::rank() const {return d_rank;}
  inline bits::Lflags StandardSchubertContext::rascent(const coxtypes::CoxNbr& x) const
    {return ~rdescent(x)&constants::leqmask[d_rank-1];}
  inline bits::Lflags StandardSchubertContext::rdescent(const coxtypes::CoxNbr& x) const
    {return d_descent[x] & constants::leqmask[d_rank-1];}
  inline coxtypes::CoxNbr StandardSchubertContext::rshift
    (const coxtypes::CoxNbr& x, const coxtypes::Generator& s)
    const {return d_shift[x][s];}
  inline bits::Lflags StandardSchubertContext::S() const
    {return constants::leqmask[d_rank-1];}
  inline coxtypes::CoxNbr StandardSchubertContext::shift
    (const coxtypes::CoxNbr& x, const coxtypes::Generator& s)
    const {return d_shift[x][s];}
  inline coxtypes::CoxNbr StandardSchubertContext::size() const {return d_size;}
  inline coxtypes::CoxNbr StandardSchubertContext::star(coxtypes::CoxNbr x,
					      const Ulong& r) const
    {return d_star[x][r];}
  inline const type::Type& StandardSchubertContext::type() const
    {return d_graph.type();}

  inline const bits::SubSet& ClosureIterator::operator()() const {return d_subSet;}
  inline ClosureIterator::operator bool() const {return d_valid;}
  inline const coxtypes::CoxNbr& ClosureIterator::current() const
    {return d_current;}
};

#endif

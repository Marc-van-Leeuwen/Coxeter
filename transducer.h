/*
  This is transducer.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef TRANSDUCER_H  /* guarantee single inclusion */
#define TRANSDUCER_H

#include "globals.h"

/******** type declarations *************************************************/

namespace transducer {
  class FiltrationTerm;
  class SubQuotient;
  class Transducer;
};

/******** constants *********************************************************/

#include "coxtypes.h"

namespace transducer {
  static const coxtypes::ParNbr undef_parnbr = coxtypes::PARNBR_MAX + 1;
};

/******** type definitions **************************************************/

#include "graph.h"
#include "list.h"
#include "memory.h"

class transducer::SubQuotient {
 private:
/* typedef in class scope */
  typedef list::List<coxtypes::ParNbr> Sub_set;
/* data  */
  coxtypes::Rank d_rank;
  Ulong d_size;
  graph::CoxGraph& d_graph;
  list::List<coxtypes::ParNbr> d_shift;
  list::List<coxtypes::Length> d_length;
  coxtypes::ParNbr& shiftref(coxtypes::ParNbr x, coxtypes::Generator s)
    {return d_shift[x*d_rank+s];}
  coxtypes::Length& lengthref(coxtypes::ParNbr x) {return d_length[x];}
 public:
/* constructors  */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(SubQuotient));}
  SubQuotient(graph::CoxGraph& G, coxtypes::Rank l);
  ~SubQuotient();
/* manipulators */
  coxtypes::ParNbr extend(coxtypes::ParNbr x, coxtypes::Generator s);
  void fill(const graph::CoxGraph& G);
/* accessors */
  coxtypes::Generator firstDescent(const coxtypes::ParNbr& x) const;
  coxtypes::Length length(const coxtypes::ParNbr& x) const;   /* inlined */
  coxtypes::Rank rank() const;                                /* inlined */
  coxtypes::CoxWord& reduced(coxtypes::CoxWord& g, coxtypes::ParNbr x) const;
  coxtypes::ParNbr shift
    (const coxtypes::ParNbr& x, const coxtypes::Generator& s) const;
                                                              /* inlined */
  Ulong size() const;                                         /* inlined */
  void schubertClosure(Sub_set& Q, coxtypes::ParNbr x);
};

class transducer::FiltrationTerm {
  SubQuotient *d_X;
  FiltrationTerm *d_next;
  list::List<coxtypes::CoxWord> d_np;
  void fillNormalPieces();
 public:
/* constructors  */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(FiltrationTerm));}
  FiltrationTerm() {};
  FiltrationTerm(graph::CoxGraph& G, coxtypes::Rank l, FiltrationTerm* p = 0);
  ~FiltrationTerm();
/* manipulators */
  coxtypes::ParNbr extend
    (const coxtypes::ParNbr& x, const coxtypes::Generator& s);    /* inlined */
  void fill(const graph::CoxGraph& G);                            /* inlined */
/* accessors */
  coxtypes::CoxLetter firstDescent
    (const coxtypes::ParNbr& x) const;                            /* inlined */
  coxtypes::Length length(const coxtypes::ParNbr& x) const;       /* inlined */
  FiltrationTerm *next() const;                                   /* inlined */
  const coxtypes::CoxWord& np(const coxtypes::ParNbr& x) const;   /* inlined */
  coxtypes::Rank rank() const;                                    /* inlined */
  coxtypes::CoxWord& reduced
    (coxtypes::CoxWord& g, coxtypes::ParNbr x) const;             /* inlined */
  coxtypes::ParNbr shift
    (const coxtypes::ParNbr& x, const coxtypes::Generator& s) const;/* inlined */
  Ulong size() const;                                             /* inlined */
};

class transducer::Transducer {
 private:
  list::List<FiltrationTerm> d_filtration;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(Transducer));}
  Transducer(graph::CoxGraph& G);
  ~Transducer();
/* manipulators */
  FiltrationTerm* transducer(const coxtypes::Rank& l) {return d_filtration.ptr()+l;}
/* accessors */
  const FiltrationTerm* transducer(const coxtypes::Rank& l) const;
};

/******** inline implementations ********************************************/

namespace transducer {

/* SubQuotient */

inline coxtypes::Length SubQuotient::length(const coxtypes::ParNbr& x) const {return d_length[x];}
inline coxtypes::Rank SubQuotient::rank() const {return d_rank;}
inline coxtypes::ParNbr SubQuotient::shift(const coxtypes::ParNbr& x, const coxtypes::Generator& s) const
  {return d_shift[x*d_rank+s];}
inline Ulong SubQuotient::size() const {return d_size;}

/* FiltrationTerm */

inline coxtypes::ParNbr FiltrationTerm::extend
  (const coxtypes::ParNbr& x, const coxtypes::Generator& s)
  {return d_X->extend(x,s);}
inline void FiltrationTerm::fill(const graph::CoxGraph& G)
  {d_X->fill(G); fillNormalPieces();}
inline coxtypes::CoxLetter FiltrationTerm::firstDescent
  (const coxtypes::ParNbr& x) const
  {return d_X->firstDescent(x);}
inline coxtypes::Length FiltrationTerm::length(const coxtypes::ParNbr& x) const
  {return d_X->length(x);}
inline FiltrationTerm* FiltrationTerm::next() const {return d_next;}
inline const coxtypes::CoxWord& FiltrationTerm::np
  (const coxtypes::ParNbr& x) const
  {return d_np[x];}
inline coxtypes::Rank FiltrationTerm::rank() const {return d_X->rank();}
inline coxtypes::CoxWord& FiltrationTerm::reduced(coxtypes::CoxWord& g, coxtypes::ParNbr x) const
  {return d_X->reduced(g,x);}
inline coxtypes::ParNbr FiltrationTerm::shift
  (const coxtypes::ParNbr& x, const coxtypes::Generator& s) const
  {return d_X->shift(x,s);}
inline Ulong FiltrationTerm::size() const {return d_X->size();}

/* Transducer */

inline const FiltrationTerm* Transducer::transducer
  (const coxtypes::Rank& l) const
  {return d_filtration.ptr()+l;}

};

#endif

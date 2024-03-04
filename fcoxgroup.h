/*
  This is fcoxgroup.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/****************************************************************************

  This module defines the hierarchy of finite Coxeter groups; see coxgroup.h
  for the layout of the general hierarchy of Coxeter groups. As explained
  there, many of these classes are barely implemented. We also explained
  that only leaf classes in the hierarchy are concrete; in preparation
  for possible derivation of Finite*RankCoxGroup, we have shadowed them
  by the concrete classes F*RCoxGroup; the same holds for SmallCoxGroup.

 ****************************************************************************/

#ifndef FCOXGROUP_H  /* guard against multiple inclusions */
#define FCOXGROUP_H

#include "globals.h"
#include "coxgroup.h"

/******** type declarations *************************************************/

namespace fcoxgroup {
  class FiniteCoxGroup;
  class FiniteBigRankCoxGroup;
  class GeneralFBRCoxGroup;
  class FiniteMedRankCoxGroup;
  class GeneralFMRCoxGroup;
  class FiniteSmallRankCoxGroup;
  class GeneralFSRCoxGroup;
  class SmallCoxGroup;
  class GeneralSCoxGroup;
  typedef coxtypes::CoxNbr DenseArray;
};

/******** function declarations *********************************************/

namespace fcoxgroup {
  bool isFiniteType(coxgroup::CoxGroup *W);
  coxtypes::Rank maxSmallRank(const type::Type& x);
};

/******** type definitions **************************************************/

namespace fcoxgroup {

class FiniteCoxGroup : public coxgroup::CoxGroup
{
 protected:

  coxtypes::CoxArr d_longest_coxarr;
  coxtypes::CoxWord d_longest_coxword;
  coxtypes::Length d_maxlength;
  coxtypes::CoxSize d_order;
  transducer::Transducer *d_transducer;
  bits::Partition d_lcell;
  bits::Partition d_rcell;
  bits::Partition d_lrcell;
  bits::Partition d_luneqcell;
  bits::Partition d_runeqcell;
  bits::Partition d_lruneqcell;
  bits::Partition d_ldescent;
  bits::Partition d_rdescent;
  bits::Partition d_ltau;
  bits::Partition d_rtau;
  bits::Partition d_lstring;
  bits::Partition d_rstring;
  list::List<coxtypes::CoxNbr> d_duflo;

 public:

/* constructors */

  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(FiniteCoxGroup));}

  FiniteCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~FiniteCoxGroup();

/* accessors */

  bool isFullContext() const;
  const coxtypes::CoxArr& longest_coxarr() const;                /* inlined */
  const coxtypes::CoxWord& longest_coxword() const;              /* inlined */
  coxtypes::Length maxLength() const;                            /* inlined */
  void modify(interface::ParseInterface& P, const interface::Token& tok) const;
  coxtypes::CoxSize order() const;                               /* inlined */
  bool parseModifier(interface::ParseInterface& P) const;
  const transducer::FiltrationTerm *transducer(const coxtypes::Rank& l = 0) const;     /* inlined */

/* modifiers */

  transducer::FiltrationTerm *transducer(const coxtypes::Rank& l = 0)
    {return d_transducer->transducer(l);}

/* array operations */

  const coxtypes::CoxArr& assign(coxtypes::CoxArr& a, const coxtypes::CoxArr& b) const;        /* inlined */
  virtual const coxtypes::CoxArr& inverseArr(coxtypes::CoxArr& a) const;
  coxtypes::Length length(const coxtypes::CoxArr& a) const;
  const coxtypes::CoxArr& powerArr(coxtypes::CoxArr& a, const Ulong& m) const;
  int prodArr(coxtypes::CoxArr& a, const coxtypes::CoxArr& b) const;
  GenSet rDescent(const coxtypes::CoxArr& a) const;
  const coxtypes::CoxWord& reducedArr(coxtypes::CoxWord& g, const coxtypes::CoxArr& a) const;
  const coxtypes::CoxArr& setZero(coxtypes::CoxArr& a) const;                        /* inlined */

/* mixed operations */

  const coxtypes::CoxArr& assign(coxtypes::CoxArr& a, const coxtypes::CoxWord& g) const;
  int prodArr(coxtypes::CoxArr& a, coxtypes::Generator s) const;
  int prodArr(coxtypes::CoxArr& a, const coxtypes::CoxWord& g) const;

// manipulators

  void fullContext();                                            /* inlined */

/* kazhdan-lusztig cells and realted partitions */

  const bits::Partition& lCell();
  const bits::Partition& lrCell();
  const bits::Partition& rCell();
  const list::List<coxtypes::CoxNbr>& duflo();
  const bits::Partition& lUneqCell();
  const bits::Partition& lrUneqCell();
  const bits::Partition& rUneqCell();
  const bits::Partition& lDescent();
  const bits::Partition& rDescent();
  const bits::Partition& lString();
  const bits::Partition& rString();
  const bits::Partition& lTau();
  const bits::Partition& rTau();

};

class FiniteBigRankCoxGroup : public FiniteCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(FiniteBigRankCoxGroup));}
  FiniteBigRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~FiniteBigRankCoxGroup();
/* accessors */
  kl::KLContext& kl() const;                                     /* inlined */
};

 class GeneralFBRCoxGroup : public FiniteBigRankCoxGroup { /* leaf class */
  public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralFBRCoxGroup));}
  GeneralFBRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralFBRCoxGroup();
};

class FiniteMedRankCoxGroup : public FiniteCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(FiniteMedRankCoxGroup));}
  FiniteMedRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~FiniteMedRankCoxGroup();
/* accessors */
  kl::KLContext& kl() const;                                     /* inlined */
};

 class GeneralFMRCoxGroup : public FiniteMedRankCoxGroup { /* leaf class */
  public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralFMRCoxGroup));}
  GeneralFMRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralFMRCoxGroup();
};

class FiniteSmallRankCoxGroup : public FiniteMedRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(FiniteSmallRankCoxGroup));}
  FiniteSmallRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~FiniteSmallRankCoxGroup();
/* accessors */
  kl::KLContext& kl() const;                                     /* inlined */
};

 class GeneralFSRCoxGroup : public FiniteSmallRankCoxGroup { /* leaf class */
  public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralFSRCoxGroup));}
  GeneralFSRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralFSRCoxGroup();
};

class SmallCoxGroup : public FiniteSmallRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(SmallCoxGroup));}
  SmallCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~SmallCoxGroup();
/* accessors */
  const coxtypes::CoxArr& assign(coxtypes::CoxArr& a, const DenseArray& x) const;
  void assign(DenseArray& x, const coxtypes::CoxArr& a) const;
  bool parseDenseArray(interface::ParseInterface& P) const;
  virtual bool parseGroupElement(interface::ParseInterface& P) const;
  int prodD(coxtypes::CoxWord& g, const DenseArray& x) const;
  int prodD(DenseArray& x, const coxtypes::CoxWord& g) const;
};

 class GeneralSCoxGroup : public SmallCoxGroup { /* leaf class */
  public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralSCoxGroup));}
  GeneralSCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralSCoxGroup();
};

};

/******** Inline implementations ******************************************/

namespace fcoxgroup {

inline const coxtypes::CoxArr& FiniteCoxGroup::assign(coxtypes::CoxArr& a, const coxtypes::CoxArr& b) const
  {memmove(a,b,rank()*sizeof(coxtypes::ParNbr)); return a;}
inline void FiniteCoxGroup::fullContext()
  {extendContext(d_longest_coxword);}
inline const coxtypes::CoxArr& FiniteCoxGroup::longest_coxarr() const
  {return d_longest_coxarr;}
inline const coxtypes::CoxWord& FiniteCoxGroup::longest_coxword() const
  {return d_longest_coxword;}
inline coxtypes::Length FiniteCoxGroup::maxLength() const {return d_maxlength;}
inline coxtypes::CoxSize FiniteCoxGroup::order() const {return d_order;}
inline const coxtypes::CoxArr& FiniteCoxGroup::setZero(coxtypes::CoxArr& a) const
  {memset(a,0,rank()*sizeof(coxtypes::ParNbr)); return a;}
inline const transducer::FiltrationTerm* FiniteCoxGroup::transducer
  (const coxtypes::Rank& l) const
  {return d_transducer->transducer(l);}

};

#endif

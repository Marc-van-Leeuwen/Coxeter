/*
  This is invkl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef INVKL_H  /* guard against multiple inclusions */
#define INVKL_H

#include "globals.h"
#include "containers.h"
#include "coxtypes.h"
#include "klsupport.h"
#include "kl.h"
#include "hecke.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "memory.h"
#include "search.h"

/******** type declarations *************************************************/

namespace invkl {
  class KLContext;
  class MuFilter;

  typedef list::List<const kl::KLPol*> KLRow;
  typedef list::List<kl::MuData> MuRow;
  using HeckeElt = containers::vector<hecke::HeckeMonomial<kl::KLPol> >;
};

/******** type definitions **************************************************/


namespace invkl {

class KLContext {
  klsupport::KLSupport* d_klsupport;
  list::List<KLRow*> d_klList;
  list::List<MuRow*> d_muList;
  search::BinaryTree<kl::KLPol> d_klTree;
  kl::KLStats* d_status;
  struct KLHelper; /* provides helper functions */
  KLHelper* d_help;
  friend struct KLHelper;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLContext));}
  KLContext(klsupport::KLSupport* kls);
  ~KLContext();
/* accessors */
  const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y) const;               /* inlined */
  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const;                        /* inlined */
  const bitmap::BitMap& involution() const { return d_klsupport->involution(); }
  bool isExtrAllocated(const coxtypes::CoxNbr& x) const;                  /* inlined */
  bool isFullKL() const;                                        /* inlined */
  bool isFullMu() const;                                        /* inlined */
  bool isKLAllocated(const coxtypes::CoxNbr& x) const;                    /* inlined */
  bool isMuAllocated(const coxtypes::CoxNbr& x) const;                    /* inlined */
  const KLRow& klList(const coxtypes::CoxNbr& y) const;                   /* inlined */
  const klsupport::KLSupport& klsupport() const;                           /* inlined */
  coxtypes::Generator last(const coxtypes::CoxNbr& y) const;                        /* inlined */
  const MuRow& muList(const coxtypes::CoxNbr& y) const;                   /* inlined */
  coxtypes::Rank rank() const;                                            /* inlined */
  const schubert::SchubertContext& schubert() const
  { return d_klsupport->schubert(); }
  Ulong size() const;                                         /* inlined */
  const search::BinaryTree<kl::KLPol>& tree() const;                /* inlined */
/* manipulators */
  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a); /* inlined */
  void clearFullKL();                                            /* inlined */
  void clearFullMu();                                            /* inlined */
  void fillKL();
  void fillMu();
  const kl::KLPol& klPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
		     const coxtypes::Generator& s = coxtypes::undef_generator);
  klsupport::KLCoeff mu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
	     const coxtypes::Generator& s = coxtypes::undef_generator);
  void permute(const bits::Permutation& a);
  void row(HeckeElt& h, const coxtypes::CoxNbr& y);
  void revertSize(const Ulong& n);
  void setFullKL();                                              /* inlined */
  void setFullMu();                                              /* inlined */
  void setSize(const Ulong& n);
};

};

/******** inline definitions ************************************************/

namespace invkl {

inline const klsupport::ExtrRow& KLContext::extrList(const coxtypes::CoxNbr& y) const
  {return klsupport().extrList(y);}
inline coxtypes::CoxNbr KLContext::inverse(const coxtypes::CoxNbr& x) const
  {return d_klsupport->inverse(x);}
inline bool KLContext::isExtrAllocated(const coxtypes::CoxNbr& x) const
  {return d_klsupport->isExtrAllocated(x);}
inline bool KLContext::isFullKL() const
{return d_status->flags&kl::KLStats::kl_done;}
inline bool KLContext::isFullMu() const
{return d_status->flags&kl::KLStats::mu_done;}
inline bool KLContext::isKLAllocated(const coxtypes::CoxNbr& x) const
  {return d_klList[x] != 0;}
inline bool KLContext::isMuAllocated(const coxtypes::CoxNbr& x) const
  {return d_muList[x] != 0;}
inline const KLRow& KLContext::klList(const coxtypes::CoxNbr& y) const
  {return *d_klList[y];}
inline const klsupport::KLSupport& KLContext::klsupport() const {return *d_klsupport;}
inline coxtypes::Generator KLContext::last(const coxtypes::CoxNbr& y) const
  {return d_klsupport->last(y);}
inline const MuRow& KLContext::muList(const coxtypes::CoxNbr& y) const
  {return *d_muList[y];}
inline coxtypes::Rank KLContext::rank() const {return d_klsupport->rank();}
inline Ulong KLContext::size() const {return d_klList.size();}
  inline const search::BinaryTree<kl::KLPol>& KLContext::tree() const
  {return d_klTree;}

inline void KLContext::applyIPermutation(const coxtypes::CoxNbr& y,
					 const bits::Permutation& a)
  {return rightRangePermute(*d_klList[y],a);}
  inline void KLContext::clearFullKL() {d_status->flags &= ~kl::KLStats::kl_done;}
  inline void KLContext::clearFullMu() {d_status->flags &= ~kl::KLStats::mu_done;}
  inline void KLContext::setFullKL() {d_status->flags |= kl::KLStats::kl_done;}
  inline void KLContext::setFullMu() {d_status->flags |= kl::KLStats::mu_done;}

};

#endif

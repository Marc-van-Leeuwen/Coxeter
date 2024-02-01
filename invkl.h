/*
  This is invkl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef INVKL_H  /* guard against multiple inclusions */
#define INVKL_H

#include "globals.h"
#include "coxtypes.h"
#include "klsupport.h"
#include "hecke.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "memory.h"
#include "search.h"

/******** type declarations *************************************************/

namespace invkl {
  class KLContext;
  class KLPol;
  struct KLStatus;
  struct MuData;
  class MuFilter;

  typedef list::List<const KLPol*> KLRow;
  typedef list::List<MuData> MuRow;
  typedef list::List<hecke::HeckeMonomial<KLPol> > HeckeElt;
};

/******** function declarations *********************************************/

namespace invkl {

  const KLPol& one();

};

/******** type definitions **************************************************/


namespace invkl {

class KLPol:public Polynomial<klsupport::KLCoeff> {
public:
  static klsupport::PolynomialType polType() {return klsupport::INV_KLPOL;}
  KLPol() {};
  KLPol(const Ulong& n):Polynomial<klsupport::KLCoeff>(n) {};
  KLPol(const klsupport::KLCoeff& c, const_tag):Polynomial<klsupport::KLCoeff>(c,const_tag()) {};
  ~KLPol() {};
  KLPol& add(const KLPol& p, const klsupport::KLCoeff& mu, const Ulong& n);
  KLPol& subtract(const KLPol& p, const Ulong& n);
};

struct MuData {
  coxtypes::CoxNbr x;
  klsupport::KLCoeff mu;
  coxtypes::Length height;
/* constructors */
  MuData() {};
  MuData(const coxtypes::CoxNbr& d_x, const klsupport::KLCoeff& d_mu, const coxtypes::Length& d_h)
    :x(d_x), mu(d_mu), height(d_h) {};
  ~MuData() {};
/* comparison */
  bool operator> (const MuData& m) const;                        /* inlined */
};

struct KLStatus {
  static const bits::Lflags kl_done = 1L;
  static const bits::Lflags mu_done = (1L << 1);
  bits::Lflags flags;
  Ulong klrows;
  Ulong klnodes;
  Ulong klcomputed;
  Ulong murows;
  Ulong munodes;
  Ulong mucomputed;
  Ulong muzero;
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLStatus));}
  KLStatus();
  ~KLStatus();
};

class KLContext {
  klsupport::KLSupport* d_klsupport;
  list::List<KLRow*> d_klList;
  list::List<MuRow*> d_muList;
  search::BinaryTree<KLPol> d_klTree;
  KLStatus* d_status;
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
  const bits::BitMap& involution() const;                             /* inlined */
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
  const schubert::SchubertContext& schubert() const;                      /* inlined */
  Ulong size() const;                                         /* inlined */
  const search::BinaryTree<KLPol>& tree() const;                /* inlined */
/* manipulators */
  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a); /* inlined */
  void clearFullKL();                                            /* inlined */
  void clearFullMu();                                            /* inlined */
  void fillKL();
  void fillMu();
  const KLPol& klPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
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

inline bool MuData::operator> (const MuData& m) const {return x > m.x;}

inline const klsupport::ExtrRow& KLContext::extrList(const coxtypes::CoxNbr& y) const
  {return klsupport().extrList(y);}
inline coxtypes::CoxNbr KLContext::inverse(const coxtypes::CoxNbr& x) const
  {return d_klsupport->inverse(x);}
inline const bits::BitMap& KLContext::involution() const
  {return d_klsupport->involution();}
inline bool KLContext::isExtrAllocated(const coxtypes::CoxNbr& x) const
  {return d_klsupport->isExtrAllocated(x);}
inline bool KLContext::isFullKL() const
  {return d_status->flags&KLStatus::kl_done;}
inline bool KLContext::isFullMu() const
  {return d_status->flags&KLStatus::mu_done;}
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
inline const schubert::SchubertContext& KLContext::schubert() const
  {return d_klsupport->schubert();}
inline Ulong KLContext::size() const {return d_klList.size();}
inline const search::BinaryTree<KLPol>& KLContext::tree() const
  {return d_klTree;}

inline void KLContext::applyIPermutation(const coxtypes::CoxNbr& y,
					 const bits::Permutation& a)
  {return rightRangePermute(*d_klList[y],a);}
inline void KLContext::clearFullKL() {d_status->flags &= ~KLStatus::kl_done;}
inline void KLContext::clearFullMu() {d_status->flags &= ~KLStatus::mu_done;}
inline void KLContext::setFullKL() {d_status->flags |= KLStatus::kl_done;}
inline void KLContext::setFullMu() {d_status->flags |= KLStatus::mu_done;}

};

#endif

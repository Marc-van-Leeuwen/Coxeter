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

  using HeckeElt = containers::vector<hecke::HeckeMonomial<kl::KLPol> >;
};

/******** type definitions **************************************************/


namespace invkl {

class KLContext {
  klsupport::KLSupport* d_klsupport;
  list::List<kl::KLRow*> d_klList;
  list::List<kl::MuRow*> d_muList;
  search::BinaryTree<kl::KLPol> d_klTree;
  kl::KLStats* d_status;
  struct KLHelper; // implements everything
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
  const klsupport::ExtrRow& extrList(coxtypes::CoxNbr y) const
    { return klsupport().extrList(y);}
  coxtypes::CoxNbr inverse(coxtypes::CoxNbr x) const
    { return d_klsupport->inverse(x);}
  const bitmap::BitMap& involution() const { return d_klsupport->involution(); }
  bool isExtrAllocated(coxtypes::CoxNbr x) const
    { return d_klsupport->isExtrAllocated(x);}
  bool isFullKL() const {return d_status->flags&kl::KLStats::kl_done;}
  bool isFullMu() const {return d_status->flags&kl::KLStats::mu_done;}
  bool isKLAllocated(coxtypes::CoxNbr x) const { return d_klList[x] != nullptr;}
  bool isMuAllocated(coxtypes::CoxNbr x) const { return d_muList[x] != nullptr;}
  const kl::KLRow& klList(coxtypes::CoxNbr y) const {return *d_klList[y];}
  const klsupport::KLSupport& klsupport() const {return *d_klsupport;}
  coxtypes::Generator last(coxtypes::CoxNbr y) const
    { return d_klsupport->last(y);}
  const kl::MuRow& muList(coxtypes::CoxNbr y) const
    { return *d_muList[y];}
  coxtypes::Rank rank() const {return d_klsupport->rank();}
  const schubert::SchubertContext& schubert() const
    { return d_klsupport->schubert(); }
  Ulong size() const {return d_klList.size();}
  const search::BinaryTree<kl::KLPol>& tree() const { return d_klTree; }
  void applyInverse(coxtypes::CoxNbr y);
  void applyIPermutation(coxtypes::CoxNbr y, const bits::Permutation& a)
  {return right_permute(*d_klList[y],a);}
  void clearFullKL() { d_status->flags &= ~kl::KLStats::kl_done; }
  void clearFullMu() { d_status->flags &= ~kl::KLStats::mu_done; }
  void fillKL();
  void fillMu();
  const kl::KLPol& klPol(coxtypes::CoxNbr x, coxtypes::CoxNbr y,
		     const coxtypes::Generator& s = coxtypes::undef_generator);
  klsupport::KLCoeff mu(coxtypes::CoxNbr x, coxtypes::CoxNbr y,
	     const coxtypes::Generator& s = coxtypes::undef_generator);
  void permute(const bits::Permutation& a);
  void row(HeckeElt& h, coxtypes::CoxNbr y);
  void revertSize(const Ulong& n);
  void setFullKL() { d_status->flags |= kl::KLStats::kl_done; }
  void setFullMu() { d_status->flags |= kl::KLStats::mu_done; }
  void setSize(const Ulong& n);
}; // |class KLContext|

};

#endif

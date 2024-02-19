/*
  This is kl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef KL_H  /* guard against multiple inclusions */
#define KL_H

#include <limits.h>

#include "globals.h"
#include "coxtypes.h"
#include "klsupport.h"
#include "hecke.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "io.h"
#include "list.h"


/******** type declarations *************************************************/

namespace kl {
  class KLContext;
  struct KLStatus;
  struct MuData;
  class MuFilter;

  class KLPol;
  typedef list::List<const KLPol*> KLRow;
  typedef list::List<MuData> MuRow;
  using HeckeElt = containers::vector<hecke::HeckeMonomial<KLPol> >;
};

/******** function declarations *********************************************/

namespace kl {
  void cBasis(HeckeElt& h, const coxtypes::CoxNbr& y, KLContext& kl);
  void extractDufloInvolutions(const KLContext& kl, const bits::Partition& pi,
			       bits::BitMap& b);
  void genericSingularities(HeckeElt& h, const coxtypes::CoxNbr& y, KLContext& kl);
  void ihBetti(schubert::Homology& h, const coxtypes::CoxNbr& y, KLContext& kl);
  const KLPol& one();
  bool isSingular(const HeckeElt& h);
  bool isSingular(const KLRow& row);
  void print(FILE* file, const schubert::Homology& h);
  void printMuTable(FILE* file, const KLContext& kl, const interface::Interface& I);
  void showKLPol(FILE* file, KLContext& kl, const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
		 const interface::Interface& I, const coxtypes::Generator& s = coxtypes::undef_generator);
  void showMu(FILE* file, KLContext& kl, const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
	      const interface::Interface& I);
  void sortByPol(KLRow& row);
};

/******** type definitions **************************************************/

#include "interface.h"
#include "search.h"

namespace kl {

class KLPol : public polynomials::Polynomial<klsupport::KLCoeff>
{
public:
  static klsupport::PolynomialType polType() {return klsupport::KLPOL;}
  KLPol() {};
  KLPol(const Ulong& n):Polynomial<klsupport::KLCoeff>(n) {};
  KLPol(const klsupport::KLCoeff& c, const_tag):Polynomial<klsupport::KLCoeff>(c,const_tag()) {};
  ~KLPol() {};
  // KLPol& add(const KLPol& p, const long& n);
  // KLPol& subtract(const KLPol& p, const MuPol& mp, const Ulong& n);
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

class MuFilter {
 private:
  const schubert::SchubertContext& d_p;
  coxtypes::Length d_l;
 public:
  MuFilter(const schubert::SchubertContext& p, const coxtypes::Length& l);
  MuFilter(const schubert::SchubertContext& p, const coxtypes::CoxNbr& y);
  ~MuFilter();
  bool operator() (const coxtypes::CoxNbr& x) const;                       /* inlined */
};

class KLContext {
 private:
  klsupport::KLSupport* d_klsupport;
  list::List<KLRow*> d_klList;
  list::List<MuRow*> d_muList;
  search::BinaryTree<KLPol> d_klTree;
  KLStatus* d_status;
  struct KLHelper; /* provides helper functions */
  KLHelper* d_help;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLContext));}
  KLContext(klsupport::KLSupport* kls);
  ~KLContext();
/* accessors */
  const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y) const
  { return klsupport().extrList(y); }

  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const;                        /* inlined */
  const bits::BitMap& involution() const;                             /* inlined */
  bool isExtrAllocated(const coxtypes::CoxNbr& x) const
  { return d_klsupport->isExtrAllocated(x); }
  bool isFullKL() const;                                        /* inlined */
  bool isFullMu() const;                                        /* inlined */
  bool isKLAllocated(const coxtypes::CoxNbr& x) const;                    /* inlined */
  bool isMuAllocated(const coxtypes::CoxNbr& x) const;                    /* inlined */
  const KLRow& klList(const coxtypes::CoxNbr& y) const;
  const klsupport::KLSupport& klsupport() const;
  coxtypes::Generator last(const coxtypes::CoxNbr& y) const;                        /* inlined */
  const MuRow& muList(const coxtypes::CoxNbr& y) const;
  coxtypes::Rank rank() const;                                            /* inlined */
  const schubert::SchubertContext& schubert() const
  { return d_klsupport->schubert(); }
  Ulong size() const;
  const search::BinaryTree<KLPol>& tree() const;                        /* inlined */
/* manipulators */
  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a);/* inlined */
  void clearFullKL();                                           /* inlined */
  void clearFullMu();                                           /* inlined */
  void fillKL();
  void fillMu();
  const KLPol& klPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
		     const coxtypes::Generator& s = coxtypes::undef_generator);
  klsupport::KLCoeff mu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);
  void row(HeckeElt& h, const coxtypes::CoxNbr& y);
  void setFullKL();                                             /* inlined */
  void setFullMu();                                             /* inlined */
  void setSize(const Ulong& n);
/* input/output */
  // String& append(String& str, const coxtypes::CoxNbr& x) const;
  void print(FILE* file, const coxtypes::CoxNbr& x, const interface::Interface& I) const;
                                                                /* inlined */
  void printStatus(FILE* file) const;
/* to be taken out! */
  void compareMu();
};

};

/******** inlined definitions **********************************************/

namespace kl {

inline bool MuData::operator> (const MuData& m) const {return x > m.x;}

inline bool MuFilter::operator() (const coxtypes::CoxNbr& x) const
  {coxtypes::Length l = d_p.length(x); return ((d_l-l)%2) && ((d_l-l) > 1);}

inline coxtypes::CoxNbr KLContext::inverse(const coxtypes::CoxNbr& x) const
  {return d_klsupport->inverse(x);}
inline const bits::BitMap& KLContext::involution() const
  {return d_klsupport->involution();}
inline bool KLContext::isFullKL() const
  {return d_status->flags&KLStatus::kl_done;}
inline bool KLContext::isFullMu() const
  {return d_status->flags&KLStatus::mu_done;}
inline bool KLContext::isKLAllocated(const coxtypes::CoxNbr& x) const
  {return d_klList[x] != 0;}
inline bool KLContext::isMuAllocated(const coxtypes::CoxNbr& x) const
  {return d_muList[x] != 0;}
inline coxtypes::Generator KLContext::last(const coxtypes::CoxNbr& y) const
  {return d_klsupport->last(y);}
inline coxtypes::Rank KLContext::rank() const {return d_klsupport->rank();}
inline const search::BinaryTree<KLPol>& KLContext::tree() const
  {return d_klTree;}
inline void KLContext::print
  (FILE* file, const coxtypes::CoxNbr& x, const interface::Interface& I)
  const {schubert().print(file,x,I);}

inline void KLContext::applyIPermutation
  (const coxtypes::CoxNbr& y, const bits::Permutation& a)
  {return rightRangePermute(*d_klList[y],a);}
inline void KLContext::clearFullKL() {d_status->flags &= ~KLStatus::kl_done;}
inline void KLContext::clearFullMu() {d_status->flags &= ~KLStatus::mu_done;}
inline void KLContext::setFullKL() {d_status->flags |= KLStatus::kl_done;}
inline void KLContext::setFullMu() {d_status->flags |= KLStatus::mu_done;}

};

#endif

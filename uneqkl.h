/*
  This is uneqkl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef UNEQKL_H  /* guard against multiple inclusions */
#define UNEQKL_H

#include "globals.h"
#include "coxtypes.h"
#include "hecke.h"
#include "klsupport.h"
#include "list.h"
#include "polynomials.h"
#include "bits.h"
#include "memory.h"
#include "search.h"

/******** type declarations **************************************************/

namespace uneqkl {
  class KLContext;
  class KLPol;
  class MuPol;
  struct KLStatus;
  struct MuData;

  typedef list::List<const KLPol*> KLRow;

  typedef list::List<MuData> MuRow;
  typedef list::List<MuRow*> MuTable;
  typedef list::List<HeckeMonomial<KLPol> > HeckeElt;
};

/******** function declarations **********************************************/

namespace uneqkl {
  void cBasis(HeckeElt& h, const coxtypes::CoxNbr& y, KLContext& kl);
  const MuPol& errorMuPol();
  const KLPol& errorPol();
  const KLPol& one();
  const MuPol& zero();
};

/******** type definitions ***************************************************/

namespace uneqkl {

class KLPol : public Polynomial<klsupport::SKLcoeff> {
  static const klsupport::SKLcoeff min_coeff = klsupport::SKLCOEFF_MIN;
  static const klsupport::SKLcoeff max_coeff = klsupport::SKLCOEFF_MAX;
public:
  static klsupport::PolynomialType polType() {return klsupport::UNEQ_KLPOL;}
  KLPol() {};
  KLPol(const Ulong& n):Polynomial<klsupport::SKLcoeff>(n) {};
  KLPol(const klsupport::SKLcoeff& c, const_tag)
    : Polynomial<klsupport::SKLcoeff>(c,const_tag()) {};
  ~KLPol() {};
  KLPol& add(const KLPol& p, const long& n);
  KLPol& subtract(const KLPol& p, const MuPol& mp, const Ulong& n);
};

class MuPol : public LaurentPolynomial<klsupport::SKLcoeff> {
public:
  struct const_tag {};
  MuPol():LaurentPolynomial<klsupport::SKLcoeff>() {}; // zero
  MuPol(const SDegree& d, const SDegree& o = 0)
    :LaurentPolynomial<klsupport::SKLcoeff>(d,o) {};
  MuPol(const klsupport::SKLcoeff& c, const_tag)
    : LaurentPolynomial<klsupport::SKLcoeff>
    (std::vector<klsupport::SKLcoeff>{c},0)
    {}
  ~MuPol() {};

  static MuPol zero() { return MuPol(); }
};

struct MuData {
  coxtypes::CoxNbr x;
  const MuPol* pol;
/* constructors anc destructors*/
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(MuData));}
  MuData() {};
  MuData(const coxtypes::CoxNbr& x, const MuPol* pol):x(x),pol(pol) {};
/* comparison */
  bool operator> (const MuData& m) const;                        /* inlined */
  bool operator< (const MuData& m) const;                        /* inlined */
  bool operator== (const MuData& m) const;                        /* inlined */
};

struct KLStatus {
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
  KLStatus() {};
  ~KLStatus() {};
};

class KLContext {
  klsupport::KLSupport* d_klsupport;
  list::List<KLRow*> d_klList;
  list::List<MuTable*> d_muTable;
  list::List<coxtypes::Length> d_L; /* lengths of generators */
  list::List<coxtypes::Length> d_length; /* lengths of context elements */
  search::BinaryTree<KLPol> d_klTree;
  search::BinaryTree<MuPol> d_muTree;
  KLStatus* d_status;
  struct KLHelper; /* provides helper functions */
  KLHelper* d_help;
  friend struct KLHelper;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLContext));}
  KLContext(klsupport::KLSupport* kls, const graph::CoxGraph& G, const interface::Interface& I);
  ~KLContext();
/* accessors */
  const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y) const;                /* inlined */
  Ulong genL(const coxtypes::Generator& s) const;                        /* inlined */
  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const;                         /* inlined */
  bool isKLAllocated(const coxtypes::CoxNbr& y) const;                     /* inlined */
  bool isMuAllocated(const coxtypes::Generator& s, const coxtypes::CoxNbr& y) const; /* inlined */
  const KLRow& klList(const coxtypes::CoxNbr& y) const;                    /* inlined */
  const klsupport::KLSupport& klsupport() const;                            /* inlined */
  coxtypes::Generator last(const coxtypes::CoxNbr& x) const;                         /* inlined */
  Ulong length(const coxtypes::CoxNbr& x) const;                         /* inlined */
  const MuRow& muList(const coxtypes::Generator& s, const coxtypes::CoxNbr& y) const;/* inlined */
  coxtypes::Rank rank() const;                                             /* inlined */
  const schubert::SchubertContext& schubert() const;                       /* inlined */
  Ulong size() const;                                          /* inlined */
/* modifiers */
  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a); /* inlined */
  void fillKL();
  void fillMu();
  void fillMu(const coxtypes::Generator& s);
  const KLPol& klPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  const MuPol mu(const coxtypes::Generator& s,
		 const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  void row(HeckeElt& h, const coxtypes::CoxNbr& y);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);
  void setSize(const Ulong& n);
};

};

/******** inline definitions ************************************************/

namespace uneqkl {

inline bool MuData::operator> (const MuData& m) const {return x > m.x;}
inline bool MuData::operator< (const MuData& m) const {return x < m.x;}
inline bool MuData::operator== (const MuData& m) const {return x == m.x;}

inline const klsupport::ExtrRow& KLContext::extrList
  (const coxtypes::CoxNbr& y) const
  {return klsupport().extrList(y);}
inline Ulong KLContext::genL(const coxtypes::Generator& s) const
  {return d_L[s];}
inline coxtypes::CoxNbr KLContext::inverse(const coxtypes::CoxNbr& x) const
  {return klsupport().inverse(x);}
inline bool KLContext::isKLAllocated(const coxtypes::CoxNbr& y) const
  {return d_klList[y] != 0;}
inline bool KLContext::isMuAllocated
  (const coxtypes::Generator& s, const coxtypes::CoxNbr& y) const
  {return (*d_muTable[s])[y] != 0;}
inline const KLRow& KLContext::klList(const coxtypes::CoxNbr& y) const
  {return *d_klList[y];}
inline const klsupport::KLSupport& KLContext::klsupport() const
  {return *d_klsupport;}
inline coxtypes::Generator KLContext::last(const coxtypes::CoxNbr& x) const
  {return klsupport().last(x);}
inline Ulong KLContext::length(const coxtypes::CoxNbr& x) const
  {return d_length[x];}
inline const MuRow& KLContext::muList
  (const coxtypes::Generator& s, const coxtypes::CoxNbr& y)
  const {return d_muTable[s][0][y][0];}
inline coxtypes::Rank KLContext::rank() const {return d_klsupport->rank();}
inline const schubert::SchubertContext& KLContext::schubert() const
  {return klsupport().schubert();}
inline Ulong KLContext::size() const {return d_klList.size();}

inline void KLContext::applyIPermutation(const coxtypes::CoxNbr& y,
					 const bits::Permutation& a)
  {return rightRangePermute(*d_klList[y],a);}

};

#endif

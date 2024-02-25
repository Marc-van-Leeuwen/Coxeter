/*
  This is uneqkl.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef UNEQKL_H  /* guard against multiple inclusions */
#define UNEQKL_H

#include "globals.h"
#include "containers.h"
#include "coxtypes.h"
#include "hecke.h"
#include "klsupport.h"
#include "polynomials.h"
#include "bits.h"
#include "memory.h"

/******** type declarations **************************************************/

namespace uneqkl {
  class KLContext;
  class KLPol;
  class MuPol;
  struct MuData;

  using KLRow = containers::vector<const KLPol*>; // non-owning pointers
  using KLRowPtr = std::unique_ptr<KLRow>;
  using KLTable = containers::vector<KLRowPtr>;

  using MuRow = containers::vector<MuData>;
  using MuRowPtr = std::unique_ptr<MuRow>; // half of the time |nullptr|
  using MuTable = containers::vector<MuRowPtr>;

  using HeckeElt = containers::vector<hecke::HeckeMonomial<KLPol> >;
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

class KLPol : public polynomials::Polynomial<klsupport::SKLcoeff>
{
  static const klsupport::SKLcoeff min_coeff = klsupport::SKLCOEFF_MIN;
  static const klsupport::SKLcoeff max_coeff = klsupport::SKLCOEFF_MAX;
public:
  static klsupport::PolynomialType polType() { return klsupport::UNEQ_KLPOL; }
  KLPol() {};
  KLPol(const Ulong& n):Polynomial<klsupport::SKLcoeff>(n) {};
  KLPol(const klsupport::SKLcoeff& c, const_tag)
    : Polynomial<klsupport::SKLcoeff>(c,const_tag()) {};
  ~KLPol() {};
  KLPol& add(const KLPol& p, const long& n);
  KLPol& subtract(const KLPol& p, const MuPol& mp, const Ulong& n);
}; // |class KLPol|

struct MuPol : public polynomials::LaurentPolynomial<klsupport::SKLcoeff>
{ // a very thin wrapper around its base, mainly to provide an extra constructor
public:
  struct const_tag {}; // to enable the constructor of constant polynomials
  MuPol() : LaurentPolynomial<klsupport::SKLcoeff>() {} // zero
  MuPol(const polynomials::SDegree& d, const polynomials::SDegree& o = 0)
    : LaurentPolynomial<klsupport::SKLcoeff>(d,o) {} // zero-terms in range:
  MuPol(const klsupport::SKLcoeff& c, const_tag) // the new constructor
    : LaurentPolynomial<klsupport::SKLcoeff>
    (containers::vector<klsupport::SKLcoeff>{c},0)
    {}
  ~MuPol() {};

  static MuPol zero() { return MuPol(); }
}; // |struct MuPol|

struct MuData
{
  coxtypes::CoxNbr x;
  const MuPol* pol; // non owning pointer (polynomials owned by |KLContext|)
/* constructors and destructors*/
  MuData() {}
  MuData(const coxtypes::CoxNbr& x, const MuPol* pol) : x(x), pol(pol) {}
/* comparison */
  bool operator>  (const MuData& m) const { return x >  m.x; }
  bool operator<  (const MuData& m) const { return x <  m.x; }
  bool operator== (const MuData& m) const { return x == m.x; }
}; // |struct MuData|

struct KLStats // holds statistics
{
  Ulong klrows;
  Ulong klnodes;
  Ulong klcomputed;
  Ulong murows;
  Ulong munodes;
  Ulong mucomputed;
  Ulong muzero;
/* constructors and destructors */
  KLStats() : klrows(0), klnodes(0), klcomputed(0)
	    , murows(0), munodes(0), mucomputed(0), muzero(0)
            {}
}; // |struct KLStats|

// the main class for this module
class KLContext
{
  struct KLHelper; /* provides helper functions */
  KLHelper* d_help; // pointer level hides implementation
 public:
/* constructors and destructors */
  void* operator new(size_t size)
    { return memory::arena().alloc(size); }
  void operator delete(void* ptr)
    { return memory::arena().free(ptr,sizeof(KLContext)); }
  KLContext(klsupport::KLSupport& kls, const graph::CoxGraph& G,
	    const interface::Interface& I);
  ~KLContext();
/* accessors */
  Ulong size() const;
  const klsupport::KLSupport& klsupport() const;
  coxtypes::Rank rank() const {return klsupport().rank();}
  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const
    { return klsupport().inverse(x); }
  const schubert::SchubertContext& schubert() const
    { return klsupport().schubert(); }
  const KLRow& klList(const coxtypes::CoxNbr& y) const;
  const MuRow& muList
  (const coxtypes::Generator& s, const coxtypes::CoxNbr& y) const;
  const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y) const
    { return klsupport().extrList(y); }

// manipulators
  void setSize(const Ulong& n);
  void revertSize(const Ulong& n);
  // the following are 5 methods entry points, called from |CoxGroup| methods
  // all are deemed modifiers because they extend the storage pools
  void fillKL();
  const KLPol& klPol(coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  void fillMu();
  const MuPol mu(const coxtypes::Generator& s,
		 const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  void row(HeckeElt& h, const coxtypes::CoxNbr& y);

  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a);
  void permute(const bits::Permutation& a);

}; // |class KLContext|

}; // |namespace uneqkl|

#endif

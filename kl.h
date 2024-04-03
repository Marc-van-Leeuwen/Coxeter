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
#include "containers.h"
#include "klsupport.h"
#include "hecke.h"
#include "polynomials.h"
#include "bits.h"
#include "io.h"


/******** type declarations *************************************************/

namespace kl {
  class KLContext;
  struct KLStatus;
  struct MuData;
  class MuFilter;

  class KLPol;
  using KLRow = containers::vector<const KLPol*>;
  using KLRowPtr = std::unique_ptr<KLRow>;
  using KLTable = containers::vector<KLRowPtr>;

  using MuRow = containers::vector<MuData>;
  using MuRowPtr = std::unique_ptr<MuRow>; // half of the time |nullptr|
  using MuTable = containers::vector<MuRowPtr>;
  using HeckeElt = containers::vector<hecke::HeckeMonomial<KLPol> >;
};

/******** type definitions **************************************************/

#include "interface.h"

namespace kl {

class KLPol : public polynomials::Polynomial<klsupport::KLCoeff>
{
public:
  static klsupport::PolynomialType polType() {return klsupport::KLPOL;}
  KLPol() {}
  KLPol(const Ulong& n):Polynomial<klsupport::KLCoeff>(n) {}
  KLPol(const klsupport::KLCoeff& c, const_tag)
    : Polynomial<klsupport::KLCoeff>(c,const_tag()) {}

  // KLPol& add(const KLPol& p, const long& n);
  // KLPol& subtract(const KLPol& p, const MuPol& mp, const Ulong& n);
};

struct KLStats {
  static constexpr unsigned char kl_done = 1L;
  static constexpr unsigned char mu_done = (1L << 1);
  Ulong klrows;
  Ulong klnodes;
  Ulong klcomputed;
  Ulong murows;
  Ulong munodes;
  Ulong mucomputed;
  Ulong muzero;
  unsigned char flags;
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    { return memory::arena().free(ptr,sizeof(KLStats)); }
  KLStats();
};

struct MuData {
  coxtypes::CoxNbr x;
  klsupport::KLCoeff mu;
  coxtypes::Length height;
/* constructors */
  MuData(coxtypes::CoxNbr d_x,
	 const klsupport::KLCoeff& d_mu,
	 const coxtypes::Length& d_h)
    :x(d_x), mu(d_mu), height(d_h)
  {}
/* comparison */
  bool operator< (const MuData& m) const { return x < m.x; }

};

class MuFilter {
 private:
  const schubert::SchubertContext& d_p;
  coxtypes::Length d_l;
 public:
  MuFilter(const schubert::SchubertContext& p, const coxtypes::Length& l);
  MuFilter(const schubert::SchubertContext& p, coxtypes::CoxNbr y);

  bool operator() (coxtypes::CoxNbr x) const
  { coxtypes::Length d = d_l-d_p.length(x); return d%2!=0 and d > 1; }
};

class KLContext {
 private:
  struct KLHelper; /* provides helper functions */
  KLHelper* d_help;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLContext));}
  KLContext(klsupport::KLSupport* kls);

  // relay methods
  const klsupport::KLSupport& klsupport() const; // uses reference in helper
  const klsupport::ExtrRow& extrList(coxtypes::CoxNbr y) const
  { return klsupport().extrList(y); }
  coxtypes::Generator last(coxtypes::CoxNbr y) const
  { return klsupport().last(y); }
  coxtypes::CoxNbr inverse(coxtypes::CoxNbr x) const
  { return klsupport().inverse(x); }
  const bitmap::BitMap& involution() const { return klsupport().involution(); }
  const schubert::SchubertContext& schubert() const
  { return klsupport().schubert(); }

  // sometimes we need non-const access to the support class or |SchubertContext|
  klsupport::KLSupport& klsupport(); // uses reference in helper
  schubert::SchubertContext& schubert() { return klsupport().schubert(); }

  void print
    (FILE* file, coxtypes::CoxNbr x, const interface::Interface& I) const
  { schubert::print(schubert(),file,x,I); }

  const KLStats& stats() const;

  // accessors
  Ulong size() const;
  coxtypes::Rank rank() const { return klsupport().rank(); }
  const KLRow& klList(coxtypes::CoxNbr y) const;
  const MuRow& muList(coxtypes::CoxNbr y) const;

  // manipulators that may expand/shrink tables as only side effect
  const KLPol& klPol
    (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s);
  const KLPol& klPol (coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  klsupport::KLCoeff mu(coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  void row(HeckeElt& h, coxtypes::CoxNbr y);
  void fillKL();
  void fillMu();
  void setSize(const Ulong& n);
  void revertSize(const Ulong& n);

  // manipulators that drastically alter the state
  void applyInverse(coxtypes::CoxNbr y); // more row of |x| to its inverse
  void applyIPermutation(coxtypes::CoxNbr y, const bits::Permutation& a);
  void permute(const bits::Permutation& a);

}; // |class KLContext|

/******** function declarations *********************************************/

  HeckeElt cBasis(coxtypes::CoxNbr y, KLContext& kl);
  void genericSingularities
    (HeckeElt& h, coxtypes::CoxNbr y, KLContext& kl);
  schubert::Homology ihBetti(coxtypes::CoxNbr y, KLContext& kl);
  const KLPol& one();
  bool isSingular(const HeckeElt& h);
  bool isSingular(const KLRow& row);
  void print_stats(const KLContext& kl, FILE* f); // unused
  void print(FILE* file, const schubert::Homology& h);
  void printMuTable
    (FILE* file, const KLContext& kl, const interface::Interface& I);
  void showKLPol
    (FILE* file, KLContext& kl,
     coxtypes::CoxNbr x, coxtypes::CoxNbr y,
     const interface::Interface& I,
     coxtypes::Generator s = coxtypes::undef_generator);
  void showMu
    (FILE* file,
     KLContext& kl, coxtypes::CoxNbr x, coxtypes::CoxNbr y,
     const interface::Interface& I);
  void sortByPol(KLRow& row);

}; // |namespace kl|

#endif

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

struct KLStats {
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
    { return memory::arena().free(ptr,sizeof(KLStats)); }
  KLStats();
};

struct MuData {
  coxtypes::CoxNbr x;
  klsupport::KLCoeff mu;
  coxtypes::Length height;
/* constructors */
  MuData() {}
  MuData(const coxtypes::CoxNbr& d_x,
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
  MuFilter(const schubert::SchubertContext& p, const coxtypes::CoxNbr& y);
  ~MuFilter();
  bool operator() (const coxtypes::CoxNbr& x) const
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

/* accessors */
  const klsupport::KLSupport& klsupport() const; // uses reference in helper
  KLStats& stats();
  const KLStats& stats() const;

  const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y) const
  { return klsupport().extrList(y); }

  coxtypes::Generator last(const coxtypes::CoxNbr& y) const
  { return klsupport().last(y); }
  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const
  { return klsupport().inverse(x); }
  const bits::BitMap& involution() const { return klsupport().involution(); }
  bool isExtrAllocated(const coxtypes::CoxNbr& x) const
  { return klsupport().isExtrAllocated(x); }
  coxtypes::Rank rank() const { return klsupport().rank(); }
 const schubert::SchubertContext& schubert() const
  { return klsupport().schubert(); }
  const KLRow& klList(const coxtypes::CoxNbr& y) const;
  const MuRow& muList(const coxtypes::CoxNbr& y) const;
  Ulong size() const;
/* manipulators */
  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a);
  void fillKL();
  void fillMu();
  const KLPol& klPol
    (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s);
  const KLPol& klPol (coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  klsupport::KLCoeff mu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);
  void row(HeckeElt& h, const coxtypes::CoxNbr& y);
  void setSize(const Ulong& n);

  bool isFullKL() const;
  bool isFullMu() const;
  void clearFullKL();
  void clearFullMu();
/* input/output */
  // String& append(String& str, const coxtypes::CoxNbr& x) const;
  void print
    (FILE* file, const coxtypes::CoxNbr& x, const interface::Interface& I)
  const { schubert().print(file,x,I); }
  void printStatus(FILE* file) const;
};

};

/******** inlined definitions **********************************************/

namespace kl {



};

#endif

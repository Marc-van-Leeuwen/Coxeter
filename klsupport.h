/*
  This is klsupport.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef KLSUPPORT_H  /* guard against multiple inclusions */
#define KLSUPPORT_H

#include <memory>
#include "globals.h"

#include "coxtypes.h"
#include "list.h"
#include "bits.h" // for |bits::BitMap|
#include "bitmap.h" // for |bitmap::BitMap|
#include "polynomials.h"
#include "schubert.h"

/******** type declarations **************************************************/

namespace klsupport {
  class KLSupport;

  typedef unsigned short KLCoeff;
  typedef short SKLcoeff;
  using ExtrRow = containers::vector<coxtypes::CoxNbr>;
};

/******** constants **********************************************************/

namespace klsupport {
  enum PolynomialType {KLPOL, UNEQ_KLPOL, INV_KLPOL, NUM_POLTYPES};

  const KLCoeff KLCOEFF_MAX = USHRT_MAX-1; /* top value is reserved */
  const KLCoeff undef_klcoeff = KLCOEFF_MAX + 1;
  const KLCoeff KLCOEFF_MIN = 0;
  const SKLcoeff SKLCOEFF_MIN = SHRT_MIN+1;
  const SKLcoeff SKLCOEFF_MAX = -SKLCOEFF_MIN;
  const SKLcoeff undef_sklcoeff = SKLCOEFF_MIN-1;
};

/******** function declarations **********************************************/

namespace klsupport {
  KLCoeff& safeAdd(KLCoeff& a, const KLCoeff& b);
  SKLcoeff& safeAdd(SKLcoeff& a, const SKLcoeff& b);
  KLCoeff& safeMultiply(KLCoeff& a, const KLCoeff& b);
  SKLcoeff& safeMultiply(SKLcoeff& a, const SKLcoeff& b);
  KLCoeff& safeSubtract(KLCoeff& a, const KLCoeff& b);
};

/******** type definitions ***************************************************/


namespace klsupport {

class KLSupport {
 private:
  std::unique_ptr<schubert::SchubertContext> d_schubert; // owned
  containers::vector<std::unique_ptr<ExtrRow> > d_extrList;
  containers::vector<coxtypes::CoxNbr> d_inverse;
  // store final letter in standard word: a (nicely distibuted) right descent:
  containers::vector<coxtypes::Generator> d_last;
  bitmap::BitMap d_involution;
  bitmap::BitMap recursively_allocated;
 public:
  containers::vector<std::unique_ptr<ExtrRow> >& extrList()
    { return d_extrList; } // this should go
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLSupport));}
  KLSupport(std::unique_ptr<schubert::SchubertContext> p); // takes ownership
  ~KLSupport();
/* accessors */
  // relay accessors
  coxtypes::Length length(coxtypes::CoxNbr x) const
    { return d_schubert->length(x); }
  coxtypes::Rank rank() const {return d_schubert->rank(); }
  coxtypes::CoxNbr size() const { return d_schubert->size(); }

  // inlined accessors
  const schubert::SchubertContext& schubert() const { return *d_schubert; }
  coxtypes::Generator last(coxtypes::CoxNbr x) const
  { return d_last[x]; } // rightmost letter in ShortLex minimal word for |x|
  const ExtrRow& extr_list(coxtypes::CoxNbr y)
    { if (d_extrList[y]==nullptr) generate_extr_list(y); return *d_extrList[y]; }
  const ExtrRow& extrList(coxtypes::CoxNbr y) const
    { return *d_extrList[y]; } // |const| version, which supposes presence
  bool isExtrAllocated(coxtypes::CoxNbr x) const
    { return d_extrList[x] != nullptr; }
  coxtypes::CoxNbr inverse(coxtypes::CoxNbr x) const
    { return d_inverse[x]; }
  const bitmap::BitMap& involution() const { return d_involution; }
  bool isInvolution(coxtypes::CoxNbr x) const
    { return d_involution.is_member(x); } // whether |x==inverse(x)|

  containers::vector<coxtypes::Generator> standard_path (coxtypes::CoxNbr x)
    const;
/* manipulators */
  void ensure_extr_rows_for(coxtypes::CoxNbr y);
  void move_extr_list_from_inverse(coxtypes::CoxNbr y);
  void applyIPermutation(coxtypes::CoxNbr y, const bits::Permutation& a)
    { bits::right_permute (*d_extrList[y],a); }
  coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);

private:
  void generate_extr_list(coxtypes::CoxNbr y);
}; // |class KLSupport|

}; // namespace sklupport|

#endif

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
  schubert::SchubertContext* d_schubert;
  containers::vector<std::unique_ptr<ExtrRow> > d_extrList;
  list::List<coxtypes::CoxNbr> d_inverse;
  list::List<coxtypes::Generator> d_last;
  bits::BitMap d_involution;
 public:
  containers::vector<std::unique_ptr<ExtrRow> >& extrList()
    { return d_extrList; } // this should go
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(KLSupport));}
  KLSupport(schubert::SchubertContext* p);
  ~KLSupport();
/* accessors */
  const ExtrRow& extrList(const coxtypes::CoxNbr& y) const;       /* inlined */
  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& x) const;                         /* inlined */
  coxtypes::CoxNbr inverseMin(const coxtypes::CoxNbr& x) const;
  const bits::BitMap& involution() const;                         /* inlined */
  bool isExtrAllocated(const coxtypes::CoxNbr& x) const
  { return d_extrList[x] != nullptr; }
  bool isInvolution(const coxtypes::CoxNbr& x) const;             /* inlined */
  coxtypes::Generator last(const coxtypes::CoxNbr& x) const;      /* inlined */
  coxtypes::Length length(const coxtypes::CoxNbr& x) const;       /* inlined */
  coxtypes::Rank rank() const;                                    /* inlined */
  const schubert::SchubertContext& schubert() const { return *d_schubert; }
  coxtypes::CoxNbr size() const { return d_schubert->size(); }
  void sortIRow(const coxtypes::CoxNbr& y, bits::Permutation& a) const; /* inlined */
  void standardPath(list::List<coxtypes::Generator>& g, const coxtypes::CoxNbr& x) const;
/* manipulators */
  void allocExtrRow(const coxtypes::CoxNbr& y);
  void allocRowComputation(const coxtypes::CoxNbr& y);
  void applyInverse(const coxtypes::CoxNbr& y);
  void applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a); /* inlined */
  coxtypes::CoxNbr extendContext(const coxtypes::CoxWord& g);
  void permute(const bits::Permutation& a);
  void revertSize(const Ulong& n);
  schubert::SchubertContext& schubert() { return *d_schubert; }
};

};

/******** inlined definitions ************************************************/

namespace klsupport {

inline const ExtrRow& KLSupport::extrList(const coxtypes::CoxNbr& y) const
  {return *d_extrList[y];}
inline coxtypes::CoxNbr KLSupport::inverse(const coxtypes::CoxNbr& x) const {return d_inverse[x];}
inline const bits::BitMap& KLSupport::involution() const {return d_involution;}
inline bool KLSupport::isInvolution(const coxtypes::CoxNbr& x) const
   {return d_involution.getBit(x);}
inline coxtypes::Generator KLSupport::last(const coxtypes::CoxNbr& x) const {return d_last[x];}
inline coxtypes::Length KLSupport::length(const coxtypes::CoxNbr& x) const
  {return d_schubert->length(x);}
inline coxtypes::Rank KLSupport::rank() const {return d_schubert->rank();}

inline void KLSupport::applyIPermutation(const coxtypes::CoxNbr& y, const bits::Permutation& a)
{ bits::right_permute (*d_extrList[y],a);}

};

#endif

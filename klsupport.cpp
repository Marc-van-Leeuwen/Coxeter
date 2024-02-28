/*
  This is klsupport.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "klsupport.h"

/*
  This module contains code for the operations that can be "factored out"
  amongst the various K-L tables : the ordinary one, the inverse one, and
  the one with unequal parameters. Foremost, this is the extremal list.
  In all these cases, K-L polynomials can be readily reduced to the case
  of "extremal pairs".
*/

/*****************************************************************************

        Chapter I -- The KLSupport class.

  This class can be seen as an extension of the schubert context, oriented
  towards K-L computations. Its main function is to construct and maintain
  extrList. Recall that this is a list of rows of CoxNbr. The row for y is
  never allocated if inverse[y] < y (this becomes a bit cumbersome, and
  I'm really tempted to drop it! the main reason for keeping it is that
  it really speeds up and reduces the computations by pushing y down faster,
  by an important factor (at least two in my trials.))

  If the row is allocated, then it always contains the full list of extremal
  pairs for y, i.e. those x <= y for which LR(x) contains LR(y).

  The other stuff is trivial : inverse is the (partially defined) table of
  inverses, last is the last term in the internal normal form (also kept
  because it provides a better descent strategy), involution flags the
  involutions in the context.

  The interface field is here because it was in KLContext, but it should
  probably go away, and be present as a parameter in the functions that
  need it.

  The idea in this program is that schubert is really an auxiliary to
  klsupport; the updating of schubert should go only through klsupport,
  and the two will always be kept consistent.

  The following functions are provided :

   - constructors and destructors :

    - KLSupport(SchubertContext*, const Interface*);
    - ~KLSupport();

   - accessors :

    - inverseMin(const CoxNbr& x) : returns min(x,x_inverse);
    - standardPath(List<Generator>& g, const CoxNbr& x) : returns in g the
      standard descent path from x

   - manipulators :

    - applyInverse(const CoxNbr& y) : writes the row from inverse in y;
    - applyIPermutation(const CoxNbr& y, const Permutation& a) : applies
      the (inverse) permutation to extrList(y); (inlined)
    - extendContext(const CoxWord& g) : extends the context to hold g;
    - permute(const Permutation& a) : applies the permutation to the data;
    - revertSize(const Ulong& n) : reverts to previous size after a
      failed extension;

 *****************************************************************************/

namespace klsupport {

KLSupport::KLSupport(std::unique_ptr<schubert::SchubertContext> p)
  : d_schubert(std::move(p))
  , d_extrList(1)
  , d_inverse(1)
  , d_last(1)
  , d_involution(1)

{
  /* make first row of d_extrList, for the identity element */

  d_extrList[0].reset(new ExtrRow{0});

  d_inverse.setSizeValue(1);

  d_last.setSizeValue(1);
  d_last[0] = coxtypes::undef_generator;

  d_involution.setBit(0);
}


  // |KLSupport| has no owning raw pointers
  KLSupport::~KLSupport() {}

/******** accessors *********************************************************/

coxtypes::CoxNbr KLSupport::inverseMin(const coxtypes::CoxNbr& x) const

/*
  Returns the minimum of x and x_inverse.
*/

{
  if (x <= inverse(x))
    return x;
  else
    return inverse(x);
}

void KLSupport::standardPath(list::List<coxtypes::Generator>& g, const coxtypes::CoxNbr& x) const

{
  const schubert::SchubertContext& p = schubert();

  /* find sequence of shifts */

  coxtypes::Length j = p.length(x);
  g.setSize(j);
  coxtypes::CoxNbr x1 = x;

  while (j) {
    --j;
    if (inverse(x1) < x1) { /* left shift */
      coxtypes::Generator s = last(inverse(x1));
      g[j] = s + rank();
      x1 = p.lshift(x1,s);
    }
    else {
      coxtypes::Generator s = last(x1);
      g[j] = last(x1);
      x1 = p.rshift(x1,s);
    }
  }

  return;
}

/******** manipulators ******************************************************/


/*
  Allocates one row in d_extrList. The row contains the list of elements
  x <= y s.t. LR(y) is contained in LR(x), i.e., which cannot be taken
  further up by the application of a generator in LR(y).

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
*/
void KLSupport::ensure_extr_row_exists(const coxtypes::CoxNbr& y)
{
  if (d_extrList[y]!=nullptr)
    return;
  const schubert::SchubertContext& p = schubert();
  bits::BitMap b(size()); // take full scurrent size of the |KLSupport| class

  p.extractClosure(b,y); // select the Bruhat interval up to |y|
  if (error::ERRNO)
    return;

  schubert::select_maxima_for(p,b,p.descent(y));

  d_extrList[y].reset // fill with the "contents" of the |BitMap b|, increasingly
    (new ExtrRow(b.begin(),b.end()));
}


/*
  This function makes sure that all the extremal rows in the standard
  descent path from y are allocated. The idea is that all these rows
  will come up when the full row for y is computed, so one might as
  well fill them anyway; doing them all at the same time will save
  many Bruhat closure computations, which are relatively expensive.
  Still, this function looks like overkill to me. I'm leaving it in
  because it is working and it was a pain to write! [dixit Fokko]

  Things wouldn't be so bad if there wasn't also the passage to inverses!
*/
void KLSupport::allocRowComputation(const coxtypes::CoxNbr& y)
{
  static list::List<coxtypes::Generator> e(0);
  const schubert::SchubertContext& p = schubert();

  /* find sequence of shifts */

  standardPath(e,y);

  bits::SubSet q(size());

  q.reset();
  q.add(0);
  if (error::ERRNO)
    goto abort;

  {

    coxtypes::CoxNbr y1 = 0;

    for (Ulong j = 0; j < e.size(); ++j) {

      coxtypes::Generator s = e[j];
      p.extendSubSet(q,s);  /* extend the subset */
      if (error::ERRNO)
	goto abort;

      y1 = p.shift(y1,s);
      coxtypes::CoxNbr y2 = inverseMin(y1);

      if (not isExtrAllocated(y2)) { /* allocate row */

	/* copy bitmap of subset to buffer */

	bits::BitMap b = q.bitMap();
	if (error::ERRNO)
	  goto abort;

	/* extremalize */

	schubert::select_maxima_for(p,b,p.descent(y1));
	d_extrList[y1].reset(new ExtrRow(b.begin(),b.end()));

	/* go over to inverses if necessary */

	if (s >= rank()) {
	  applyInverse(y2);
	  std::sort(d_extrList[y2]->begin(),d_extrList[y2]->end());
	}
      }
    }

  }

  return;
 abort:
  error::Error(error::ERRNO);
  error::ERRNO = error::ERROR_WARNING;
  return;
}


/*
  This function moves the contents of |d_extrList[yi]| to |d_extrList[y]|
  (where |yi| is the inverse of |y|) while taking the inverses of all entries.
  The resulting row is not sorted (this can be done with sortIRow).
*/
void KLSupport::applyInverse(const coxtypes::CoxNbr& y)
{
  coxtypes::CoxNbr yi = inverse(y);
  d_extrList[y] = std::move(d_extrList[yi]);

  ExtrRow& e = *d_extrList[y];
  for (Ulong j = 0; j < e.size(); ++j) {
    e[j] = inverse(e[j]);
  }
}


/*
  Extends the context to accomodate g.

  The return value is the context number of g in case of success, and
  undef_coxnbr in case of failure.

  Forwards the error EXTENSION_FAIL in case of error.
*/
coxtypes::CoxNbr KLSupport::extendContext(const coxtypes::CoxWord& g)
{
  coxtypes::CoxNbr prev_size = size();
  schubert::SchubertContext& p = *d_schubert;

  coxtypes::CoxNbr x = p.extendContext(g); // this increases |size()|

  if (error::ERRNO) /* error::ERRNO is EXTENSION_FAIL */
    return coxtypes::undef_coxnbr;

  error::CATCH_MEMORY_OVERFLOW = true;

  d_extrList.resize(size()); // extend with values |std::unique_ptr<>(nullptr)|
  if (error::ERRNO)
    goto revert;
  d_inverse.setSize(size());
  if (error::ERRNO)
    goto revert;
  d_last.setSize(size());
  if (error::ERRNO)
    goto revert;
  d_involution.setSize(size());
  if (error::ERRNO)
    goto revert;

  error::CATCH_MEMORY_OVERFLOW = false;

  /* extend the list of inverses */

  for (coxtypes::CoxNbr x = 0; x < prev_size; ++x) {
    if (inverse(x) == coxtypes::undef_coxnbr) { /* try to extend */
      coxtypes::Generator s = p.firstRDescent(x);
      coxtypes::CoxNbr xs = p.rshift(x,s);
      if (inverse(xs) == coxtypes::undef_coxnbr)
	d_inverse[x] = coxtypes::undef_coxnbr;
      else
	d_inverse[x] = p.lshift(inverse(xs),s);
    }
  }

  // play it again, Sam
  for (coxtypes::CoxNbr x = prev_size; x < size(); ++x) {
    coxtypes::Generator s = p.firstRDescent(x);
    coxtypes::CoxNbr xs = p.rshift(x,s);
    if (inverse(xs) == coxtypes::undef_coxnbr)
      d_inverse[x] = coxtypes::undef_coxnbr;
    else
      d_inverse[x] = p.lshift(inverse(xs),s);
  }

  for (coxtypes::CoxNbr x = prev_size; x < size(); ++x) {
    if (inverse(x) == x)
      d_involution.setBit(x);
  }

  /* extend list of last elements */

  for (coxtypes::CoxNbr x = prev_size; x < size(); ++x) {
    coxtypes::Generator s = p.firstLDescent(x);
    coxtypes::CoxNbr sx = p.lshift(x,s);
    if (sx!=0)
      d_last[x] = d_last[sx];
    else /* x = s */
      d_last[x] = s;
  }

  return x;

 revert:
  error::CATCH_MEMORY_OVERFLOW = false;
  revertSize(prev_size);

  return coxtypes::undef_coxnbr;
}

void KLSupport::permute(const bits::Permutation& a)

/*
  Applies the permutation a to the data in the context. The meaning of a
  is that it takes element number x in the context to element number a(x).

  The procedure is explained in full in kl.h.
*/

{
  /* permute schubert context */

  d_schubert->permute(a);

  /* permute values */

    for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
    if (d_extrList[y] == nullptr)
      continue;
    ExtrRow& e = *d_extrList[y];
    for (Ulong j = 0; j < e.size(); ++j) {
      e[j] = a[e[j]];
    }
  }

  for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) != coxtypes::undef_coxnbr)
      d_inverse[y] = a[inverse(y)];
  }

  /* permute ranges */

  bits::BitMap b(a.size());

  for (coxtypes::CoxNbr x = 0; x < size(); ++x) {
    if (b.getBit(x))
      continue;
    if (a[x] == x) {
      b.setBit(x);
      continue;
    }
    for (coxtypes::CoxNbr y = a[x]; y != x; y = a[y]) {
      d_extrList[x].swap(d_extrList[y]);

      /* back up values for y */
      coxtypes::CoxNbr inverse_buf = inverse(y);
      coxtypes::Generator last_buf = last(y);
      bool involution_buf = isInvolution(y);
      /* put values for x in y */
      d_inverse[y] = inverse(x);
      d_last[y] = last(x);
      if (isInvolution(x))
	d_involution.setBit(y);
      else
	d_involution.clearBit(y);
      /* store backup values in x */
      d_inverse[x] = inverse_buf;
      d_last[x] = last_buf;
      if (involution_buf)
	d_involution.setBit(x);
      else
	d_involution.clearBit(x);
      /* set bit*/
      b.setBit(y);
    }

    b.setBit(x);
  }

  return;
}


/*
  This function reverts the size of the context to a previous (smaller) value n.
  Note that the allocated sizes of the lists are not changed; we simply preserve
  the consistency of the various size values.
*/

void KLSupport::revertSize(const Ulong& n)
{
  d_schubert->revertSize(n);
  d_extrList.resize(n);
  d_inverse.setSize(n);
  d_last.setSize(n);

  return;
}

};

/*****************************************************************************

        Chapter II -- Utilities.

  This section defines some utility functions declared in klsupport.h :

    - safeAdd(const KLCoeff&, const KLCoeff&) : safe addition;
    - safeAdd(SKLcoeff&, const SKLcoeff&) : safe addition;

 *****************************************************************************/

namespace klsupport {

KLCoeff& safeAdd(KLCoeff& a, const KLCoeff& b)

/*
  This function increments a with b, if the result does not exceed
  KLCOEFF_MAX; otherwise it sets the error KLCOEFF_OVERFLOW and leaves
  a unchanged.
*/

{
  if (b <= klsupport::KLCOEFF_MAX - a)
    a += b;
  else
    error::ERRNO = error::KLCOEFF_OVERFLOW;

  return a;
}

SKLcoeff& safeAdd(SKLcoeff& a, const SKLcoeff& b)

/*
  This function increments a with b if the result lies in the interval
  [SKLCOEFF_MIN,SKLCOEFF_MAX]; sets the error SKLCOEFF_OVERFLOW if we
  exceed SKLCOEFF_MAX, and SKLCOEFF_UNDERFLOW if we are less than
  SKLCOEFF_MIN.

  Note that overflow can occur only if b is positive, underflow only
  if b is negative.
*/

{
  if ((b > 0) && (a > klsupport::SKLCOEFF_MAX - b)) {
    error::ERRNO = error::SKLCOEFF_OVERFLOW;
  }
  else if ((b < 0) && (a < klsupport::SKLCOEFF_MIN - b)) {
    error::ERRNO = error::SKLCOEFF_UNDERFLOW;
  }
  else
    a += b;

  return a;
}

KLCoeff& safeMultiply(KLCoeff& a, const KLCoeff& b)

/*
  This function multiplies a with b, if the result does not exceed
  KLCOEFF_MAX; otherwise it sets the error KLCOEFF_OVERFLOW and leaves
  a unchanged.
*/

{
  if (a == 0)
    return a;

  if (b <= klsupport::KLCOEFF_MAX/a)
    a *= b;
  else
    error::ERRNO = error::KLCOEFF_OVERFLOW;

  return a;
}

SKLcoeff& safeMultiply(SKLcoeff& a, const SKLcoeff& b)

/*
  This function multiplies a with b, if the result lies between SKLCOEFF_MIN
  and SKLCOEFF_MAX. Otherwise it sets the error SKLCOEFF_UNDERFLOW or
  SKLCOEFF_OVERFLOW as appropriate.
*/

{
  if (a == 0)
    return a;

  if (a > 0) {
    if (b > klsupport::SKLCOEFF_MAX/a)
      error::ERRNO = error::SKLCOEFF_OVERFLOW;
    else if (b < klsupport::SKLCOEFF_MIN/a)
      error::ERRNO = error::SKLCOEFF_UNDERFLOW;
    else
      a *= b;
  }
  else {
    if (b > klsupport::SKLCOEFF_MIN/a)
      error::ERRNO = error::SKLCOEFF_UNDERFLOW;
    else if (b < klsupport::SKLCOEFF_MAX/a)
      error::ERRNO = error::SKLCOEFF_OVERFLOW;
    else
      a *= b;
  }

  return a;
}

KLCoeff& safeSubtract(KLCoeff& a, const KLCoeff& b)

/*
  This function subtracts b from a, if the result is non-negative; sets
  the error KLCOEFF_UNDERFLOW otherwise, and leaves a unchanged.
*/

{
  if (b <= a)
    a -= b;
  else
    error::ERRNO = error::KLCOEFF_UNDERFLOW;

  return a;
}

};

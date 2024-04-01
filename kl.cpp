/*
  This is kl.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "kl.h"

#include <cassert>

#include "sl_list.h"
#include "error.h"
#include "iterator.h"

namespace kl {
  using namespace error;


/****************************************************************************

  This file contains the code for the computation of the (ordinary) Kazhdan-
  Lusztig polynomials for W.

  Our storage class |KLContext| holds two main tables (actually, like all its
  data, inside its helper class |KLHelper|), |KL_table| for the Kazhdan- Lusztig
  polynomials $P(x,y)$ themselves, and |mu_Table| for the nonnegative integers
  $\mu(x,y)$ determined by them; in fact in both cases only a necessary part of
  all thes values are actually stored. The basic range is all pairs $(x,y)$ of
  Coceter group elements currently covered in the |schubert()| context for which
  $x\leq y$ iin the Bruhat order. Each table is a |std::vector| indexed by |y|
  of pointers to "rows" indexed by |x| (so they are actually columns in a sparse
  matrix); their type are |KLRow| and a |MuRow| respectively). Resizing function
  are provided so the new space can be added when the Schubert context grows.

  The row |KL_table[y]|, when non null, contains one entry for each extremal
  $x\leq y$; the list of those |x| is available as |klsupport().extrList|, so
  that these lists are shared among the different |KLContext| tbales. Each entry
  in |*KL_list[y]| is a non-owned pointer to a (constant) polynomial $P(x,y)$.
  The main complication is that I've (rightly or wrongly) decided to allocate
  |KL_table[y]| only if $y\leq inverse(y)$ (numerically); if not, the
  polynomials $P(x,y)$ are looked up as $P(x^{-1},y^{-1}$; this saves space, but
  makes lookup more complicated. More importantly though, systematically going
  over to |inverse(y)| when it is smaller thatn |y|, seems to make the algorithm
  quite a bit faster and leaner (in the sense that fewer rows need to be
  computed); this, more than the memory saving, is my main reason for keeping
  this complication in.

  The situation for the |mu_Table| is a little bit more delicate. First of all,
  we do not here go over to inverses; we fill all the rows that are needed
  (using the mu-value for the inverses if it is available.) Further, we only
  look at pairs where the length difference is at least three (it not the $\mu$
  vale is implied by the Bruhat comparison); it is known that in those cases the
  nonzero $\mu$ only occur at extremal pairs. But in fact, even among extremal
  pairs, most values are zero as well. So we do not necessarily allocate one
  entry for each pair; what is guaranteed is that all entries correspond to
  extremal pairs, and all non-zero mu-values have an entry. Ideally, we would
  wish to allocate only the non-zero mu's; however, this would require computing
  a full row as soon as one coefficient is needed, so we have refrained from
  that. [What Fokko appears to imply here is that if were to maintain only
  nonzero $\mu$ values, we would either have to compute a whole row at once to
  decide which ones are nonzero, and then we can subsequently deduce that absent
  values must be zero, or else insert nonzero values as we go, but this would
  require growing the row continuously, and also would require recomputing the
  values that are not found even if they turn out to be zero, because in this
  scenario we cannot infer zeroness from absence. MvL] Anyway, to find out the
  value for $\mu(x,y)$, we extremalize $x$, check that the length difference is
  odd and at least 3, if so look x up in |*mu_Table[y]|, and return the
  corresponding mu-value if a slot with a well defined value (which could be 0)
  is found, compute the $\mu$ value if a slot is found with an undefined value,
  and 0 if no slot is found at all (namely for non extrmal |x|).

  The idea is to compute everything upon request: we compute exactly what
  is needed in the recursive computation of what is required.

  The requests that we mainly focus on are the following :

    - compute a single K-L polynomial (or mu-coefficient);
    - compute all P_{x,y}'s (or mu(x,y)'s) for a given y;
    - compute the full table of K-L polynomials (or mu coefficients);

  It turns out that the dominant factor in the computation is the Bruhat order
  comparison. That is why computing a whole row can be done much more
  efficently than computing each entry in the row separately: using our
  closure function from the schubert context, we can basically factor out
  most of the Bruhat things for a whole row. In any case, the row computation
  can be done without a single call to the slow ruhat |inOrder| function!

 ****************************************************************************/

struct KLContext::KLHelper
{
// data
  klsupport::KLSupport& d_klsupport; // unowned, |CoxGroup| owns it

  KLTable KL_table;
  MuTable mu_Table;
  containers::bag<KLPol> KL_pool;
  KLStats d_stats;

// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
      {return memory::arena().free(ptr,sizeof(KLHelper));}
  KLHelper(klsupport::KLSupport& kls,KLContext* kl);

// methods
// relay methods, mostly |const|:
  const klsupport::KLSupport& klsupport() const // method used when |const|
    { return d_klsupport; }
  const KLStats& stats() const { return d_stats; } // method used when |const|
  coxtypes::Rank rank() const { return klsupport().rank();}
  const schubert::SchubertContext& schubert() const
    { return klsupport().schubert(); }
  coxtypes::Generator last (const coxtypes::CoxNbr& x) const
    { return klsupport().last(x); }
  const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y) const
    { return klsupport().extrList(y); }
  coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& y) const
    { return klsupport().inverse(y); }


  Ulong size() const { return KL_table.size(); }

  KLRow& KL_row(coxtypes::CoxNbr y) { return *KL_table[y]; }
  MuRow& mu_row(coxtypes::CoxNbr y) { return *mu_Table[y]; }

  bool row_needs_creation(coxtypes::CoxNbr x) const
    { return KL_table[x] == nullptr; }
  void create_KL_row(coxtypes::CoxNbr y);
  bool row_is_incomplete(coxtypes::CoxNbr y);
  void take_mu_row_from_inverse(const coxtypes::CoxNbr& y);
  void ensure_correction_terms
    (coxtypes::CoxNbr y, coxtypes::Generator s, bitmap::BitMap* done);
  containers::vector<KLPol> initial_polys
    (coxtypes::CoxNbr y,coxtypes::Generator s);
  void add_second_terms(containers::vector<KLPol>& pols, coxtypes::CoxNbr y);
  void mu_correct_row(containers::vector<KLPol>& pols, coxtypes::CoxNbr y);
  void coatom_correct_row(containers::vector<KLPol>& pols,coxtypes::CoxNbr y);
  void compute_KL_row(coxtypes::CoxNbr y, bitmap::BitMap* done = nullptr);
  void copy_mu_row_from_KL(coxtypes::CoxNbr y);
  void fill_KL_table ();

  const KLPol* compute_KL_pol
    (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s);

  void allocMuRow(const coxtypes::CoxNbr& y);
  void allocMuTable();

  bool mu_row_is_complete(coxtypes::CoxNbr y);
    klsupport::KLCoeff computeMu
      (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
    void fillMuRow(MuRow& row, const coxtypes::CoxNbr& y);
  void writeKLRow(const coxtypes::CoxNbr& y, containers::vector<KLPol>& pols);
    bool isExtrAllocated(const coxtypes::CoxNbr& y)
      { return klsupport().isExtrAllocated(y); }
  bool isMuAllocated(const coxtypes::CoxNbr& x) const
  { return mu_Table[x] != nullptr; }
    void makeMuRow(const coxtypes::CoxNbr& y);
  void mu_correct_KL_pol
    (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
     const coxtypes::Generator& s, KLPol& pol);
    klsupport::KLCoeff recursiveMu
      (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
       const coxtypes::Generator& s);
    void writeMuRow(const MuRow& row, const coxtypes::CoxNbr& y);

  void coatom_correct_KL_pol
    (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
     const coxtypes::Generator& s,  KLPol& pol);

  void grow(Ulong prev, Ulong n);
  void shrink(const Ulong& n);

  void fill_mu_table ();

  const KLPol& klPol(coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  const KLPol& klPol(coxtypes::CoxNbr x, coxtypes::CoxNbr y,
		     coxtypes::Generator s);
  klsupport::KLCoeff mu(coxtypes::CoxNbr x, coxtypes::CoxNbr y);

  bool KL_is_filled() const { return stats().flags&KLStats::kl_done; }
  bool mu_is_filled() const { return stats().flags&KLStats::mu_done; }
  void set_KL_filled() { d_stats.flags |= KLStats::kl_done;}
  void set_mu_filled() { d_stats.flags |= KLStats::mu_done;}
  void clear_KL_filled() { d_stats.flags &= ~KLStats::kl_done;}
  void clear_mu_filled() { d_stats.flags &= ~KLStats::mu_done;}

  HeckeElt KL_row_as_HeckeElt(coxtypes::CoxNbr y);
  void move_KL_row_from_inverse(coxtypes::CoxNbr x);
  void permute(const bits::Permutation& a);
}; // |struct KLContext::KLHelper|

namespace {
  using namespace kl;

  struct PolPtrF {
    typedef const KLPol* valueType;
    const KLPol* operator()(const KLPol* pol) {return pol;}
  };

  MuData* find(MuRow& row, const coxtypes::CoxNbr& x);
  KLPol& safeAdd(KLPol& p, const KLPol& q, const polynomials::Degree& n);
  KLPol& safeSubtract(KLPol& p, const KLPol& q, const klsupport::KLCoeff& mu,
		      const coxtypes::Length& h);
  void showSimpleMu(FILE* file, KLContext& kl,
		    coxtypes::CoxNbr x, coxtypes::CoxNbr y,
		    klsupport::KLCoeff r, const interface::Interface& I);
  void showRecursiveMu(FILE* file, KLContext& kl,
		       coxtypes::CoxNbr x, coxtypes::CoxNbr y,
		       klsupport::KLCoeff r, const interface::Interface& I);
  KLPol& zeroPol();
}; // |namespace|

/****************************************************************************

        Chapter I --- The KLContext class

  This section defines the functions for the KLContext class. The following
  functions are defined :

    constructors and destructors :

      - KLContext(KLSupport* kls);
      - ~KLContext();

    accessors :

    manipulators :

      - applyInverse(y) : exchanges rows in klList for y and y_inverse;
      - fillKL() : fills the full K-L table;
      - fillMu() : fills the full mu-table;
      - klPol(x,y) : for x <= y, returns the Kazhdan-Lusztig polynomial;
      - klPol(x,y,s) : same as above, using s as descent;
      - lcell() : returns the left cell partition;
      - lrcell() : returns the two-sided cell partition;
      - mu(x,y) : for x <= y, returns the mu-coefficient;
      - rcell() : returns the right cell partition;
      - row(h,y,I) : returns the row of the coxtypes::CoxNbr y in h;
      - row(e_row,kl_row,y) : returns the row of the coxtypes::CoxNbr y;
      - revertSize(n) : reverts the size to size n;
      - setSize(n) : extends the context to size n;

    input/output :

      - printStatus(file) : prints the status;

 ****************************************************************************/

KLContext::KLContext(klsupport::KLSupport* kls)
  : d_help(new KLHelper(*kls,this))
{
}


/******** accessors **********************************************************/

const klsupport::KLSupport& KLContext::klsupport() const
  { return d_help->d_klsupport; }
klsupport::KLSupport& KLContext::klsupport() { return d_help->d_klsupport; }
const KLStats& KLContext::stats() const { return d_help->stats(); }
Ulong KLContext::size() const { return d_help->size(); }
const KLRow& KLContext::klList(const coxtypes::CoxNbr& y) const
  { return *d_help->KL_table[y]; }
const MuRow& KLContext::muList(const coxtypes::CoxNbr& y) const
  { return *d_help->mu_Table[y]; }

/******** manipulators *******************************************************/


/*
  This function extends the context so that it can accomodate |n| elements.
  The idea is to have the same size as the basic schubert context.
*/
void KLContext::setSize(const Ulong& n)
{
  d_help->grow(size(),n);
}

void KLContext::KLHelper::grow(Ulong prev, Ulong n)
{
  coxtypes::CoxNbr prev_size = size();

  CATCH_MEMORY_OVERFLOW = true;

  KL_table.resize(n);
  if (ERRNO)
    goto revert;
  mu_Table.resize(n);
  if (ERRNO)
    goto revert;

  CATCH_MEMORY_OVERFLOW = false;

  clear_KL_filled();
  clear_mu_filled();

  return;

 revert:
  CATCH_MEMORY_OVERFLOW = false;
  shrink(prev_size);
} // |KLContext::KLHelper::grow|




// This function fills all the rows in klList, in a straightforward way.
void KLContext::fillKL() { d_help->fill_KL_table(); }

/*
  This function fills all the rows in the mu-list, in a straightforward way.
  Sets the error ERROR_WARNING in case of error.

  NOTE : error handling should be improved!
*/
void KLContext::fillMu() { d_help->fill_mu_table(); }

void KLContext::KLHelper::fill_mu_table ()
{
  if (mu_is_filled())
    return;

  allocMuTable();

  if (ERRNO)
    goto abort;

  for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) < y)
      take_mu_row_from_inverse(y);
    fillMuRow(*mu_Table[y],y);
    if (ERRNO)
      goto abort;
  }

  set_mu_filled();
  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}


/*
  This function returns the Kazhdan-Lusztig polynomial P_{x,y}. It is
  assumed that the condition x <= y has already been checked, and that
  x and y are valid context numbers.
*/
const KLPol& KLContext::klPol (coxtypes::CoxNbr x, const coxtypes::CoxNbr y)
{ return d_help->klPol(x,y); }


const KLPol& KLContext::klPol
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s)
{ if (s==coxtypes::undef_generator) // this used to be a default, no longer is
    return d_help->klPol(x,y);
  return d_help->klPol(x,y,s);
}

const KLPol& KLContext::KLHelper::klPol
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s)
{ // we modify |x| and |y| as local variables, but $P(x,y)$ remains the same

  assert(s != coxtypes::undef_generator);
  // we are asked for an explicit right descent |s|, so we cannot use inverse

  const schubert::SchubertContext& p = schubert();

  // make |x| extramal for |y|
  x = p.maximize(x,p.descent(y));

  // check for (now) trivially small cases
  if (p.length(y) - p.length(x) < 3) { /* result is 1 */
    return one();
  }

  // make sure |KL_table[y]| is allocated (and hence |extrList(y)| as well)
  if (row_needs_creation(y)) {
    create_KL_row(y);
    if (ERRNO)
      return zeroPol();
  }

  // find |x| in |extrList[y]|
  const auto& eL = extrList(y);
  Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();
  const KLPol*& pol = KL_row(y)[m];

  if (pol == nullptr) { /* we have to compute the polynomial */
    pol = compute_KL_pol(x,y,s);
    if (ERRNO)
      return zeroPol();
  }

  return *pol;
} // |KLContext::KLHelper::klPol(x,y,s)|



/*
  This function returns the mu-coefficient mu(x,y). It is assumed that
  the condition x <= y has already been checked, and that x and y are
  valid context numbers.

  The return value is zero if the length difference is even.

  If an error occurs, it forwards the error value and returns the
  value undef_klcoeff for mu.
*/
klsupport::KLCoeff KLContext::mu
  (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)
{ return d_help->mu(x,y); }

klsupport::KLCoeff KLContext::KLHelper::mu
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{
  const schubert::SchubertContext& p = schubert();

  coxtypes::Length d = p.length(y) - p.length(x);

  if (d%2 == 0)
    return 0;

  if (d == 1) // x is a coatom of y
    return 1;

  // check if x is in extremal position w.r.t. y
  if (x != p.maximize(x,p.descent(y)))
    return 0;

  // allocate |*mu_Table[y]| if necessary
  if (not isMuAllocated(y)) {
    allocMuRow(y);
    if (ERRNO)
      return klsupport::undef_klcoeff;
  }

  // find x in |*mu_Table[y]|
  MuRow& m = *mu_Table[y];
  MuData* md = find(m,x);
  if (md == nullptr)
    return 0;

  if (md->mu == klsupport::undef_klcoeff) { // we need to compute $\mu(x,y)$
    md->mu = computeMu(x,y);
    if (ERRNO)
      return klsupport::undef_klcoeff;
  }

  return md->mu;
} // |KLContext::KLHelper::mu(x,y)|

/*
  This function returns in h the data for the full row of y in the K-L table,
  sorted in the context number order.
*/
void KLContext::row(HeckeElt& h, const coxtypes::CoxNbr& y)
{ h = d_help->KL_row_as_HeckeElt(y); }

HeckeElt KLContext::KLHelper::KL_row_as_HeckeElt(coxtypes::CoxNbr y)
{
  HeckeElt h;
  if (row_is_incomplete(y))
  { bitmap::BitMap done(y+1); // start with no completion information
    compute_KL_row(y,&done);
  }

  if (ERRNO) {
    Error(ERRNO);
    ERRNO = ERROR_WARNING;
    return h;
  }

  if (y <= inverse(y)) {
    const klsupport::ExtrRow& e = extrList(y);
    h.reserve(e.size());
    const KLRow& klr = KL_row(y);
    for (Ulong j = 0; j < e.size(); ++j)
      h.emplace_back(e[j],klr[j]);
  }
  else { /* go over to inverses */
    coxtypes::CoxNbr yi = inverse(y);
    const klsupport::ExtrRow& e = extrList(yi);
    h.reserve(e.size());
    const KLRow& klr = KL_row(yi);
    for (Ulong j = 0; j < e.size(); ++j)
      h.emplace_back(inverse(e[j]),klr[j]);

    std::sort(h.begin(),h.end()); // make sure list is ordered
  }
  return h;
} // |::KLHelper::KL_row_as_HeckeElt|


/*
  Move row from |inverse(x)| to |x| in |KL_table|, assuming that the boths rows
  are within the bounds of |KL_table|. The row is unchanged (!)
*/
void KLContext::applyInverse(const coxtypes::CoxNbr& x)
{ d_help->move_KL_row_from_inverse(x);
}

// this function moes in the direction of what its name suggests!
void KLContext::KLHelper::move_KL_row_from_inverse(coxtypes::CoxNbr x)
{
  coxtypes::CoxNbr xi = inverse(x);
  KL_table[x] = std::move(KL_table[xi]);
}

void KLContext::applyIPermutation
  (const coxtypes::CoxNbr& y, const bits::Permutation& a)
{ return right_permute(*d_help->KL_table[y],a); }

/*
  This function permutes the context according to the permutation a. The
  following objects have to be permuted :

   - klList, muList : each row is a table with values in the context; the
     list itself has range in the context; don't forget to sort the permuted
     rows!
   - inverse : range and values in the context, with some annoyance caused
     by undef_coxnbr values;
   - last : range in the context;

  Applying a permutation to something goes as follows. The permutation
  a is the table of the function which assigns to each x in the context its
  number in the new enumeration. Applying this to the context means the
  following :

    (a) for tables with _values_ in the context, we need to replace t[j] with
        a[t[j]];
    (b) for tables with _range_ in the context, we need to put t[x] at a[x];
        this amounts to replacing t[x] by t[a_inv[x]] for all x;
    (c) for tables with both range and values in the context, we need to do
        both : replace t[x] with a[t[a_inv[x]]]; in practice this is done
        by first replacing the values as in (a), then permuting as in (b).

  The replacements in (a) are straightforward. The permutations in (b) are
  more delicate, and would seem to require additional memory for each list to
  be transformed. In fact, things can be done in place. The thing is to use
  the cycle decomposition of the permutation : put t[x] at position a[x],
  t[a[x]] at position a^2[x], ... t[a^{k-1}[x]] at position x, if a^k[x] = x.
  This requires only the use of the bitmap of a to mark off the entries that
  have been handled, and skip to the next entry.

*/
void KLContext::permute(const bits::Permutation& a)
{ d_help->permute(a); }

void KLContext::KLHelper::permute(const bits::Permutation& a)
{
  // permute values inside each row of |muTable|
  for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
    if (!isMuAllocated(y))
      continue;
    MuRow& row = *mu_Table[y];
    for (Ulong j = 0; j < row.size(); ++j)
      row[j].x = a[row[j].x];
    std::sort(row.begin(),row.end());
  }

  // permute ranges
  bitmap::BitMap seen(a.size());

  for (coxtypes::CoxNbr x = 0; x < size(); ++x) {
    if (seen.is_member(x))
      continue;
    if (a[x] == x) {
      seen.insert(x);
      continue;
    }

    for (coxtypes::CoxNbr y = a[x]; y != x; y = a[y])
    {
      /* back up values for y */
      auto kl_buf = std::move(KL_table[y]);
      auto mu_buf = std::move(mu_Table[y]);
      /* put values for x in y */
      KL_table[y] = std::move(KL_table[x]);
      mu_Table[y] = std::move(mu_Table[x]);
      /* store backup values in x */
      KL_table[x] = std::move(kl_buf);
      mu_Table[x] = std::move(mu_buf);
      /* set bit*/
      seen.insert(y);
    }  // |for(y)|

    seen.insert(x);
  } // |for(x)|
} // |KLContext::KLHelper::permute|

/*
  Revert the sizes of the lists to size |n|. This is meant to be used
  only immediately after a failing context extension, to preserve the
  consistency of the various list sizes. In particular, it will fail
  miserably if a premutation has taken place in-between.
*/
void KLContext::revertSize(const Ulong& n)
{
  d_help->shrink(n);
}

void KLContext::KLHelper::shrink(const Ulong& n)
{
  KL_table.resize(n);
  mu_Table.resize(n);
} // |shrink|


/****************************************************************************

        Chapter II -- The KLHelper class

  The purpose of the KLHelper class is to hide from the public eye a number
  of helper functions, used in the construction and maintenance of the
  K-L context. This unclutters kl.h quite a bit.

  The following functions are defined :

   - create_KL_row(const coxtypes::CoxNbr& y) : allocate row in the K-L list;
   - initial_polys(coxtypes::CoxNbr y, coxtypes::Generator s):
     another preliminary to the computation of a row;
   - add_second_terms(containers::vector<KLPol>& pols, coxtypes::CoxNbr y) :
     incorporates the second terms P_{x,ys} in the computation of a full row;
   - mu_correct_row(containers::vector<KLPol>& pols,coxtypes::CoxNbr y):
     subtracts the non-coatom mu-part in the computation of a row;
   - row_is_incomplete(const coxtypes::CoxNbr& y) : whether a K-L row is not
     yet fully computed;
   - allocMuRow(const coxtypes::CoxNbr& y) : allocate row in the mu list;
   - allocMuTable() : allocates the full mu-table;
   - coatom_correct_row( containers::vector<KLPol>& pol,coxtypes::CoxNbr y):
     subtractd the terms for coatoms in the mu-correction, for a full row;
   - coatom_correct_KL_pol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s,
     containers::vector<KLPol>& pol, const Ulong& a) : same, for a single polynomial
   - computeMu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y) : computes a mu-coefficient;
   - fillMuRow(MuRow& row, const coxtypes::CoxNbr& y) : fills a row in the mu-table;
   - compute_KL_pol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s =
     coxtypes::undef_generator) : fills in one polynomial, using s as descent;
   - compute_KL_row(coxtypes::CoxNbr y) : fills in one row in the K-L table;
     coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& y) : returns the inverse of y;
   - take_mu_row_from_inverse(const coxtypes::CoxNbr& y) :
     construct the mu-row for y from that of the inverse of y;
   - makeMuRow(const coxtypes::CoxNbr& y);
   - mu_correct_KL_pol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s,
     containers::vector<KLPol>& pol, const Ulong& a) : subtracts the non-coatom mu-part,
     for the computation of a single polynomial;
   - mu_row(const coxtypes::CoxNbr& y) : returns the row for y in muList;
   - ensure_correction_terms(coxtypes::CoxNbr y, coxtypes::Generator s,
           BitMap* done) : a preliminary to the computation of a row;
   - copy_mu_row_from_KL(const coxtypes::CoxNbr& y) :
     fills in the mu-row from the K-L row;
   - recursiveMu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s) :
     computes mu(x,y) using the general recursive formula;
   - writeKLRow(const coxtypes::CoxNbr& y, containers::vector<KLPol>& pol) : transfers the
     polynomials from |pol| to |KL_row(y)|;
   - writeMuRow(const MuRow& row, const coxtypes::CoxNbr& y) transfers the
     mu-coefficients from row to muList;

 ****************************************************************************/

KLContext::KLHelper::KLHelper(klsupport::KLSupport& kls,KLContext* kl)
  : d_klsupport(kls)
  , KL_table(kls.size())
  , mu_Table(kls.size())
  , KL_pool()
  , d_stats()
{
  KL_table[0].reset(new KLRow(1));
  (*KL_table[0])[0] = KL_pool.find(one());

  mu_Table[0].reset(new MuRow);
  d_stats.klnodes++;
  d_stats.klrows++;
  d_stats.klcomputed++;
  d_stats.murows++;
}


/*
  Allocate one row of the kl_list. The row contains one
  entry for each x <= y which is extremal w.r.t. the descent set of y.
*/
void KLContext::KLHelper::create_KL_row(coxtypes::CoxNbr y)
{
  Ulong n = d_klsupport.extr_list(y).size(); // might generate that lsit

  KL_table[y].reset(new KLRow(n));
  d_stats.klnodes += n;
  d_stats.klrows++;
}

// Whether the row for |y| (or for |inverse(y)|) in |klList| has been filled.
bool KLContext::KLHelper::row_is_incomplete(coxtypes::CoxNbr y)
{
  if (inverse(y) < y)
    y = inverse(y);

  if (row_needs_creation(y)) /* row is not allocated */
    return true;

  KLRow& kl_row = KL_row(y);

  for (const auto pol :  kl_row)
    if (pol == nullptr)
      return true;

  return false;
}



/*
  Allocate one row in the muList. There is one entry for each x < y which is
  extremal w.r.t. y, and has odd length-difference > 1 with y. As with
  |create_KL_row|, this function is not designed for maximal efficiency; row
  allocations for big computations should be handled differently.
*/
void KLContext::KLHelper::allocMuRow(const coxtypes::CoxNbr& y)
{
  using EI = iterator::FilteredIterator
    <coxtypes::CoxNbr,klsupport::ExtrRow::const_iterator,MuFilter>;
  using BI = iterator::FilteredIterator
    <Ulong,bitmap::BitMap::const_iterator,MuFilter>;

  const schubert::SchubertContext& p = schubert();
  klsupport::ExtrRow e;
  MuFilter f(p,y);

  if (isExtrAllocated(y)) {
    EI first(extrList(y).begin(),extrList(y).end(),f);
    EI last(extrList(y).end(),extrList(y).end(),f);
    e.assign(first,last);
  }
  else
  {
    auto b = p.closure(y);
    if (ERRNO)
      return;
    schubert::select_maxima_for(p,p.descent(y),b);
    BI first(b.begin(),b.end(),f);
    BI last(b.end(),b.end(),f);
    e.assign(first,last);
  }

  coxtypes::Length ly = p.length(y);

  mu_Table[y].reset(new MuRow);
  if (ERRNO)
  {
    Error(ERRNO);
    ERRNO = ERROR_WARNING;
    return;
  }

  auto& dest = mu_row(y);
  dest.reserve(e.size());

  for (Ulong j = 0; j < e.size(); ++j)
  {
    coxtypes::CoxNbr x = e[j];
    coxtypes::Length lx = p.length(x);
    dest.emplace_back(x,klsupport::undef_klcoeff,(ly-lx-1)/2);
  }

  d_stats.munodes += e.size();
  d_stats.murows++;
}


/*
  Allocate the full muList, using the closure iterator. It is not as
  satisfactory as it should be, because the rows are not obtained in order, and
  therefore we can not fill them right away to eliminate zero entries. So I am
  currently not too happy with the situation. Even if we recover the memory
  later on, it will be highly fragmented.

  So currently this is just an exercise in formal elegance.

*/
void KLContext::KLHelper::allocMuTable()
{
  const schubert::SchubertContext& p = schubert();
  klsupport::ExtrRow buffer;

  for (schubert::ClosureIterator cit(p); cit; ++cit)
  {
    coxtypes::CoxNbr y = cit.current();
    if (inverse(y) < y)
      continue;
    if (isMuAllocated(y))
      continue;

    const coxtypes::Length ly = p.length(y);
    /* find extremal list */
    bitmap::BitMap b = cit.closure();
    if (ERRNO) {
      printf("error! y = %lu\n",static_cast<Ulong>(y));
      Error(ERRNO);
      ERRNO = ERROR_WARNING;
      return;
    }

    schubert::select_maxima_for(p,p.descent(y),b);

    buffer.clear();
    for (coxtypes::CoxNbr x : b)
    { coxtypes::Length dl = ly-p.length(x);
      if (dl%2!=0 and dl!=1)
	buffer.push_back(x);
    }

    /* transfer to muList */

    mu_Table[y].reset(new MuRow);
    auto& dest = mu_row(y);
    dest.reserve(buffer.size());

    for (coxtypes::CoxNbr x : buffer)
      dest.emplace_back(x,klsupport::undef_klcoeff,(ly-p.length(x))/2); // round down

    d_stats.murows++;
    d_stats.munodes += dest.size();
  } // |for(cit)|
} // |KLHelper::acclocMuTable|



// Whether the row for y in muList has been filled.
bool KLContext::KLHelper::mu_row_is_complete(coxtypes::CoxNbr y)
{
  if (!isMuAllocated(y)) /* row is not allocated */
    return false;

  const MuRow& row = mu_row(y);

  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu == klsupport::undef_klcoeff)
      return false;
  }

  return true;
}


/*
  This function subtracts the coatom correction from pol, which at this
  point should contain the value P_{xs,ys}+q.P_{x,ys}-(mu-correction)
  (although we can apply CoatomCorrection and MuCorrection in any order.)

  This means that we subtract the correcting terms corresponding to the
  coatoms z of ys s.t. zs < s; the corresponding subtraction is q.P_{x,z}.

*/
void KLContext::KLHelper::coatom_correct_KL_pol
  (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
   const coxtypes::Generator& s, KLPol& pol)
{
  const schubert::SchubertContext& p = schubert();
  coxtypes::CoxNbr ys = p.shift(y,s);
  const schubert::CoxNbrList& c = p.hasse(ys);

  for (Ulong j = 0; j < c.size(); ++j) {

    coxtypes::CoxNbr z = c[j];

    if (p.shift(z,s) > z) /* z is not considered */
      continue;

    if (not p.Bruhat_leq(x,z)) /* z is not in [x,ys] */
      continue;

    /* at this point we have to do an actual subtraction */

    const KLPol& p_xz = klPol(x,z);
    if (ERRNO)
      return;
    safeSubtract(pol,p_xz,1,1);
    if (ERRNO) {
      Error(ERRNO,this,x,y);
      ERRNO = ERROR_WARNING;
      return;
    }
  }
}


/*
  This function gets a previously uncomputed entry in the muList. It
  is based on the following remark (essentially Lusztig's "star-operation"
  situation). We already assume that LR(x) contains LR(y). Now let s be
  in LR(y), such that LR(ys) is _not_ contained in LR(x) (such an s exists
  iff LR(x) does not contain twoDescent(y).) Then let t be in LR(ys), not in
  LR(x) (and hence not in LR(y)). Then it must be so that s and t do not
  commute; in particular they act on the same side; assume this is on the
  right. So we have yst < ys < y < yt, xs < x < xt. Assume that x <= ys
  (otherwise mu(x,y) = mu(xs,ys)). Then from the fact that xt > x, yst < ys,
  exactly as in the proof of thm. 4.2. in the original K-L paper, one sees
  that at most four terms survive in the recursion formula : we have

  mu(x,y) = mu(xs,ys) + mu(xt,ys) - mu(x,yst)(if ysts < yst)
            - mu(xt,ys)(if xts < xt)

  So in all these cases we get an elementary recursion.

  Sets the error MU_FAIL, and returns the value undef_klcoeff, in case of
  failure (this can be due to memory overflow, or to coefficient over- or
  underflow.)
*/
klsupport::KLCoeff KLContext::KLHelper::computeMu
  (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y)
{
  if (inverse(y) < y)
    return computeMu(inverse(x),inverse(y));

  const schubert::SchubertContext& p = schubert();

  Lflags f = p.twoDescent(y);

  if ((p.descent(x)&f) == f) { /* x is super-extremal w.r.t. y */
    return recursiveMu(x,y,last(y));
  }

  coxtypes::Generator s=-1, t=-1; // init silences compiler warnings

  /* choose s s.t. LR(ys) not contained in LR(x) */
  Lflags desc;
  for (desc = p.descent(y); desc; desc &= desc-1)
  {
    coxtypes::Generator u = constants::firstBit(desc);
    coxtypes::CoxNbr yu = p.shift(y,u);
    Lflags fu = p.descent(yu);
    if ((p.descent(x)&fu) != fu) {
      s = u;
      t = constants::firstBit(fu & ~p.descent(x));
      break;
    }
  }
  assert(desc!=0);

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);

  klsupport::KLCoeff r1 = mu(xs,ys);
  if (ERRNO)
    goto abort;

  if (not p.Bruhat_leq(x,ys)) // check whether x <= ys
  { // if not, value $\mu=0$ is found
    d_stats.mucomputed++;
    if (r1 == 0)
      d_stats.muzero++;
    return r1;
  }

  {
    coxtypes::CoxNbr xt = p.shift(x,t);
    coxtypes::CoxNbr yst = p.shift(ys,t);

    /* consider four cases */

    if (p.isDescent(xt,s)) { /* xts < xt */
      if (p.isDescent(yst,s)) { /* ysts < yst */
	klsupport::KLCoeff r3 = mu(x,yst);
	if (ERRNO)
	  goto abort;
	if (r1 < r3) { /* negative mu-coefficient */
	  ERRNO = MU_NEGATIVE;
	  goto abort;
	}
	d_stats.mucomputed++;
	if (r1 == r3)
	  d_stats.muzero++;
	return r1-r3;
      }
      else { /* ysts > yst */
	d_stats.mucomputed++;
	if (r1 == 0)
	  d_stats.muzero++;
	return r1;
      }
    }
    else { /* xts > xt */
      if (p.isDescent(yst,s)) { /* ysts < yst */
	klsupport::KLCoeff r2 = mu(xt,ys);
	if (ERRNO)
	  goto abort;
	klsupport::KLCoeff r3 = mu(x,yst);
	if (ERRNO)
	  goto abort;
	if ((r1+r2) < r3) { /* negative mu-coefficient */
	  ERRNO = MU_NEGATIVE;
	  goto abort;
	}
	d_stats.mucomputed++;
	if (r1+r2 == r3)
	  d_stats.muzero++;
	return r1+r2-r3;
      }
      else { /* ysts > yst */
	klsupport::KLCoeff r2 = mu(xt,ys);
	if (ERRNO)
	  goto abort;
	d_stats.mucomputed++;
	if (r1+r2 == 0)
	  d_stats.muzero++;
	return r1+r2;
      }
    }
  }

 abort:
  if (ERRNO != MEMORY_WARNING)
    ERRNO = MU_FAIL;
  return klsupport::undef_klcoeff;
} // |KLContext::KLHelper::computeMu|


// ----------------- K-L computation for an entire row  ----------------------



/*
  This function fills the mu-row from the corresponding KL-row. If the
  row has not been allocated yet, it allocates for nonzero values only
*/
void KLContext::KLHelper::copy_mu_row_from_KL(coxtypes::CoxNbr y)
{
  assert(not row_is_incomplete(y));
  const KLRow& kl_row = KL_row(y); // so this row is complete
  const klsupport::ExtrRow& e_row = extrList(y); // these are its indices

  if (isMuAllocated(y))
  { // a row already exists, write its |mu| fields with coefs from |kl_row|
    Ulong i = 0;
    for (MuData& mu : mu_row(y))
    {
      coxtypes::CoxNbr x = mu.x;
      while(e_row[i] < x) // search for corresponding |i| in |extrList(y)|
	++i;
      assert(e_row[i]==x); // nonzero |mu| implies extremal

      const KLPol& pol = *kl_row[i];
      coxtypes::Length d = mu.height;
      if (pol.deg() == d)
	mu.mu = pol[d];
      else {
	mu.mu = 0;
	d_stats.muzero++;
      }
      d_stats.mucomputed++;
    } // |for(mu)|
  }
  else
  { // make row from scratch, use only those values that are nonzero
    const schubert::SchubertContext& p = schubert();
    coxtypes::Length ly = p.length(y);

    containers::sl_list<MuData> mus;

    for (Ulong j = 0; j < kl_row.size(); ++j) // |j| also indexes |e_row|
    {
      coxtypes::CoxNbr x = e_row[j];
      coxtypes::Length lx = p.length(x);
      auto dl = ly-lx;

      if (dl == 1 or dl%2 == 0)
	continue;

      polynomials::Degree d = dl/2; // namely |(dl-1)/2|, but rounding is down
      const KLPol& pol = *kl_row[j];

      if (pol.deg() < d)
	continue;

      mus.emplace_back(x,pol[d],d);
    } // |for(j)|

    mu_Table[y].reset(new MuRow(mus.to_vector()));
    if (ERRNO)
      goto abort;

    d_stats.munodes += mus.size();
    d_stats.mucomputed += mus.size();
    d_stats.murows++;
  } // |else|

  return;

 abort:
  Error(ERRNO);
  ERRNO = MEMORY_WARNING;
  return;
} // |KLContext::KLHelper::copy_mu_row_from_KL|

/*
  Construct the mu-row for |y| from that of |inverse(y)|, where it is assumed
  the two are distinct and that the row for |y| is filled in. We delete any
  partial old row for |y|, since it is faster to reconstruct all its values
  from those of |inverse(y)| than to only extract those not yet present.
*/
void KLContext::KLHelper::take_mu_row_from_inverse(const coxtypes::CoxNbr& y)
{
  coxtypes::CoxNbr yi = inverse(y);
  assert(mu_row_is_complete(yi));

  if (isMuAllocated(y))  // then update |d_stats| for coming destruction
  {
    const MuRow& m = mu_row(y);
    for (const auto& entry : m)
    {
      klsupport::KLCoeff mu = entry.mu;
      if (mu != klsupport::undef_klcoeff)
	d_stats.mucomputed--;
      if (mu == 0)
	d_stats.muzero--;
    }
    d_stats.munodes -= m.size();
  }

  mu_Table[y].reset(new MuRow(mu_row(yi))); // make a copy

  MuRow& m = mu_row(y);
  for (auto& entry : m)
    entry.x = inverse(entry.x);

  std::sort(m.begin(),m.end());

  // update status
  for (const auto& entry : m)
  {
    klsupport::KLCoeff mu = entry.mu;
    if (mu != klsupport::undef_klcoeff)
      d_stats.mucomputed++;
    if (mu == 0)
      d_stats.muzero++;
  }

  d_stats.munodes += m.size();
} // |KLContext::KLHelper::take_mu_row_from_inverse|

/*
  Prepare for the filling of row |y| in |klList| and |muList|, by making sure
  that all the correction terms are available. It is assumed that the caller
  already did |compute_KL_row(ys)|, so that every $P(x,ys)$ can be looked up.

  Deals with the error if it occurs, and sets the error ERROR_WARNING;
*/
void KLContext::KLHelper::ensure_correction_terms
  (coxtypes::CoxNbr y, coxtypes::Generator s, bitmap::BitMap* done)
{
  const schubert::SchubertContext& p = schubert();
  coxtypes::CoxNbr ys = p.rshift(y,s);

  assert(y<=inverse(y));
  assert(not row_is_incomplete(ys));

  if (not mu_row_is_complete(ys))
  {
    if (ys <= inverse(ys))
      copy_mu_row_from_KL(ys);
    else
    { auto ysi = inverse(ys);
      copy_mu_row_from_KL(ysi);
      take_mu_row_from_inverse(ys);
    }
  }

  containers::sl_list<coxtypes::CoxNbr> zs;

  // find |z| with $\mu(z,ys)>0$ and for which |s| is a descent
  for (const auto& entry : mu_row(ys))
  {
    coxtypes::CoxNbr z = entry.x;
    if (entry.mu != 0 and p.isDescent(z,s)) // only |z| for which |s| is descent
      zs.push_back(z);
  }

  // find coatoms of |ys| that for which |s| is a descent
  for (coxtypes::CoxNbr z :  p.hasse(ys))
    if (p.isDescent(z,s))
      zs.push_back(z);

  for (const auto z : zs)
    if (done==nullptr ? row_is_incomplete(z)
	: done->is_member(z) ? false
	: row_is_incomplete(z) ? true
	: (done->insert(z),false)
       )
  {
    compute_KL_row(z,done); // this enters into recursion
    if (ERRNO)
      goto abort;
  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
} // |KLContext::KLHelper::ensure_correction_terms|


/*
  Return a row of one polynomial for each |x| in |KL_row(y)|, initialized to
  |klPol(xs,ys)|.
*/
containers::vector<KLPol>
  KLContext::KLHelper::initial_polys(coxtypes::CoxNbr y, coxtypes::Generator s)
{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  containers::vector<KLPol> pols; pols.reserve(e.size());

  try {
    pols.reserve(e.size());
  }
  catch(...) {
    error::Error(error::MEMORY_WARNING);
    error::ERRNO = error::ERROR_WARNING;
    return pols;
  }

  // initialize with values $P_{xs,ys}$
  { // avoid jumping across variable definitions
    coxtypes::CoxNbr ys = p.rshift(y,s);

    for (Ulong j = 0; j < e.size(); ++j) {
      coxtypes::CoxNbr xs = p.rshift(e[j],s);
      pols.push_back(klPol(xs,ys)); // no recursion, no new memory needed here
    }
  }

  return pols;
} // |KLContext::KLHelper::initial_polys|


/*
  Take care of the "second term" q.P_{x,ys} in the recursion formula for
  P_{x,y}. It is assumed that y <= inverse(y) and that the descent strategy is
  via |last|.

  Here all the memory allocations have been made successfully; the only cause of
  error would be an overflow condition. In that case, the error is treated and
  ERROR_WARNING is set.

  Since we want to avoid all calls to InOrder, the method here is to extract
  [e,ys], extremalize it w.r.t. the descent set of y, and run through it and
  make the correction.
*/
void KLContext::KLHelper::add_second_terms
  (containers::vector<KLPol>& pols, coxtypes::CoxNbr y)
{
  const schubert::SchubertContext& p = schubert();
  coxtypes::CoxNbr ys = p.rshift(y,last(y));

  auto b = p.closure(ys); // Bruhat interval [e,ys]

  // retain two-sided extrema for |y| only
  schubert::select_maxima_for(p,p.descent(y),b);

  Ulong i = 0;
  const klsupport::ExtrRow& e = extrList(y);

  for (coxtypes::CoxNbr x : b)
  {
    while(e[i] < x) // linearly advance |i| until |e[i]==x|
      ++i;
    safeAdd(pols[i],klPol(x,ys),1);
    if (ERRNO) {
      Error(ERRNO,this,x,y);
      ERRNO = ERROR_WARNING;
      return;
    }
  }

  return;
} // |KLContext::KLHelper::add_second_terms|


/*
  This function subtracts the term P_{x,z}mu(z,ys) from the appropriate
  entry in pol. The idea is to extract [e,z] and run through it, so that
  we avoid calls to inOrder entirely.
*/
void KLContext::KLHelper::mu_correct_row
  (containers::vector<KLPol>& pols, coxtypes::CoxNbr y)
{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  coxtypes::Generator s = last(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);

  const MuRow& row = mu_row(ys);

  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu == 0)
      continue;
    coxtypes::CoxNbr z = row[j].x;
    coxtypes::Length h = row[j].height;
    klsupport::KLCoeff mu_value = row[j].mu;
    if (p.shift(z,s) > z)
      continue;

    auto b = p.closure(z);
    schubert::select_maxima_for(p,p.descent(y),b);

    Ulong i = 0;

    for (coxtypes::CoxNbr x : b)
    {
      while (e[i] < x)
	++i;
      safeSubtract(pols[i],klPol(x,z),mu_value,h+1);
      if (ERRNO) {
	Error(ERRNO,this,x,y);
	ERRNO = ERROR_WARNING;
	return;
      }
    }
  }

  return;
} // |KLContext::KLHelper::mu_correct_row|


/*
  Subtract the coatom correction from the list of polynomials.
  For each coatom |z| of |ys|, this is done for all the relevant x-es, as in
  |mu_correct_row|.
*/
void KLContext::KLHelper::coatom_correct_row
  (containers::vector<KLPol>& pol,coxtypes::CoxNbr y)
{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  coxtypes::Generator s = last(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);
  const schubert::CoxNbrList& c = p.hasse(ys);

  for (Ulong j = 0; j < c.size(); ++j)
  {
    coxtypes::CoxNbr z = c[j];
    if (p.shift(z,s) > z)
      continue;

    bitmap::BitMap b = p.closure(z);
    schubert::select_maxima_for(p,p.descent(y),b);

    Ulong i = 0;

    for (coxtypes::CoxNbr x : b)
    {
      while (e[i] < x)
	++i;
      safeSubtract(pol[i],klPol(x,z),1,1);
      if (ERRNO) {
	Error(ERRNO,this,x,y);
	ERRNO = ERROR_WARNING;
	return;
      }
    }
  }

  return;
} // |KLContext::KLHelper::coatom_correct_row|


/*
  Fill the |rows klList(y)| and |muList(y)|. It needs to be called when
  filling the table if |row_is_incomplete(y)| holds.

  This is one of the big functions in the program, of course. It is typically
  called when an element of the K-L basis of the Hecke algebra is required,
  or when we want to study the singularities of a Schubert variety. We have
  tried to optimize it for speed rather than memory efficiency.

  The outline of the function is as follows. Our descent strategy is to
  always chop off the last element in the normal form of the element under
  consideration, as provided by d_last. This is much better than taking,
  say, the first descent.

   - it is no longer assumed that the KL and mu rows for |y| (or for its inverse
     if that is numerically smaller) have been previously allocated: thus
     function starts by allocating the KL row it will compute if necessary, so
     everything will be allocated in due time during the recursion

   - then, we fill the corresponding rows, starting from the identity;
     the nice thing here is that it will always be guaranteed that
     P_{xs,ys}, P_{x,ys} and the mu(z,ys) are available; this is a recursive
     call to |compute_KL_row|, except that we can be sure that the rows have
     been allocated already.

   - then, we run through mu_row(ys), which has been filled, and determine
     the z for which mu(z,ys) != 0, and zs < z; these are the terms which
     cause corrections. For each such z, we call |compute_KL_row(z)|.

   - we do the same for the coatom list of y;

   - now, we are certain that all the terms we need for all the P_{x,y} in
     KL_row(y) are available. To avoid having to test the condition x <= z
     in the correction terms, we compute them all at once; again we run
     through the list of z's as above, and for all x <= z which is extremal
     w.r.t. the descent set of y, do the required subtraction. So we will
     work with a large vector of KLCoeff's, and need a table of offsets
     to the corresponding polynomials.

   One nice thing of this setup is that the work on the rows for y starts only
   when all the preparations are finished; it is guaranteed that there will not
   be any recursive calls in the duration of this filling. So, contrary to the
   compute_KL_pol function, workspace doesn't have to be managed as a stack.
   [That was only due to Fokko's insistance of using static variables for work
   space, and is no longer true. The extent to which this paragraph points to
   something that is still a serious advantage is not so clear. MvL]

  Returns the error ERROR_WARNING in case of failure, after printing an
  error message.
*/
void KLContext::KLHelper::compute_KL_row
  (coxtypes::CoxNbr y, bitmap::BitMap* done)
{
  if (y == 0)
    { if (done!=nullptr) done->insert(0); return; }

  if (inverse(y) < y) // then actually fill the row for |inverse(y)|
    y = inverse(y);

  if (row_needs_creation(y))
    create_KL_row(y);

  const schubert::SchubertContext& p = schubert();

  /* recursively fill in rows in the descent path */

  coxtypes::Generator s = last(y); // get last right descent for |y|
  coxtypes::CoxNbr ys = p.rshift(y,s);

  if (done==nullptr ? row_is_incomplete(ys)
      : done->is_member(ys) ? false
      : row_is_incomplete(ys) ? true
      : (done->insert(ys),false)
     )
  {
    compute_KL_row(ys,done); // recursive call
    if (ERRNO)
      goto abort;
  }

  // make sure the correcting terms are available
  ensure_correction_terms(y,s,done); // more recursive calls are hidden here
  if (ERRNO)
    goto abort;

  {
    auto pols = initial_polys(y,s);

    // add $q.P_{x,ys}$ when appropriate
    add_second_terms(pols,y);
    if (ERRNO)
      goto abort;

    // subtract correcting terms
    mu_correct_row(pols,y);
    if (ERRNO)
      goto abort;
    coatom_correct_row(pols,y);
    if (ERRNO)
      goto abort;

    // copy results to row in KL table
    writeKLRow(y,pols);
    if (ERRNO)
      goto abort;

    if (done!=nullptr)
      done->insert(y); // finally mark this row as fully computed
    return;
  }

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
} // |KLContext::KLHelper::compute_KL_row|

// Finally the driver function for filling the whole table row by row
void KLContext::KLHelper::fill_KL_table ()
{
  if (KL_is_filled())
    return;

  bitmap::BitMap done(size());
  for (coxtypes::CoxNbr y = 0; y < size(); ++y)
    if (y<=inverse(y))
    {
      compute_KL_row(y,&done);
      copy_mu_row_from_KL(y);
    }

  set_KL_filled();
} // |KLContext::KLHelper::fill_KL_table|

// ---------------------------------------------------------------------------

// The following function, when called from above "row" code, returns rapidly

const KLPol& KLContext::KLHelper::klPol (coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{ // we modify |x| and |y| as local variables, but $P(x,y)$ remains the same

  // Use the fact that $P(x,y)=P(x^{-1},y^{-1})$ to reduce necessary storage
  if (inverse(y) < y) {
    y = inverse(y);
    x = inverse(x);
  }

  const schubert::SchubertContext& p = schubert();

  // make |x| extramal for |y|
  x = p.maximize(x,p.descent(y));

  // check for (now) trivially small cases
  if (p.length(y) - p.length(x) < 3)
    return one();


  // make sure |KL_table[y]| is allocated (and hence |extr_list(y)| as well)
  if (row_needs_creation(y)) {
    create_KL_row(y);
    if (ERRNO)
      return zeroPol();
  }

  // find |x| in |extrList[y]|
  const auto& eL = extrList(y);
  Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();
  const KLPol*& pol = KL_row(y)[m];

  if (pol == nullptr) // should not happen after |compute_KL_row(y)| was done
  { // now  we have to compute the polynomial
    pol = compute_KL_pol(x,y,last(y));
    if (ERRNO)
      return zeroPol();
  }

  return *pol;
} // |KLHelper::klPol(x,y)|

/*
  This function performs one of the main tasks in the program, namely
  filling in a Kazhdan-Lusztig polynomial. In this function we are
  not concerned with maximal efficiency, but mostly with computing
  as few polynomials as possible in the recursion formula.

  It is assumed that x is extremal for y, and the $x\leq y$ in the Bruhat order

  Sets the error KL_FAIL, and returns the null pointer, in case of failure.

  NOTE : since this function potentially triggers the allocation of many
  rows in the klList and in the muList (not to mention the storage of the
  computed polynomials), it should labor under the constraint of
  CATCH_MEMORY_OVERFLOW.
*/
const KLPol* KLContext::KLHelper::compute_KL_pol
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s)
{
  const schubert::SchubertContext& p = schubert();

  assert(s != coxtypes::undef_generator); // choice of |s| is imposed upon us
  assert(p.length(y) - p.length(x) >= 3); // short cases were already taken out

  coxtypes::CoxNbr ys = p.shift(y,s);
  coxtypes::CoxNbr xs = p.shift(x,s);

  // unless $x<ys$, we simply have $P(x,y)=P(xs,ys)$, a single recursive call
  if (not p.Bruhat_leq(x,ys))
  { d_stats.klcomputed++; // count that call as one extra computation
    return &klPol(xs,ys); // recursion is (or may be) here, one call down
  }


  // "temporarily" be forgiving about running out of memory
  CATCH_MEMORY_OVERFLOW = true; // but this is not properly managed as stack

  KLPol pol = klPol(xs,ys); // potential recursion on call down, as above
  // N.B. quite likely CATCH_MEMORY_OVERFLOW has been switched off again here

  if (ERRNO) // maybe recursion ran out of memory
    goto abort;

  { // add $q.P_{x,ys}$
    const KLPol& p_xys = klPol(x,ys); // another recursive call
    if (ERRNO)
      goto abort;
    safeAdd(pol,p_xys,1);
    if (ERRNO)
      goto abort;
  }

  // subtract correction terms
  coatom_correct_KL_pol(x,y,s,pol);
  if (ERRNO)
    goto abort;

  mu_correct_KL_pol(x,y,s,pol);
  if (ERRNO)
    goto abort;

  { // look up and return pointer tot |pol|
    const KLPol* result = KL_pool.find(pol);
    if (ERRNO)
      goto abort;

    CATCH_MEMORY_OVERFLOW = false; // regardless of its value upon entry

    d_stats.klcomputed++;
    return result;
  }

 abort: /* an error occurred */

  CATCH_MEMORY_OVERFLOW = false;
  if (ERRNO != MEMORY_WARNING)
    ERRNO = KL_FAIL;
  return 0;

} // |KLContext::KLHelper::compute_KL_pol|


/*
  Fills row with the values for the mu(x,y). It is assumed that the row
  has been correctly allocated (i.e., that it holds only extremal entries
  with odd length difference > 1 w.r.t. y, and at least one entry for each
  non-zero mu.)
*/
void KLContext::KLHelper::fillMuRow(MuRow& row, const coxtypes::CoxNbr& y)
{
  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu == klsupport::undef_klcoeff) {
      coxtypes::CoxNbr x = row[j].x;
      row[j].mu = computeMu(x,y);
      if (ERRNO)
	return;
    }
  }

  return;
}


void KLContext::KLHelper::makeMuRow(const coxtypes::CoxNbr& y)

/*
  This function makes a row in the mu-table from scratch (i.e., it is assumed
  that isMuAllocated(y) = false). The idea is to do the minimal allocation
  (i.e., write down non-zero mu's only), and to compute as much as possible
  within the mu-table itself.

  This function is crucial to the optimal determination of the W-graph of
  the group, so a lot of care has been expended onto it. One crucial property
  of the mu-function is its stability under *-operations. Otherwise, one
  should remark that in the recursive formula for mu(x,y), there is just one
  term (the one coming from q.P_{x,ys}) that is not already a mu-value;
  and it will be a mu-value if x moves up under the descent set of ys.

  So the only entries that have to be computed are the ones that are extremal
  w.r.t. *-operations, and we will only need the highest-degree term of some
  kl-pols with even length-difference.

  NOTE : to be implemented!
*/

{}


/*
  Subtract from pol the correction terms

               P_{x,z}.mu(z,ys)q^{(l(y)-l(z))/2}

  corresponding to the z in [x,ys] s.t. mu(z,ys) != 0, zs < z, and l(ys)-l(z)
  >= 3; the correction for the coatoms in [x,ys] is dealt with in the
  function coatom_correct_KL_pol.

  NOTE : as usual, we watch out for memory overflow --- it is assumed that
  CATCH_MEMORY_OVERFLOW is turned on. Overflow in the coefficients cannot
  occur during subtraction; however, we set an error if a negative coefficient
  is found, as this would be a major discovery!

*/
void KLContext::KLHelper::mu_correct_KL_pol
  (const coxtypes::CoxNbr& d_x, const coxtypes::CoxNbr& y,
   const coxtypes::Generator& d_s, KLPol& pol)
{
  const schubert::SchubertContext& p = schubert();

  coxtypes::CoxNbr x = d_x;
  coxtypes::Generator s = d_s;
  coxtypes::CoxNbr ys = p.shift(y,s);

  if (!isMuAllocated(ys)) { /* allocate row */
    allocMuRow(ys);
    if (ERRNO)
      goto abort;
  }

  {
    MuRow& m = mu_row(ys);

    {
      coxtypes::Length ly = p.length(y);

      for (Ulong j = 0; j < m.size(); ++j) {

	coxtypes::CoxNbr z = m[j].x;

	if (p.shift(z,s) > z)
	  continue;
	if (not p.Bruhat_leq(x,z))
	  continue;

	/* compute the mu-coefficient if it was not already computed */

	if (m[j].mu == klsupport::undef_klcoeff) {
	  m[j].mu = computeMu(z,ys);
	  if (ERRNO)
	    goto abort;
	}

	/* subtract the correction if mu(z,ys) != 0 */

	if (m[j].mu) {
	  coxtypes::Length h = (ly - p.length(m[j].x))/2;

	  const KLPol& p_xz = klPol(x,z);
	  if (ERRNO)
	    goto abort;

	  safeSubtract(pol,p_xz,m[j].mu,h);
	  if (ERRNO)
	    goto abort;
	}
      }
    }
  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}



/*
  This function computes mu(x,y) using the general recursive formula for
  the descent s. In practice, this will be used only if the descent set
  of x contains the _second_ descent set of y. It is assumed that x < y
  x extremal w.r.t. y, and l(y)-l(x) odd > 1.

  Sets the error MU_FAIL, and returns the value undef_klcoeff, in case of
  failure (this can be due to memory overflow, or to coefficient over- or
  underflow.)
*/
klsupport::KLCoeff KLContext::KLHelper::recursiveMu
  (const coxtypes::CoxNbr& d_x, const coxtypes::CoxNbr& y,
			       const coxtypes::Generator& d_s)
{
  const schubert::SchubertContext& p = schubert();

  coxtypes::CoxNbr x = d_x;
  coxtypes::Generator s = d_s;

  coxtypes::Length l = p.length(y) - p.length(x); /* l is odd > 1 */

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);

  klsupport::KLCoeff r = mu(xs,ys);
  if (ERRNO)
    goto abort;

  if (not p.Bruhat_leq(x,ys)) { /* value is found */
    d_stats.mucomputed++;
    if (r == 0)
      d_stats.muzero++;
    return r;
  }

  /* special case when the length difference is three */

  if (l == 3) { /* P_{x,ys} = 1 */
    klsupport::safeAdd(r,1);
    if (ERRNO) { /* overflow; highly unlikely! */
      Error(MU_OVERFLOW,this,x,y);
      goto abort;
    }
    goto coatom_correction;
  }

  /* get the term from q.P_{x,ys} */

  {
    coxtypes::CoxNbr x1 = p.maximize(x,p.descent(ys));

    if (x != x1) { // no polynomial is needed; does not happen in double
                   // extremal case
      if (p.length(x1) == p.length(x) + 1) // no correction otherwise
	r += mu(x1,ys);
    }
    else { /* we need a polynomial */
      const KLPol& pol = klPol(x,ys);
      if (ERRNO)
	goto abort;
      polynomials::Degree d = (l-1)/2 - 1;
      if (pol.deg() == d) {
	klsupport::safeAdd(r,pol[d]);
	if (ERRNO) { /* overflow; highly unlikely! */
	  Error(MU_OVERFLOW,this,x,y);
	  goto abort;
	}
      }
    }
  }

  /* subtract correction terms where l(ys) - l(z) > 1 */

  {
    if (!isMuAllocated(ys)) { /* allocate row */
      allocMuRow(ys);
      if (ERRNO)
	goto abort;
    }

    MuRow& m = mu_row(ys);

    for (Ulong j = 0; j < m.size(); ++j) {

      coxtypes::CoxNbr z = m[j].x;

      if (z == x)
	continue;
      if (p.shift(z,s) > z)
	continue;
      if (not p.Bruhat_leq(x,z))
	continue;

      /* fill in the mu-coefficient if it was not already computed */

      if (m[j].mu == klsupport::undef_klcoeff) {
	m[j].mu = computeMu(z,ys);
	if (ERRNO)
	  goto abort;
      }

      /* subtract the correction if mu(z,ys) != 0 */

      if (m[j].mu) {
	klsupport::KLCoeff r1 = mu(x,z);
	if (ERRNO)
	  goto abort;
	klsupport::safeMultiply(r1,m[j].mu);
	if (ERRNO) { /* overflow; highly unlikely! */
	  Error(MU_OVERFLOW,this,x,y);
	  goto abort;
	}
	klsupport::safeSubtract(r,r1);
	if (ERRNO) { /* negative coefficient */
	  Error(MU_NEGATIVE,this,x,y);
	  goto abort;
	}
      }
    }
  }

 coatom_correction:

  /* subtract the coatom correction */

  {
    const schubert::CoxNbrList& c = p.hasse(ys);

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      coxtypes::CoxNbr zs = p.shift(z,s);
      if (zs > z)
	continue;
      if (not p.Bruhat_leq(x,z))
	continue;
      klsupport::KLCoeff r1 = mu(x,z);
      if (ERRNO)
	goto abort;
      klsupport::safeSubtract(r,r1);
      if (ERRNO) { /* negative coefficient */
	Error(MU_NEGATIVE,this,x,y);
	goto abort;
      }
    }
  }

  d_stats.mucomputed++;
  if (r == 0)
    d_stats.muzero++;

  return r;

 abort:
  if (ERRNO != MEMORY_WARNING)
    ERRNO = MU_FAIL;
  return klsupport::undef_klcoeff;
} // |KLContext::KLHelper::recursiveMu|


/*
  Write the polynomials from the list pol to KL_row(y); more precisely, it finds
  their adresses in klTree(), and writes those to KL_row(y). First it has to put
  the true degrees in the pol[j].

  It is assumed that y <= inverse(y).

  The only error that can occur here is memory overflow because of the
  allocation for new polynomials in KL_pool. In that case, the error is
  reported, and ERROR_WARNING is set.
*/
void KLContext::KLHelper::writeKLRow
  (const coxtypes::CoxNbr& y, containers::vector<KLPol>& pol)
{
  KLRow& kl_row = KL_row(y);

  for (Ulong j = 0; j < kl_row.size(); ++j)
    if (kl_row[j]==nullptr) // dont't overwrite already stored K-L polynomials
    {
      pol[j].snap_degree();
      kl_row[j] = KL_pool.find(pol[j]);
      if (kl_row[j] == nullptr)
      { /* an error occurred */
	Error(ERRNO);
	ERRNO = ERROR_WARNING;
	return; // give up storing remaining  polynomials
      }
      d_stats.klcomputed++;
    } // |for(j)| and |if(..)|
}


/*
  Copy nonzero entries of |row| to the corresponding row in the mu-list.
*/
void KLContext::KLHelper::writeMuRow
  (const MuRow& row, const coxtypes::CoxNbr& y)
{
  /* count non-zero entries */

  Ulong count = 0;

  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu != 0)
      count++;
  }

  MuRow& y_row = mu_row(y); // existing but cleared |MuRow| (?!)
  y_row.clear();
  y_row.reserve(count);
  if (ERRNO) {
    Error(ERRNO);
    ERRNO = ERROR_WARNING;
    return;
  }

  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu != 0)
      y_row.emplace_back(row[j].x,row[j].mu,row[j].height);
  }

  d_stats.munodes += count;
  d_stats.murows++;

} // |KLHelper::writeMuRow|


/****************************************************************************

        Chapter III -- The MuFilter class

  The MuFilter class is a small functor, useful for adapting iterators
  when constructing mu-lists. It filters away elements according with even
  length difference (with respect to supplied length), or length difference 1.
  The intention is to use it with the FilteredIterator adaptor.

 ****************************************************************************/


MuFilter::MuFilter(const schubert::SchubertContext& p, const coxtypes::Length& l)
  :d_p(p), d_l(l)
{}

MuFilter::MuFilter(const schubert::SchubertContext& p, const coxtypes::CoxNbr& y)
  :d_p(p), d_l(p.length(y))
{}


/****************************************************************************

        Chapter IV -- KLStats

  This section defines the functions declared for the KLStats structure :

   - KLStats() : constructor;
   - ~KLStats() : destructor;

 ****************************************************************************/


KLStats::KLStats()
  :klrows(0), klnodes(0), klcomputed(0), murows(0), munodes(0), mucomputed(0),
   muzero(0)
{}



/****************************************************************************

        Chapter V -- Input/Output

  This section defines the input/output functions declared in kl.h :

   - print(file,h) : prints a homology vector;
   - printMuTable(file,kl) : prints the full mu-table;
   - showKLPol(file,kl,x,y) : maps out the computation of P_{x,y};
   - showMu(file,kl,x,y) : maps out the computation of mu(x,y);
   - showSimpleMu(file,kl,x,y) : auxiliary to ShowMu;
   - showRecursiveMu(file,kl,x,y) : auxiliary to ShowMu;

 ****************************************************************************/

/*
  Print the stats of |kl|. This is a data structure that monitors precisely the
  computations that have been effected, and the memory that has been allocated.
*/
void print_stats(const KLContext& kl, FILE* file)
{
  const auto stat = kl.stats();
  fprintf(file,"klrows = %lu\n",stat.klrows);
  fprintf(file,"klnodes = %lu\n",stat.klnodes);
  fprintf(file,"klcomputed = %lu\n",stat.klcomputed);
  fprintf(file,"murows = %lu\n",stat.murows);
  fprintf(file,"munodes = %lu\n",stat.munodes);
  fprintf(file,"mucomputed = %lu\n",stat.mucomputed);
  fprintf(file,"muzero = %lu\n",stat.muzero);
}

void print(FILE* file, const schubert::Homology& h)

/*
  Prints the homology vector as a single long line.
*/

{
  if (h.size()) /* print first term */
    fprintf(file," h[0] = %lu",h[0]);

  for (Ulong j = 1; j < h.size(); ++j) {
    fprintf(file," h[%lu] = %lu",j,h[j]);
  }

  return;
}

void printMuTable(FILE* file, const KLContext& kl, const interface::Interface& I)

/*
  This function prints ou the contents of the mu-table. It prints only the
  entries for which mu(x,y) != 0 (and of course, for which l(y)-l(x)>1; the
  others are gotten from the coatom table.)
*/

{
  const schubert::SchubertContext& p = kl.schubert();

  for (coxtypes::CoxNbr y = 0; y < p.size(); ++y) {
    kl.print(file,y,I);
    fprintf(file," : ");
    const MuRow& row = kl.muList(y);
    Ulong count = 0;
    for (Ulong j = 0; j < row.size(); ++j) {
      const MuData& mu = row[j];
      if (mu.mu == 0)
	continue;
      if (count)
	fprintf(file,",");
      count++;
      fprintf(file,"{");
      fprintf(file,"x = ");
      kl.print(file,mu.x,I);
      fprintf(file,", mu = %lu, height = %lu",static_cast<Ulong>(mu.mu),
	      static_cast<Ulong>(mu.height));
      fprintf(file,"}");
    }
    fprintf(file,"\n");
  }

  return;
}


/*
  Print out the various terms appearing in the computation of the
  Kazhdan-Lusztig polynomial P_{x,y} through the standard recursion formula,
  using the generator |s| as descent generator.

  It is assumed that $x \eq y$ in the Bruhat order, and that $s$ is indeed a
  descent generator for $y$.
*/
void showKLPol
  (FILE* file, KLContext& kl,
   coxtypes::CoxNbr x, coxtypes::CoxNbr y,
   const interface::Interface& I, coxtypes::Generator s)
{
  const schubert::SchubertContext& p = kl.schubert();

  const KLPol& pol = kl.klPol(x,y,s);
  if (ERRNO) {
    Error (ERRNO);
    return;
  }

  const coxtypes::CoxNbr x_orig = x;

  unsigned long ls = io::LINESIZE;

  std::string buf;

  buf.append("x = ");
  p.append(buf,x,I);
  buf.append("; y = ");
  p.append(buf,y,I);
  buf.append(" L:");
  append(buf,p.ldescent(y),I);
  buf.append(" R:");
  append(buf,p.rdescent(y),I);
  io::foldLine(file,buf,ls,0,"yL");
  fprintf(file,"\n\n");

  if (kl.inverse(y) < y) { // go over to inverses
    x = kl.inverse(x);
    y = kl.inverse(y);
    fprintf(file,"inverse(y) < y\n");
    fprintf(file,"new x : ");
    p.print(file,x,I);
    fprintf(file,"\nnew y : ");
    p.print(file,y,I);
    fprintf(file,"\n\n");
  }

  Lflags f = p.descent(y);
  x = p.maximize(x,f);
  if (x > x_orig) {
      fprintf(file,"x is not extremal w.r.t. y\nnew x: ");
      p.print(file,x,I);
      fprintf(file,"\n\n");
  }

  coxtypes::Length d = p.length(y) - p.length(x);
  if (d < 3) { /* trivial case */
    fprintf(file,"l(y)-l(x) < 3\n\n");
    goto end;
  }

  {
    if (s == coxtypes::undef_generator)
      s = kl.last(y);
    coxtypes::CoxNbr xs = p.shift(x,s);
    coxtypes::CoxNbr ys = p.shift(y,s);

    if (not p.Bruhat_leq(x,ys)) { /* easy case */
      if (s < kl.rank())
	{ // action is on the right
	fprintf(file,"x not comparable to ys for s = %d\n",s+1);
	buf.clear();
	buf.append("xs = ");
	p.append(buf,xs,I);
	buf.append("; ys = ");
	p.append(buf,ys,I);
	io::foldLine(file,buf,ls,0,"y");
	fprintf(file,"\n\n");
	goto end;
      }
      else
      { // action is on the left
	fprintf(file,"x not comparable to sy for s = %d\n",s+1-kl.rank());
	buf.clear();
	buf.append("sx = ");
	p.append(buf,xs,I);
	buf.append("; sy = ");
	p.append(buf,ys,I);
	io::foldLine(file,buf,ls,0,"s");
	fprintf(file,"\n\n");
	goto end;
      }
    }

  /* apply recursion formula */

    if (s < kl.rank()) {
      fprintf(file,"applying recursion formula with s = %d on the right\n\n",
	      s+1);
      buf.clear();
      buf.append("xs = ");
      p.append(buf,xs,I);
      buf.append("; ys = ");
      p.append(buf,ys,I);
      io::foldLine(file,buf,ls,0,"y");
      fprintf(file,"\n\n");
    }
    else {
      fprintf(file,"applying recursion formula with s = %d on the left\n\n",
	      s+1-kl.rank());
      buf.clear();
      buf.append("sx = ");
      p.append(buf,xs,I);
      buf.append("; sy = ");
      p.append(buf,ys,I);
      io::foldLine(file,buf,ls,0,"s");
      fprintf(file,"\n\n");
    }

    /* first term */

    buf.clear();

    if (s < kl.rank()) {
      buf.append("P_{xs,ys} = ");
      append(buf,kl.klPol(xs,ys),"q");
    }
    else {
      buf.append("P_{sx,sy} = ");
      append(buf,kl.klPol(xs,ys),"q");
    }

    io::foldLine(file,buf,ls,4,"+");
    fprintf(file,"\n");

    /* second term */

    buf.clear();

    if (s < kl.rank()) {
      buf.append("P_{x,ys}  = ");
      append(buf,kl.klPol(x,ys),"q");
    }
    else {
      buf.append("P_{x,sy}  = ");
      append(buf,kl.klPol(x,ys),"q");
    }

    io::foldLine(file,buf,ls,4,"+");
    fprintf(file,"\n\n");

    /* coatom correction */

    const schubert::CoxNbrList& c = p.hasse(ys);
    bool coatomcorrection = false;

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      if (p.shift(z,s) > z)
	continue;
      if (not p.Bruhat_leq(x,z))
	continue;
      coatomcorrection = true;
      buf.clear();
      buf.append("z = ");
      p.append(buf,z,I);
      buf.append(" P_{x,z} = ");
      polynomials::append(buf,kl.klPol(x,z),"q");
      io::foldLine(file,buf,ls,4,"P+");
      fprintf(file,"\n");
    }

    if (coatomcorrection)
      fprintf(file,"\n");

    /* mu correction */

    const MuRow& m = kl.muList(ys);
    coxtypes::Length l_ys = p.length(ys);
    bool mucorrection = false;

    for (Ulong j = 0; j < m.size(); ++j) {
      coxtypes::CoxNbr z = m[j].x;
      if (p.shift(z,s) > z)
	continue;
      if (not p.Bruhat_leq(x,z))
	continue;
      if (m[j].mu) {
	mucorrection = true;
	buf.clear();
	buf.append("z = ");
	p.append(buf,z,I);
	io::pad(buf,l_ys+1); /* remember the four characters "z = " */
	buf.append(" mu = ");
	io::append(buf,m[j].mu);
	buf.append(" height = ");
	io::append(buf,m[j].height);
	buf.append(" P_{x,z} = ");
	append(buf,kl.klPol(x,z),"q");
	io::foldLine(file,buf,ls,4,"Pmh+");
	fprintf(file,"\n");
      }
    }

    if (mucorrection)
      fprintf(file,"\n");
  }

 end:

  buf.clear();
  buf.append("result : ");
  append(buf,pol,"q");
  if (2*pol.deg()+1 == d)
    buf.append(" *");
  io::foldLine(file,buf,ls,4,"+");
  fprintf(file,"\n\n");

  return;
}


/*
  Map out the computation of a mu-coefficient. See
  |KLContext::KLHelper::computeMu| for the algorithm.
*/
void showMu(FILE* file, KLContext& kl,
	    coxtypes::CoxNbr x, coxtypes::CoxNbr y,
	    const interface::Interface& I)
{
  std::string buf;

  const schubert::SchubertContext& p = kl.schubert();

  const coxtypes::CoxNbr x_orig = x;

  klsupport::KLCoeff r = kl.mu(x,y);
  if (ERRNO)
    goto abort;

  {
    unsigned long ls = io::LINESIZE;
    buf.append("x = ");
    p.append(buf,x,I);
    buf.append("  y = ");
    p.append(buf,y,I);
    buf.append(" L:");
    append(buf,p.ldescent(y),I);
    buf.append(" R:");
    append(buf,p.rdescent(y),I);
    io::foldLine(file,buf,ls,0,"yL");
    fprintf(file,"\n\n");

    Lflags fy = p.descent(y);
    x = p.maximize(x,fy);
    if (x > x_orig) {
      fprintf(file,"x is not extremal w.r.t. y\n\nresult: 0\n\n");
      return;
    }

    coxtypes::Length d = p.length(y) - p.length(x);
    if ((d%2) == 0) {
      fprintf(file,"even length difference\n\nresult: 0\n\n");
      return;
    }
    if (d == 1) { /* trivial case */
      fprintf(file,"x is coatom of y\n\nresult: 1\n\n");
      return;
    }

    Lflags f2 = p.twoDescent(y);
    Lflags fx = p.descent(x);

    if ((fx&f2) != f2) { /* recursion case */

      fprintf(file,"x is not doubly extremal w.r.t. y\n\n");
      showSimpleMu(file,kl,x,y,r,I);

      return;
    }

    fprintf(file,"x is doubly extremal w.r.t. y\n\n");
    showRecursiveMu(file,kl,x,y,r,I);

    return;

  }

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;

}

namespace {


/*
  Maps out the computation in the case where the recursion formula has to
  be used.
*/
void showRecursiveMu(FILE* file, KLContext& kl,
		     coxtypes::CoxNbr x, coxtypes::CoxNbr y,
		     klsupport::KLCoeff r, const interface::Interface& I)
{
  std::string buf;

  const schubert::SchubertContext& p = kl.schubert();
  unsigned long ls = io::LINESIZE;

  coxtypes::Generator s = kl.last(y);
  coxtypes::Length l = p.length(y) - p.length(x); /* l is odd > 1 */

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);

  if (not p.Bruhat_leq(x,ys)) { // mu(x,y) = mu(xs,ys)
    if (s < kl.rank()) { // action is on the right
      fprintf(file,"x not comparable to ys for s = %d\n",s+1);
      buf.clear();
      buf.append("xs = ");
      p.append(buf,xs,I);
      buf.append("; ys = ");
      p.append(buf,ys,I);
      io::foldLine(file,buf,ls,0,"y");
      fprintf(file,"\n\nresult : %lu\n\n",static_cast<Ulong>(r));
    }
    else { // action is on the left
      fprintf(file,"x not comparable to sy for s = %d\n",s-kl.rank()+1);
      buf.clear();
      buf.append("sx = ");
      p.append(buf,xs,I);
      buf.append("; sy = ");
      p.append(buf,ys,I);
      io::foldLine(file,buf,ls,0,"s");
      fprintf(file,"\n\n");
      fprintf(file,"\n\nresult : %lu\n\n",static_cast<Ulong>(r));
    }
    return;
  }

  // if we get to this point, w need to apply the full recursion formula

    if (s < kl.rank()) {
      fprintf(file,"applying recursion formula with s = %d on the right\n\n",
	      s+1);
      buf.clear();
      buf.append("xs = ");
      p.append(buf,xs,I);
      buf.append("; ys = ");
      p.append(buf,ys,I);
      io::foldLine(file,buf,ls,0,"y");
      fprintf(file,"\n\n");
    }
    else {
      fprintf(file,"applying recursion formula with s = %d on the left\n\n",
	      s+1-kl.rank());
      buf.clear();
      buf.append("sx = ");
      p.append(buf,xs,I);
      buf.append("; sy = ");
      p.append(buf,ys,I);
      io::foldLine(file,buf,ls,0,"s");
      fprintf(file,"\n\n");
    }

  // first term

    buf.clear();

    if (s < kl.rank()) {
      fprintf(file,"mu(xs,ys) = %lu\n",static_cast<Ulong>(kl.mu(xs,ys)));
    }
    else {
      fprintf(file,"mu(sx,sy) = %lu\n",static_cast<Ulong>(kl.mu(xs,ys)));
    }

  // second term

    buf.clear();

    const KLPol& pol = kl.klPol(x,ys);
    klsupport::KLCoeff r1 = 0;
    polynomials::Degree d = (l-1)/2 - 1;
    if (pol.deg() == d) {
      r1 = pol[d];
    }

    if (s < kl.rank()) {
      fprintf(file,"second term is %lu\n",static_cast<Ulong>(r1));
    }
    else {
      fprintf(file,"second term is %lu\n",static_cast<Ulong>(r1));
    }

    fprintf(file,"\n");

  // coatom correction

    const schubert::CoxNbrList& c = p.hasse(ys);
    bool coatomcorrection = false;

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      if (p.shift(z,s) > z)
	continue;
      if (not p.Bruhat_leq(x,z))
	continue;
      coatomcorrection = true;
      buf.clear();
      buf.append("z = ");
      p.append(buf,z,I);
      buf.append(" mu(x,z) = ");
      io::append(buf,kl.mu(x,z));
      io::foldLine(file,buf,ls,4," ");
      fprintf(file,"\n");
    }

    if (coatomcorrection)
      fprintf(file,"\n");

  // mu-correction

   const MuRow& m = kl.muList(ys);
   coxtypes::Length l_ys = p.length(ys);
   bool mucorrection = false;

   for (Ulong j = 0; j < m.size(); ++j) {

     coxtypes::CoxNbr z = m[j].x;
     if (p.shift(z,s) > z)
       continue;
     if (not p.Bruhat_leq(x,z))
       continue;

     // fill in the mu-coefficient if it was not already computed

     if (m[j].mu == klsupport::undef_klcoeff) {
       kl.mu(z,ys); // this will fill m[j].mu
       if (ERRNO) {
	 Error(ERRNO);
	 return;
       }
     }

     if (m[j].mu) {
       mucorrection = true;
       buf.clear();
       buf.append("z = ");
       p.append(buf,z,I);
       io::pad(buf,l_ys+1); // remember the four characters "z = "
       buf.append(" mu = ");
       io::append(buf,m[j].mu);
       buf.append(" height = ");
       io::append(buf,m[j].height);
       buf.append(" mu(x,z) = ");
       io::append(buf,kl.mu(x,z));
       io::foldLine(file,buf,ls,4," ");
       fprintf(file,"\n");
     }

   }

   if (mucorrection)
     fprintf(file,"\n");

   // print result :

   fprintf(file,"result : %lu\n\n",static_cast<Ulong>(r));

   return;
}


/*
  Auxiliary to ShowMu. Maps out the computation of a mu-coefficient in the
  case where there is a direct recursion. It is assumed that the computation
  proper has already been tried with success.
*/
void showSimpleMu(FILE* file, KLContext& kl, coxtypes::CoxNbr x,
		  coxtypes::CoxNbr y, klsupport::KLCoeff r,
		  const interface::Interface& I)
{
  std::string buf;

  const schubert::SchubertContext& p = kl.schubert();
  unsigned long ls = io::LINESIZE;

  coxtypes::Generator s=-1,t=-1; // init silences compiler warnings
  Lflags desc;
  for (desc = p.descent(y); desc; desc &= desc-1) {
    coxtypes::Generator u = constants::firstBit(desc);
    coxtypes::CoxNbr yu = p.shift(y,u);
    Lflags fu = p.descent(yu);
    if ((p.descent(x)&fu) != fu) {
      s = u;
      t = constants::firstBit(fu & ~p.descent(x));
      break;
    }
  }
  assert(desc!=0);

  fprintf(file,"using descent s = %lu and ascent t = %lu\n\n",
	  static_cast<Ulong>(s+1),static_cast<Ulong>(t+1));

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);
  coxtypes::CoxNbr xt = p.shift(x,t);
  coxtypes::CoxNbr yst = p.shift(ys,t);

  /* consider four cases */

  buf.clear();

  if (p.descent(xt) & constants::eq_mask[s]) { /* xts < xt */

    buf.append("xs = ");
    p.append(buf,xs,I);
    buf.append("  ys = ");
    p.append(buf,ys,I);

    if (p.descent(yst) & constants::eq_mask[s]) { /* ysts < yst */
      buf.append("  yst = ");
      p.append(buf,yst,I);
      io::foldLine(file,buf,ls,0,"xy");
      fprintf(file,"\n\n");
      fprintf(file,
	      "result is mu(xs,ys)-mu(x,yst) = %lu - %lu = %lu\n\n",
	      static_cast<Ulong>(kl.mu(xs,ys)),
	      static_cast<Ulong>(kl.mu(x,yst)),
	      static_cast<Ulong>(r));
      return;
    }
    else { /* ysts > yst */
      io::foldLine(file,buf,ls,0,"xy");
      fprintf(file,"\n\n");
      fprintf(file,"result is mu(xs,ys) = %lu\n\n",static_cast<Ulong>(r));
      return;
    }
  }
  else { /* xts > xt */
    if (p.descent(yst) & constants::eq_mask[s]) { /* ysts < yst */
      buf.append("xs = ");
      p.append(buf,xs,I);
      buf.append("  xt = ");
      p.append(buf,xt,I);
      buf.append("  ys = ");
      p.append(buf,ys,I);
      buf.append("  yst = ");
      p.append(buf,yst,I);
      io::foldLine(file,buf,ls,0,"xy");
      fprintf(file,"\n\n");
      fprintf(file,
	"result is mu(xs,ys)+mu(xt,ys)-mu(x,yst) = %lu + %lu - %lu = %lu\n\n",
	      static_cast<Ulong>(kl.mu(xs,ys)),
	      static_cast<Ulong>(kl.mu(xt,ys)),
	      static_cast<Ulong>(kl.mu(x,yst)),
	      static_cast<Ulong>(r));
      return;
    }
    else { /* ysts > yst */
      if (p.descent(xs) & constants::eq_mask[t]) {
	buf.append("xs = ");
	p.append(buf,xs,I);
	buf.append("  xt = ");
	p.append(buf,xt,I);
	buf.append("  ys = ");
	p.append(buf,ys,I);
	io::foldLine(file,buf,ls,0,"xy");
	fprintf(file,"\n\n");
	fprintf(file,"result is mu(xs,ys)+mu(xt,ys) = %lu + %lu = %lu\n\n",
		static_cast<Ulong>(kl.mu(xs,ys)),
		static_cast<Ulong>(kl.mu(xt,ys)),
		static_cast<Ulong>(r));
	return;
      }
      else { /* mu(xs,ys) = 0 */
	buf.append("xt = ");
	p.append(buf,xt,I);
	buf.append("  ys = ");
	p.append(buf,ys,I);
	io::foldLine(file,buf,ls,0,"xy");
	fprintf(file,"\n\n");
	fprintf(file,"result is mu(xt,ys) = %lu\n\n", static_cast<Ulong>(r));
	return;
      }
    }
  }

  return;
} // |showSimpleMu|

}; // |namespace|

/****************************************************************************

        Chapter VI -- Singularities

  This section provides some functions for the study of the (rational)
  singular locus of Schubert varieties. These functions can be defined
  for arbitrary Coxeter groups, but of course the geometrical interpretation
  only makes sense when the Schubert variety is known to exist.

  The following functions are defined :

   - genericSingularities(h,y,kl) : returns the singular locus;
   - isSingular(row) : checks rational singularity;

   The function singularStratification has been moved (perhaps mistakenly)
   to hecke.cpp.

 *****************************************************************************/


/*
  Return in h the singular locus of cl(X_y). The point is to
  do this while computing as few K-L polynomials as possible.
*/
void genericSingularities(HeckeElt& h, coxtypes::CoxNbr y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  auto b = p.closure(y);
  schubert::select_maxima_for(p,p.descent(y),b);

  h.clear();

  Ulong x=p.size(); // prepare for decreasing traversal
  while (b.back_up(x))
  {
    const KLPol& pol = kl.klPol(x,y);
    if (ERRNO)
      return;
    if (pol.deg() > 0)
    { // remove the almost-interval [e,x[ from |b|
      h.emplace_back(x,&pol);
      b.andnot(p.closure(x)); // never mind it removes |x|; this bird has flown
    }
  }
  std::reverse(h.begin(),h.end());

  return;
}


/*
  Whether one of the polynomials in the row is distinct from unity. This is
  equivalent to the existence of a rational rational singularity of the Schubert
  variety, when such a geometric context is defined, and the row is the extremal
  row for an element |y|.

  NOTE : conjecturally, annulation of the term corresponding to the
  extremalization of the origin ensures annulation of all the others.
*/

bool isSingular(const HeckeElt& h)
{
  for (Ulong j = 0; j < h.size(); ++j) {
    const KLPol& pol = h[j].pol();
    if (pol.deg() != 0)
      return true;
  }

  return false;
}


/*
  Whether one of the polynomials in the row is distinct from unity. This is
  equivalent to the existence of a rational singularity of the Schubert variety,
  when such a geometric context is defined, and the row is the extremal row for
  an element |y|.
*/
bool isSingular(const KLRow& row)
{
  for (Ulong j = 0; j < row.size(); ++j) {
    const KLPol* pol = row[j];
    if (pol->deg() != 0)
      return true;
  }

  return false;
}


/****************************************************************************

        Chapter VII -- schubert::Homology vectors

  This section contains code for the computation of homology vectors. For
  now, we follow a simple-minded approach, making the vector element-wise.
  Something more sophisticated could be done, using the orbits in the
  Schubert closure under the action of the descent set of y; this would
  involve recognizing types of finite groups and such things.

  The following functions are defined :

   - h = ihBetti(kl,row) : puts the IH betti numbers of row in h;

 *****************************************************************************/

/*
  This function puts the IH betti numbers of the row in h, in a simple-minded
  approach. There is a serious danger of ineteger overflow here. For now, we set
  the value to undef_coxsize in case of overflow; this should be improved of
  course if it happens too often, using long integers for instance.
*/
schubert::Homology ihBetti(coxtypes::CoxNbr y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  auto cl = p.closure(y);

  schubert::Homology result(p.length(y)+1,0);

  for (auto it = cl.begin(); it(); ++it)
  {
    const KLPol& pol = kl.klPol(*it,y);
    coxtypes::Length d = p.length(*it);
    for (Ulong i = 0; i <= pol.deg(); ++i) {
      if (result[d+i] > coxtypes::COXSIZE_MAX - pol[i])
	result[d+i] = coxtypes::undef_coxnbr;
      else
	result[d+i] += pol[i];
    }
  }
  return result;
}


/****************************************************************************

        Chapter VIII -- Kazhdan-Lustig bases.

  This section defines functions returning Kazhdan-Lusztig bases. Note that
  the coefficients of these bases should actually be Laurent polynomials;
  however we (perhaps mistakenly) leave it for now to the output functions
  to do the shifting that is required; this saves us from introducing
  a new type at this point.

  The following functions are defined :

    - cBasis(h,y,kl) : also called sometimes the C'-basis; in our opinion,
      the right choice of basis;

 ****************************************************************************/



/*
  This is what in the original Kazhdan-Lusztig paper is called the C'-basis,
  but is now usually denoted c. The C-basis from the K-L paper doesn't seem
  worth implementing.
*/
HeckeElt cBasis(coxtypes::CoxNbr y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  bitmap::BitMap b = p.closure(y);

  HeckeElt result; result.reserve(b.size());

  for (coxtypes::CoxNbr x : b)
    result.emplace_back(x,&kl.klPol(x,y));
  return result;
}


/****************************************************************************

        Chapter IX -- Utility functions

  This section defines some utility functions used in this module :

   - appendStar(str,kl,x,pol,l) : appends a star if mu != 0;
   - find(row,x) : finds x in an extremal row (for KLRow or MuRow);
   - permuteValues(kl,q) : permutes the values of tables in kl;
   - permuteRanges(kl,q) : permutes the ranges of tables in kl;
   - safeAdd(p,q,n) : adds x^n.q to p, checking for overflow;
   - safeSubtract(p,q,mu,h) : subtracts x^h.mu.q from p, checking for
     underflow;
   - zeroPol() : returns the zero-polynomial;

 ****************************************************************************/

namespace {


/*
  Find |x| in |row| and return the address of the corresponding row.
  Returns null pointer if |x| is not found. Uses binary search.
*/
MuData* find(MuRow& row, const coxtypes::CoxNbr& x)
{
  Ulong j0 = (Ulong)(-1);

  for (Ulong j1 = row.size(); j1-j0 > 1;) {
    Ulong j = j0 + (j1-j0)/2;
    if (row[j].x == x) // m was found
      return &row[j];
    if (row[j].x < x)
      j0 = j;
    else
      j1 = j;
  }

  return nullptr;
}

}; // |namespace|


const KLPol& one()

{
  static KLPol p(1,KLPol::const_tag());
  return p;
}

namespace {


/*
  Increment |p| by $q*X^n$, checking for overflow.

  Forwards the error KLCOEFF_OVERFLOW in case of error.
*/
KLPol& safeAdd(KLPol& p, const KLPol& q, const polynomials::Degree& n)
{
  p.ensure_degree(q.deg() + n);

  for (polynomials::Degree j = 0; j <= q.deg(); ++j) {
    klsupport::safeAdd(p[j+n],q[j]);
    if (ERRNO)
      return p;
  }

  return p;
}


/*
  This function subtracts mu times q shifted by x^h from p, checking for
  underflow.
  Sets the error KLCOEFF_NEGATIVE in case of problem.
*/
KLPol& safeSubtract(KLPol& p, const KLPol& q, const klsupport::KLCoeff& mu,
		    const coxtypes::Length& h)
{
  for (polynomials::Degree j = 0; j <= q.deg(); ++j) {
    klsupport::KLCoeff a = mu;
    klsupport::safeMultiply(a,q[j]);
    if (ERRNO) { /* overflow; this will cause an underflow */
      ERRNO = KLCOEFF_NEGATIVE;
      return p;
    }
    klsupport::safeSubtract(p[j+h],a);
    if (ERRNO)
      return p;
  }

  p.snap_degree();

  return p;
}


/*
  Returns the zero polynomial (usually this indicates an error condition.)
*/
KLPol& zeroPol()
{
  static KLPol z(polynomials::undef_degree);
  return z;
}

} // |namespace|

} // |namespace kl|

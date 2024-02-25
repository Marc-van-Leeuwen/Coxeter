/*
  This is kl.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "kl.h"

#include "error.h"
#include "iterator.h"

namespace kl {
  using namespace error;
  using namespace iterator;


/****************************************************************************

  This file contains the code for the computation of the (ordinary) Kazhdan-
  Lusztig polynomials for W.

  The K-L context maintains two main lists : klList and muList. Each of these
  has one pointer (to a KLRow and a MuRow, respectively) for each element in
  the current schubert context; the enlargement and permutation functions
  for the schubert context make sure that the various active K-L contexts
  are enlarged/permuted accordingly.

  The row y in klList, when allocated, contains one entry for each extremal
  x <= y; the list of those x can be gotten from the klsupport structure,
  which maintains these extremal lists. The lists of extremal pairs are shared
  among all K-L contexts. Each entry in klList(y) is simply a pointer to
  a (constant) polynomial, the K-L polynomial for (x,y). The main complication
  is that I've (rightly or wrongly) decided to not allocate klList(y) if
  inverse(y) < y; in that case, the polynomial is read from the list for
  inverse(y). This saves space, but makes lookup more complicated. More
  importantly though, the fact of systematically going over to inverse(y) when
  it is smaller, seems to make the algorithm quite a bit faster and leaner
  (in the sense that fewer rows need to be computed); this, more than the
  memory saving, is my main reason for keeping this complication in.

  The situation for mu-rows is a little bit more delicate. First of all, we
  do not here go over to inverses; we fill all the rows that are needed
  (using the mu-value for the inverses if it is available.) Further, we only
  look at pairs where the length difference is at least three; it is known
  that we can then reduce to extremal pairs, the other mu-values being zero.
  But in fact, even among extremal pairs, most values are zero as well. So
  we do not necessarily allocate one entry for each pair; what is guaranteed
  is that all entries correspond to extremal pairs, and all non-zero mu-values
  have an entry. Ideally, we would wish to allocate only the non-zero mu's;
  however, this would require computing a full row as soon as one coefficient
  is needed, so we have refrained from that. Anyway, to find out the value
  for mu(x,y), we extremalize x, check that the length difference is odd >= 3,
  look x up in muList(y), and return the corresponding mu-value if found, 0
  otherwise.

  The idea is to compute everything upon request : we compute exactly what
  is needed in the recursive computation of what is required.

  The requests that we mainly focus on are the following :

    - compute a single K-L polynomial (or mu-coefficient);
    - compute all P_{x,y}'s (or mu(x,y)'s) for a given y;
    - compute the full file of K-L polynomials (mu coefficients);

  It turns out that the dominant factor in the computation is the Bruhat order
  comparison. That is why computing a whole row can be done much more
  efficently than computing each entry in the row separately : using our
  closure function from the schubert context, we can basically factor out
  most of the Bruhat things for a whole row. In any case, the row computation
  can be done without a single call to the slow inOrder function!

 ****************************************************************************/

struct KLContext::KLHelper
{
// data
  KLContext& d_kl;
  klsupport::KLSupport& d_klsupport; // unowned, |CoxGroup| owns it

  containers::bag<KLPol> d_klTree;
  KLStats d_stats;

// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
      {return memory::arena().free(ptr,sizeof(KLHelper));}
  KLHelper(klsupport::KLSupport& kls,KLContext* kl)
    : d_kl(*kl)
    , d_klsupport(kls)
    , d_klTree()
    , d_stats()
  {}

// methods
// relay methods
  KLStats& status() { return d_stats; }

    void allocExtrRow(const coxtypes::CoxNbr& y) {klsupport().allocExtrRow(y);}
    void allocKLRow(const coxtypes::CoxNbr& y);
    void allocMuRow(const coxtypes::CoxNbr& y);
    void allocMuRow(MuRow& row, const coxtypes::CoxNbr& y);
    void allocMuTable();
    void allocRowComputation(const coxtypes::CoxNbr& y);
    bool checkKLRow(const coxtypes::CoxNbr& y);
    bool checkMuRow(const coxtypes::CoxNbr& y);
    void coatomCorrection(const coxtypes::CoxNbr& y, list::List<KLPol>& pol);
    void coatomCorrection
      (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
       const coxtypes::Generator& s,
       list::List<KLPol>& pol, const Ulong& a);
    klsupport::KLCoeff computeMu
      (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y);
    const klsupport::ExtrRow& extrList(const coxtypes::CoxNbr& y)
      { return klsupport().extrList(y); }
    void fillMuRow(MuRow& row, const coxtypes::CoxNbr& y);
    const KLPol* fillKLPol
      (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
       const coxtypes::Generator& s = coxtypes::undef_generator);
    void fillKLRow(const coxtypes::CoxNbr& y);
    void initWorkspace(const coxtypes::CoxNbr& y, list::List<KLPol>& pol);
    coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& y) {return klsupport().inverse(y);}
    void inverseMuRow(const coxtypes::CoxNbr& y);
    bool isExtrAllocated(const coxtypes::CoxNbr& y)
      { return klsupport().isExtrAllocated(y); }
    bool isKLAllocated(const coxtypes::CoxNbr& y)
      { return d_kl.isKLAllocated(y); }
    bool isMuAllocated(const coxtypes::CoxNbr& y) {return d_kl.isMuAllocated(y);}
    KLRow& klList(const coxtypes::CoxNbr& y) {return *d_kl.d_klList[y];}
    klsupport::KLSupport& klsupport() {return d_klsupport;}
    containers::bag<KLPol>& klTree() {return d_klTree;}
    coxtypes::Generator last(const coxtypes::CoxNbr& x) {return klsupport().last(x);}
    void makeMuRow(const coxtypes::CoxNbr& y);
    void muCorrection(const coxtypes::CoxNbr& y, list::List<KLPol>& pol);
    void muCorrection(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s,
		      list::List<KLPol>& pol, const Ulong& a);
    MuRow& muList(const coxtypes::CoxNbr& y) {return *d_kl.d_muList[y];}
    void prepareRow(const coxtypes::CoxNbr& y, const coxtypes::Generator& s);
    coxtypes::Rank rank() {return d_kl.rank();}
    void readMuRow(const coxtypes::CoxNbr& y);
    klsupport::KLCoeff recursiveMu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s);
    const schubert::SchubertContext& schubert() {return klsupport().schubert();}
    void secondTerm(const coxtypes::CoxNbr& y, list::List<KLPol>& pol);
    Ulong size() {return d_kl.d_klList.size(); }
    void writeKLRow(const coxtypes::CoxNbr& y, list::List<KLPol>& pol);
    void writeMuRow(const MuRow& row, const coxtypes::CoxNbr& y);

  void grow(Ulong prev, Ulong n);
  void shrink(const Ulong& n);

  void fill_KL_table ();
  void fill_mu_table ();

  const KLPol& klPol(coxtypes::CoxNbr x, coxtypes::CoxNbr y);
  const KLPol& klPol(coxtypes::CoxNbr x, coxtypes::CoxNbr y,
		     coxtypes::Generator s);
  klsupport::KLCoeff mu(coxtypes::CoxNbr x, coxtypes::CoxNbr y);

  HeckeElt KL_row_as_HeckeElt(coxtypes::CoxNbr y);
  void move_KL_row_to_inverse(coxtypes::CoxNbr x);
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
  void showSimpleMu(FILE* file, KLContext& kl, const coxtypes::CoxNbr& x,
		    const coxtypes::CoxNbr& y, const klsupport::KLCoeff& r,
		    const interface::Interface& I);
  void showRecursiveMu(FILE* file, KLContext& kl, const coxtypes::CoxNbr& x,
		       const coxtypes::CoxNbr& y, const klsupport::KLCoeff& r, const interface::Interface& I);
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
  : d_klList(kls->size()), d_muList(kls->size())
{
  d_help = new KLHelper(*kls,this);

  d_klList.setSizeValue(kls->size());
  d_klList[0] = new KLRow(1);
  d_klList[0]->setSizeValue(1);
  (*d_klList[0])[0] = d_help->d_klTree.find(one());

  stats().klnodes++;
  stats().klrows++;
  stats().klcomputed++;

  d_muList[0].reset(new MuRow(0));
}


/*
  Recall that a KLContext is assumed to own its schubert::SchubertContext. The
  only destructions that are not done automatically are those of the Schubert
  context, the status pointer, and the |klList| and |muList| pointers.
*/
KLContext::~KLContext()
{
  for (Ulong j = 0; j < size(); ++j)
    delete d_klList[j];
}

/******** accessors **********************************************************/

const klsupport::KLSupport& KLContext::klsupport() const
  { return d_help->d_klsupport; }
KLStats& KLContext::stats() { return d_help->d_stats; }
const KLStats& KLContext::stats() const { return d_help->d_stats; }
Ulong KLContext::size() const { return d_help->size(); }
const KLRow& KLContext::klList(const coxtypes::CoxNbr& y) const
  { return *d_klList[y]; }
const MuRow& KLContext::muList(const coxtypes::CoxNbr& y) const
  { return *d_muList[y]; }

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

  d_kl.d_klList.setSize(n);
  if (ERRNO)
    goto revert;
  d_kl.d_muList.resize(n);
  if (ERRNO)
    goto revert;

  CATCH_MEMORY_OVERFLOW = false;

  d_kl.clearFullKL();
  d_kl.clearFullMu();

  return;

 revert:
  CATCH_MEMORY_OVERFLOW = false;
  d_kl.revertSize(prev_size);
} // |KLContext::KLHelper::grow|




// This function fills all the rows in klList, in a straightforward way.
void KLContext::fillKL() { d_help->fill_KL_table(); }

void KLContext::KLHelper::fill_KL_table ()
{
  if (d_kl.isFullKL())
    return;

  for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) < y)
      continue;
    if (not isKLAllocated(y))
      allocKLRow(y);
    fillKLRow(y);
    readMuRow(y);
  }

  d_kl.setFullKL();
} // |KLContext::KLHelper::fill_KL_table|


/*
  This function fills all the rows in the mu-list, in a straightforward way.
  Sets the error ERROR_WARNING in case of error.

  NOTE : error handling should be improved!
*/
void KLContext::fillMu() { d_help->fill_mu_table(); }

void KLContext::KLHelper::fill_mu_table ()
{
  if (d_kl.isFullMu())
    return;

  static MuRow mu_row(0);

  allocMuTable();

  if (ERRNO)
    goto abort;

  for (coxtypes::CoxNbr y = 0; y < size(); ++y) {
    if (inverse(y) < y)
      inverseMuRow(inverse(y));
    fillMuRow(*d_kl.d_muList[y],y);
    if (ERRNO)
      goto abort;
  }

  d_kl.setFullMu();
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


  // make sure |d_klList[y]| is allocated (and hence |extrList(y)| as well)
  if (not isKLAllocated(y)) {
    allocKLRow(y);
    if (ERRNO)
      return zeroPol();
  }

  // find |x| in |extrList[y]|
  const auto& eL = extrList(y);
  Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();
  const KLPol*& pol = klList(y)[m];

  if (pol == nullptr) { // now  we have to compute the polynomial
    pol = fillKLPol(x,y,last(y));
    if (ERRNO)
      return zeroPol();
  }

  return *pol;
} // |KLHelper::klPol(x,y)|

const KLPol& KLContext::klPol
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s)
{ return d_help->klPol(x,y,s); }

const KLPol& KLContext::KLHelper::klPol
  (coxtypes::CoxNbr x, coxtypes::CoxNbr y, coxtypes::Generator s)
{ // we modify |x| and |y| as local variables, but $P(x,y)$ remains the same

  // we are asked for an explicit right descent |s|, so we cannot use inverse

  const schubert::SchubertContext& p = schubert();

  // make |x| extramal for |y|
  x = p.maximize(x,p.descent(y));

  // check for (now) trivially small cases
  if (p.length(y) - p.length(x) < 3) { /* result is 1 */
    return one();
  }

  // make sure |d_klList[y]| is allocated (and hence |extrList(y)| as well)
  if (not isKLAllocated(y)) {
    allocKLRow(y);
    if (ERRNO)
      return zeroPol();
  }

  // find |x| in |extrList[y]|
  const auto& eL = extrList(y);
  Ulong m = std::lower_bound(eL.begin(),eL.end(),x)-eL.begin();
  const KLPol*& pol = klList(y)[m];

  if (pol == nullptr) { /* we have to compute the polynomial */
    pol = fillKLPol(x,y,s);
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

  // allocate |*d_muList[y]| if necessary
  if (not isMuAllocated(y)) {
    allocMuRow(y);
    if (ERRNO)
      return klsupport::undef_klcoeff;
  }

  // find x in |*d_muList[y]|
  MuRow& m = *d_kl.d_muList[y];
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
  if (not checkKLRow(y)) {
    allocRowComputation(y);
    fillKLRow(y);
  }
  if (ERRNO) {
    Error(ERRNO);
    ERRNO = ERROR_WARNING;
    return h;
  }

  if (y <= inverse(y)) {
    const klsupport::ExtrRow& e = extrList(y);
    h.reserve(e.size());
    const KLRow& klr = klList(y);
    for (Ulong j = 0; j < e.size(); ++j)
      h.emplace_back(e[j],klr[j]);
  }
  else { /* go over to inverses */
    coxtypes::CoxNbr yi = inverse(y);
    const klsupport::ExtrRow& e = extrList(yi);
    h.reserve(e.size());
    const KLRow& klr = klList(yi);
    for (Ulong j = 0; j < e.size(); ++j)
      h.emplace_back(inverse(e[j]),klr[j]);

    std::sort(h.begin(),h.end()); // make sure list is ordered
  }
  return h;
} // |::KLHelper::KL_row_as_HeckeElt|


/*
  Exchange rows for |x| and |inverse(x)| in |d_klList|. It is assumed that the
  boths rows are within the bounds of |d_klList|.
*/
void KLContext::applyInverse(const coxtypes::CoxNbr& x)
{ d_help->move_KL_row_to_inverse(x);
}

void KLContext::KLHelper::move_KL_row_to_inverse(coxtypes::CoxNbr x)
{
  coxtypes::CoxNbr xi = inverse(x);
  d_kl.d_klList[x] = d_kl.d_klList[xi];
  d_kl.d_klList[xi] = nullptr;
}


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
    MuRow& row = *d_kl.d_muList[y];
    for (Ulong j = 0; j < row.size(); ++j)
      row[j].x = a[row[j].x];
    std::sort(row.begin(),row.end());
  }

  // permute ranges
  bits::BitMap b(a.size());

  for (coxtypes::CoxNbr x = 0; x < size(); ++x) {
    if (b.getBit(x))
      continue;
    if (a[x] == x) {
      b.setBit(x);
      continue;
    }

    for (coxtypes::CoxNbr y = a[x]; y != x; y = a[y])
    {
      /* back up values for y */
      KLRow* kl_buf = d_kl.d_klList[y];
      auto mu_buf = std::move(d_kl.d_muList[y]);
      /* put values for x in y */
      d_kl.d_klList[y] = d_kl.d_klList[x];
      d_kl.d_muList[y] = std::move(d_kl.d_muList[x]);
      /* store backup values in x */
      d_kl.d_klList[x] = kl_buf;
      d_kl.d_muList[x] = std::move(mu_buf);
      /* set bit*/
      b.setBit(y);
    }  // |for(y)|

    b.setBit(x);
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
  d_kl.d_klList.setSize(n);
  d_kl.d_muList.resize(n);
} // |shrink|




/******** input/output ******************************************************/


/*
  This function prints the status of the context. This is a data structure
  that monitors precisely the computations that have been effected, and
  the memory that has been allocated.
*/
void KLContext::printStatus(FILE* file) const
{
  fprintf(file,"klrows = %lu\n",stats().klrows);
  fprintf(file,"klnodes = %lu\n",stats().klnodes);
  fprintf(file,"klcomputed = %lu\n",stats().klcomputed);
  fprintf(file,"murows = %lu\n",stats().murows);
  fprintf(file,"munodes = %lu\n",stats().munodes);
  fprintf(file,"mucomputed = %lu\n",stats().mucomputed);
  fprintf(file,"muzero = %lu\n",stats().muzero);
}


/****************************************************************************

        Chapter II -- The KLHelper class

  The purpose of the KLHelper class is to hide from the public eye a number
  of helper functions, used in the construction and maintenance of the
  K-L context. This unclutters kl.h quite a bit.

  The following functions are defined :

   - allocExtrRow(const coxtypes::CoxNbr& y) : allocates a row in the extremal list;
   - allocKLRow(const coxtypes::CoxNbr& y) : allocates a row in the K-L list;
   - allocMuRow(const coxtypes::CoxNbr& y) : allocates a row in the mu list;
   - allocMuRow(MuRow& row, const coxtypes::CoxNbr& y) : same, for an external mu-row;
   - allocMuTable() : allocates the full mu-table;
   - allocRowComputation(const coxtypes::CoxNbr& y) : initial allocation for a
     row-computation
   - checkKLRow(const coxtypes::CoxNbr& y) : checks if a K-L row is fully computed;
   - coatomCorrection(const coxtypes::CoxNbr& y, list::List<KLPol>& pol) : subtracts the
     terms ofr coatoms in the mu-correction, for a full row;
   - coatomCorrection(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s,
     list::List<KLPol>& pol, const Ulong& a) : same, for a single polynomial
   - computeMu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y) : computes a mu-coefficient;
   - fillMuRow(MuRow& row, const coxtypes::CoxNbr& y) : fills a row in the mu-table;
   - fillKLPol(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s =
     coxtypes::undef_generator) : fills in one polynomial, using s as descent;
   - fillKLRow(const coxtypes::CoxNbr& y) : fills in one row in the K-L table;
     coxtypes::CoxNbr inverse(const coxtypes::CoxNbr& y) : returns the inverse of y;
   - initWorkspace(const coxtypes::CoxNbr& y, list::List<KLPol>& pol) : another preliminary
     to the computation of a row;
   - inverseMuRow(const coxtypes::CoxNbr& y) : constructs the mu-row for y from that
     of the inverse of y;
   - makeMuRow(const coxtypes::CoxNbr& y);
   - muCorrection(const coxtypes::CoxNbr& y, list::List<KLPol>& pol) : subtracts the non-coatom
     mu-part in the computation of a row;
   - muCorrection(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s,
     list::List<KLPol>& pol, const Ulong& a) : subtracts the non-coatom mu-part,
     for the computation of a single polynomial;
   - muList(const coxtypes::CoxNbr& y) : returns the row for y in muList;
   - prepareRow(const coxtypes::CoxNbr& y, const coxtypes::Generator& s) : a preliminary to the
     computation of a row;
   - readMuRow(const coxtypes::CoxNbr& y) : fills in the mu-row from the K-L row;
   - recursiveMu(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y, const coxtypes::Generator& s) :
     computes mu(x,y) using the general recursive formula;
   - secondTerm(const coxtypes::CoxNbr& y, list::List<KLPol>& pol) : takes care of the
     second term P_{x,ys} in the computation of a full row;
   - writeKLRow(const coxtypes::CoxNbr& y, list::List<KLPol>& pol) : transfers the
     polynomials from pol to klList;
   - writeMuRow(const MuRow& row, const coxtypes::CoxNbr& y) transfers the
     mu-coefficients from row to muList;

 ****************************************************************************/

void KLContext::KLHelper::allocKLRow(const coxtypes::CoxNbr& y)

/*
  This function allocates one row of the kl_list. The row contains one
  entry for each x <= y which is extremal w.r.t. the descent set of y.

  The algorithm is as follows : we extract the interval [e,y] as a
  bitmap, extremalize it by intersecting with the downsets for the
  various generators, and read it into the row.

  Forwards the error MEMORY_WARNING if there is a memory overflow
  and CATCH_MEMORY_OVERFLOW is turned on.
*/

{
  if (!isExtrAllocated(y))
    allocExtrRow(y);

  Ulong n = extrList(y).size();

  d_kl.d_klList[y] = new KLRow(n);
  if (ERRNO)
    return;
  klList(y).setSizeValue(n);

  status().klnodes += n;
  status().klrows++;

  return;

}


/*
  This function allocates one row in the muList. There is one entry for
  each x < y which is extremal w.r.t. y, and has odd length-difference > 1
  with y. As with allocKLRow, this function is not designed for maximal
  efficiency; row allocations for big computations should be handled
  differently.
*/
void KLContext::KLHelper::allocMuRow(const coxtypes::CoxNbr& y)
{
  using EI = FilteredIterator
    <coxtypes::CoxNbr,klsupport::ExtrRow::const_iterator,MuFilter>;
  using BI = FilteredIterator<Ulong,bits::BitMap::Iterator,MuFilter>;

  const schubert::SchubertContext& p = schubert();
  klsupport::ExtrRow e;
  MuFilter f(p,y);

  if (isExtrAllocated(y)) {
    EI first(extrList(y).begin(),extrList(y).end(),f);
    EI last(extrList(y).end(),extrList(y).end(),f);
    e.assign(first,last);
  }
  else {
    bits::BitMap b(size());
    p.extractClosure(b,y);
    if (ERRNO)
      return;
    schubert::select_maxima_for(p,b,p.descent(y));
    BI first(b.begin(),b.end(),f);
    BI last(b.end(),b.end(),f);
    e.assign(first,last);
  }

  coxtypes::Length ly = p.length(y);

  d_kl.d_muList[y].reset(new MuRow(e.size()));
  if (ERRNO)
    goto abort;
  muList(y).setSizeValue(e.size());

  for (Ulong j = 0; j < e.size(); ++j) {
    coxtypes::CoxNbr x = e[j];
    coxtypes::Length lx = p.length(x);
    new(muList(y).ptr()+j) MuData(x,klsupport::undef_klcoeff,(ly-lx-1)/2);
  }

  status().munodes += e.size();
  status().murows++;

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}


/*
  Like allocMuRow(const coxtypes::CoxNbr& y), but does the allocation in row.

  NOTE : this function does not worry about memory overflow. This should
  probably be revised.

*/
void KLContext::KLHelper::allocMuRow(MuRow& row, const coxtypes::CoxNbr& y)
{
  const schubert::SchubertContext& p = schubert();
  klsupport::ExtrRow e;

  if (isExtrAllocated(y)) { /* extremal row is already available */
    e = extrList(y);
  }
  else { /* make it from scratch */
    bits::BitMap b(size());
    p.extractClosure(b,y);
    schubert::select_maxima_for(p,b,p.descent(y));
    schubert::read_bitmap(e,b);
  }

  /* extract elements with odd length difference > 1 */

  Ulong mu_count = 0;
  coxtypes::Length ly = p.length(y);

  for (Ulong j = 0; j < e.size(); ++j) {
    coxtypes::CoxNbr x = e[j];
    coxtypes::Length lx = p.length(x);
    if ((ly-lx)%2 == 0)
      continue;
    if (ly-lx == 1)
      continue;
    e[mu_count] = x;
    mu_count++;
  }

  row.setSize(mu_count);

  for (Ulong j = 0; j < mu_count; ++j) {
    coxtypes::CoxNbr x = e[j];
    coxtypes::Length lx = p.length(x);
    new(row.ptr()+j) MuData(x,klsupport::undef_klcoeff,(ly-lx-1)/2);
  }

  return;
}

void KLContext::KLHelper::allocMuTable()

/*
  This function does the allocation of the full muList, using the closure
  iterator. It is not as satisfactory as it should be, because the rows
  are not obtained in order, and therefore we can not fill them right
  away to eliminate zero entries. So I am currently not too happy with
  the situation. Even if we recover the memory later on, it will be
  highly fragmented.

  So currently this is just an exercise in formal elegance.

*/

{
  typedef FilteredIterator<Ulong,bits::BitMap::Iterator,MuFilter> I;

  const schubert::SchubertContext& p = schubert();

  for (schubert::ClosureIterator cl(p); cl; ++cl) {

    coxtypes::CoxNbr y = cl.current();
    if (inverse(y) < y)
      continue;
    if (isMuAllocated(y))
      continue;

    /* find extremal list */

    bits::BitMap b = cl().bitMap();
    if (ERRNO) {
      printf("error! y = %lu\n",static_cast<Ulong>(y));
      goto abort;
    }

    schubert::select_maxima_for(p,b,p.descent(y));

    MuFilter f(p,y);
    I first(b.begin(),b.end(),f);
    I last(b.end(),b.end(),f);

    klsupport::ExtrRow e(first,last);
    if (ERRNO) {
      goto abort;
    }

    /* transfer to muList */

    coxtypes::Length ly = p.length(y);
    d_kl.d_muList[y].reset(new MuRow(e.size()));
    muList(y).setSizeValue(e.size());

    for (Ulong j = 0; j < e.size(); ++j) {
      coxtypes::CoxNbr x = e[j];
      coxtypes::Length lx = p.length(x);
      MuData mu_data(x,klsupport::undef_klcoeff,(ly-lx-1)/2);
      muList(y).append(mu_data);
      if (ERRNO) {
	goto abort;
      }
    }

    status().murows++;
    status().munodes += e.size();
  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}


/*
  This function does the primary memory allocation for the computation
  of a schubert row. This means that for each initial subword of the
  normal form of y, we check if the corresponding row in klList (or the
  inverse row, if appropriate) is allocated, and if not, do the allocation.

  First, we determine the sequence of left and right shifts that takes
  y to the identity element, as follows : if inverse(y) >= y, we shift
  on the right by last(y); else, we shift on the left by last(inverse(y)),
  and get the inverse of the element which will be used for the recursion
  for inverse(y). So when we reconstruct the interval using this sequence
  of shifts, we know that when we apply a left shift we are actually
  constructing the inverse of the interval we want.

  We can do this in one pass, as we extract the interval [e,y].

  Deals with the possible memory error, and returns ERROR_WARNING in case
  of error.
*/
void KLContext::KLHelper::allocRowComputation(const coxtypes::CoxNbr& y)
{
  klsupport().allocRowComputation(y);

  list::List<coxtypes::Generator> g(0);
  klsupport().standardPath(g,y);

  coxtypes::CoxNbr y1 = 0;

  for (Ulong j = 0; j < g.size(); ++j) {

    coxtypes::Generator s = g[j];
    y1 = schubert().shift(y1,s);
    coxtypes::CoxNbr y2 = klsupport().inverseMin(y1);
    const klsupport::ExtrRow& e = extrList(y2);

    if (!isKLAllocated(y2)) {
      d_kl.d_klList[y2] = new KLRow(e.size());
      if (ERRNO)
	goto abort;
      klList(y2).setSizeValue(extrList(y2).size());
      status().klnodes += extrList(y2).size();
      status().klrows++;
    }

  }

  return;
 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

bool KLContext::KLHelper::checkKLRow(const coxtypes::CoxNbr& d_y)

/*
  This function checks if the row for y (or for inverse(y) if appropriate)
  in klList has been filled.
*/

{
  coxtypes::CoxNbr y = d_y;

  if (inverse(y) < y)
    y = inverse(y);

  if (!isKLAllocated(y)) /* row is not allocated */
    return false;

  KLRow& kl_row = klList(y);

  for (Ulong j = 0; j < kl_row.size(); ++j) {
    if (kl_row[j] == 0)
      return false;
  }

  return true;
}

bool KLContext::KLHelper::checkMuRow(const coxtypes::CoxNbr& y)

/*
  This function checks if the row for y in muList has been filled.
*/

{
  if (!isMuAllocated(y)) /* row is not allocated */
    return false;

  const MuRow& mu_row = muList(y);

  for (Ulong j = 0; j < mu_row.size(); ++j) {
    if (mu_row[j].mu == klsupport::undef_klcoeff)
      return false;
  }

  return true;
}

void KLContext::KLHelper::coatomCorrection(const coxtypes::CoxNbr& y, list::List<KLPol>& pol)

/*
  This function subtracts the coatom correction from the list of polynomials.
  For each coatom z of ys, this is done for all the relevant x-es, as in
  muCorrection.
*/

{
  const schubert::SchubertContext& p = schubert();
  bits::BitMap b(size());
  const klsupport::ExtrRow& e = extrList(y);
  coxtypes::Generator s = last(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);
  const schubert::CoatomList& c = p.hasse(ys);

  for (Ulong j = 0; j < c.size(); ++j) {
    coxtypes::CoxNbr z = c[j];
    if (p.shift(z,s) > z)
      continue;

    p.extractClosure(b,z);
    schubert::select_maxima_for(p,b,p.descent(y));

    Ulong i = 0;
    bits::BitMap::Iterator b_end = b.end();

    for (bits::BitMap::Iterator k = b.begin(); k != b_end; ++k) {
      coxtypes::CoxNbr x = *k;
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
}

void KLContext::KLHelper::coatomCorrection(const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
				 const coxtypes::Generator& s, list::List<KLPol>& pol,
				 const Ulong& a)

/*
  This function subtracts the coatom correction from pol, which at this
  point should contain the value P_{xs,ys}+q.P_{x,ys}-(mu-correction)
  (although we can apply CoatomCorrection and MuCorrection in any order.)

  This means that we subtract the correcting terms corresponding to the
  coatoms z of ys s.t. zs < s; the corresponding subtraction is q.P_{x,z}.

*/

{
  const schubert::SchubertContext& p = schubert();
  coxtypes::CoxNbr ys = p.shift(y,s);  const schubert::CoatomList& c = p.hasse(ys);

  for (Ulong j = 0; j < c.size(); ++j) {

    coxtypes::CoxNbr z = c[j];

    if (p.shift(z,s) > z) /* z is not considered */
      continue;

    if (!p.inOrder(x,z)) /* z is not in [x,ys] */
      continue;

    /* at this point we have to do an actual subtraction */

    const KLPol& p_xz = klPol(x,z);
    if (ERRNO)
      return;
    safeSubtract(pol[a],p_xz,1,1);
    if (ERRNO) {
      Error(ERRNO,this,x,y);
      ERRNO = ERROR_WARNING;
      return;
    }
  }

  return;
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

  bits::Lflags f = p.twoDescent(y);

  if ((p.descent(x)&f) == f) { /* x is super-extremal w.r.t. y */
    return recursiveMu(x,y,last(y));
  }

  coxtypes::Generator s, t;

  /* choose s s.t. LR(ys) not contained in LR(x) */

  for (bits::Lflags f1 = p.descent(y); f1; f1 &= f1-1) {
    coxtypes::Generator u = constants::firstBit(f1);
    coxtypes::CoxNbr yu = p.shift(y,u);
    bits::Lflags fu = p.descent(yu);
    if ((p.descent(x)&fu) != fu) {
      s = u;
      t = constants::firstBit(fu & ~p.descent(x));
      break;
    }
  }

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);

  klsupport::KLCoeff r1 = mu(xs,ys);
  if (ERRNO)
    goto abort;

  /* check if x <= ys */

  if (!p.inOrder(x,ys)) { /* value is found */
    status().mucomputed++;
    if (r1 == 0)
      status().muzero++;
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
	status().mucomputed++;
	if (r1 == r3)
	  status().muzero++;
	return r1-r3;
      }
      else { /* ysts > yst */
	status().mucomputed++;
	if (r1 == 0)
	  status().muzero++;
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
	status().mucomputed++;
	if (r1+r2 == r3)
	  status().muzero++;
	return r1+r2-r3;
      }
      else { /* ysts > yst */
	klsupport::KLCoeff r2 = mu(xt,ys);
	if (ERRNO)
	  goto abort;
	status().mucomputed++;
	if (r1+r2 == 0)
	  status().muzero++;
	return r1+r2;
      }
    }
  }

 abort:
  if (ERRNO != MEMORY_WARNING)
    ERRNO = MU_FAIL;
  return klsupport::undef_klcoeff;
} // |KLContext::KLHelper::computeMu|


/*
  This function performs one of the main tasks in the program, namely
  filling in a Kazhdan-Lusztig polynomial. In this function we are
  not concerned with maximal efficiency, but mostly with computing
  as few polynomials as possible in the recursion formula.

  It is already assumed that x is extremal w.r.t. y.

  There is a tricky point about this function : since it is highly recursive,
  it needs to manage its workspace as a stack; this is buf.

  Sets the error KL_FAIL, and returns the null pointer, in case of failure.

  NOTE : since this function potentially triggers the allocation of many
  rows in the klList and in the muList (not to mention the storage of the
  computed polynomials), it should labor under the constraint of
  CATCH_MEMORY_OVERFLOW.
*/
const KLPol* KLContext::KLHelper::fillKLPol
  (const coxtypes::CoxNbr& x, const coxtypes::CoxNbr& y,
   const coxtypes::Generator& d_s)
{
  static list::List<KLPol> pol(0);

  const schubert::SchubertContext& p = schubert();

  /* check easy cases */

  coxtypes::Length l = p.length(y) - p.length(x);

  if (l < 3) {
    status().klcomputed++;
    return &(one());
  }

  coxtypes::Generator s = d_s;

  /* If d_s is coxtypes::undef_generator, we compute the polynomial using descent by
     last term in normal form */

  if (s == coxtypes::undef_generator)
    s = last(y);

  coxtypes::CoxNbr ys = p.shift(y,s);
  coxtypes::CoxNbr xs = p.shift(x,s);

  /* check if x is comparable to ys */

  if (!p.inOrder(x,ys)) { /* return the answer recursively */
    status().klcomputed++;
    return &klPol(xs,ys);
  }

  /* get workspace */

  CATCH_MEMORY_OVERFLOW = true;

  Ulong a = pol.size();
  pol.setSize(a+1);

  /* initialize the workspace to P_{xs,ys} */

  const KLPol& p_xsys = klPol(xs,ys);
  if (ERRNO)
    goto abort;
  pol[a] = p_xsys;

  /* add q.P_{x,ys} */

  {
    const KLPol& p_xys = klPol(x,ys);
    if (ERRNO)
      goto abort;
    safeAdd(pol[a],p_xys,1);
    if (ERRNO)
      goto abort;
  }

  /* subtract correction terms */

  coatomCorrection(x,y,s,pol,a);
  if (ERRNO)
    goto abort;

  muCorrection(x,y,s,pol,a);
  if (ERRNO)
    goto abort;

  /* find address of polynomial */

  {
    const KLPol& p_xy = *klTree().find(pol[a]);
    if (ERRNO)
      goto abort;

    /* return workspace and exit */

    CATCH_MEMORY_OVERFLOW = false;

    pol.setSize(a);
    status().klcomputed++;
    return &p_xy;
  }

 abort: /* an error occurred */

  CATCH_MEMORY_OVERFLOW = false;
  if (ERRNO != MEMORY_WARNING)
    ERRNO = KL_FAIL;
  return 0;

}

void KLContext::KLHelper::fillKLRow(const coxtypes::CoxNbr& d_y)

/*
  This function fills the rows for y in klList and muList. It is to be called
  if checkKLRow(y) has returned the value false.

  This is one of the big functions in the program, of course. It is typically
  called when an element of the K-L basis of the Hecke algebra is required,
  or when we want to study the singularities of a Schubert variety. We have
  tried to optimize it for speed rather than memory efficiency.

  The outline of the function is as follows. Our descent strategy is to
  always chop off the last element in the normal form of the element under
  consideration, as provided by d_last. This is much better than taking,
  say, the first descent.

   - it is assumed that the rows in klList and muList corresponding to
     the descent path that comes from our descent strategy (and the fact
     that sometimes we move over to inverses) have been allocated; this
     can be assomplished by calling allocRowComputation. This condition will
     then also hold for the elements in the descent path.

   - then, we fill the corresponding rows, starting from the identity;
     the nice thing here is that it will always be guaranteed that
     P_{xs,ys}, P_{x,ys} and the mu(z,ys) are available; this is a recursive
     call to fillKLRow, except that we can be sure that the rows have
     been allocated already.

   - then, we run through muList(ys), which has been filled, and determine
     the z for which mu(z,ys) != 0, and zs < z; these are the terms which
     cause corrections. For each such z, we call fillKLRow(z).

   - we do the same for the coatom list of y;

   - now, we are certain that all the terms we need for all the P_{x,y} in
     klList(y) are available. To avoid having to test the condition x <= z
     in the correction terms, we compute them all at once; again we run
     through the list of z's as above, and for all x <= z which is extremal
     w.r.t. the descent set of y, do the required subtraction. So we will
     work with a large vector of KLCoeff's, and need a table of offsets
     to the corresponding polynomials.

   One nice thing of this setup is that the work on the rows for y starts
   only when all the preparations are finished; it is guaranteed that there
   will not be any recursive calls in the duration of this filling. So,
   contrary to the fillKLPol function, workspace doesn't have to be managed
   as a stack.

  Returns the error ERROR_WARNING in case of failure, after printing an
  error message.
*/

{
  static list::List<KLPol> pol(0);
  const schubert::SchubertContext& p = schubert();
  coxtypes::CoxNbr y = d_y;

  if (y == 0)
    return;

  if (inverse(y) < y) /* fill in the row for inverse(y) */
    y = inverse(y);

  /* recursively fill in rows in the descent path */

  coxtypes::Generator s = last(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);

  if (!checkKLRow(ys)) {
    fillKLRow(ys);
    if (ERRNO)
      goto abort;
  }

  /* make sure the correcting terms are available */

  prepareRow(y,s);
  if (ERRNO)
    goto abort;

  /* prepare workspace; pol holds the row of polynomials;
     initialize workspace with P_{xs,ys} */

  initWorkspace(y,pol);

  /* add q.P_{x,ys} when appropriate */

  secondTerm(y,pol);
  if (ERRNO)
    goto abort;

  /* subtract correcting terms */

  muCorrection(y,pol);
  if (ERRNO)
    goto abort;
  coatomCorrection(y,pol);
  if (ERRNO)
    goto abort;

  /* copy results to row */

  writeKLRow(y,pol);
  if (ERRNO)
    goto abort;

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::fillMuRow(MuRow& row, const coxtypes::CoxNbr& y)

/*
  Fills row with the values for the mu(x,y). It is assumed that the row
  has been correctly allocated (i.e., that it holds only extremal entries
  with odd length difference > 1 w.r.t. y, and at least one entry for each
  non-zero mu.)
*/

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

void KLContext::KLHelper::initWorkspace(const coxtypes::CoxNbr& y, list::List<KLPol>& pol)

/*
  This function sets pol to a row of one polynomial for each x in klList(y),
  and initializes the corresponding pol[j] to klPol(xs,ys).
*/

{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  pol.setSize(e.size());
  if (ERRNO)
    goto abort;

  /* initialize with values P_{xs,ys} */

  {
    coxtypes::Generator s = last(y);
    coxtypes::CoxNbr ys = p.rshift(y,s);

    for (Ulong j = 0; j < e.size(); ++j) {
      coxtypes::CoxNbr xs = p.shift(e[j],s);
      pol[j] = klPol(xs,ys);
      if (ERRNO)
	goto abort;
    }

  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}


/*
  Construct the mu-row for |inverse(y)| from that of |y|, where it is assumed
  the two are distinct and that the row for |inverse(y)| is filled in. We delete
  the old row for |inverse(y)|, since it is faster to reconstruct all its values
  from those of |y| than to only extract the not yet present values from |y|.
*/
void KLContext::KLHelper::inverseMuRow(const coxtypes::CoxNbr& y)
{
  coxtypes::CoxNbr yi = inverse(y);

  if (isMuAllocated(yi))  // then update stats for coming destruction
  {
    const MuRow& m = muList(yi);
    for (Ulong j = 0; j < m.size(); ++j) {
      klsupport::KLCoeff mu = m[j].mu;
      if (mu != klsupport::undef_klcoeff)
	status().mucomputed--;
      if (mu == 0)
	status().muzero--;
    }
    status().munodes -= m.size();
  }

  d_kl.d_muList[yi].reset(new MuRow(muList(y))); // make a copy
  MuRow& m = muList(yi);

  for (Ulong j = 0; j < m.size(); ++j) {
    m[j].x = inverse(m[j].x);
  }

  std::sort(m.begin(),m.end());

  /* update status */

  for (Ulong j = 0; j < m.size(); ++j) {
      klsupport::KLCoeff mu = m[j].mu;
      if (mu != klsupport::undef_klcoeff)
	status().mucomputed++;
      if (mu == 0)
	status().muzero++;
  }

  status().munodes += m.size();
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

void KLContext::KLHelper::muCorrection(const coxtypes::CoxNbr& y, list::List<KLPol>& pol)

/*
  This function subtracts the term P_{x,z}mu(z,ys) from the appropriate
  entry in pol. The idea is to extract [e,z] and run through it, so that
  we avoid calls to inOrder entirely.
*/

{
  const schubert::SchubertContext& p = schubert();
  const klsupport::ExtrRow& e = extrList(y);
  coxtypes::Generator s = last(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);

  const MuRow& mu_row = muList(ys);

  for (Ulong j = 0; j < mu_row.size(); ++j) {
    if (mu_row[j].mu == 0)
      continue;
    coxtypes::CoxNbr z = mu_row[j].x;
    coxtypes::Length h = mu_row[j].height;
    klsupport::KLCoeff mu_value = mu_row[j].mu;
    if (p.shift(z,s) > z)
      continue;

    bits::BitMap b(size());
    p.extractClosure(b,z);
    schubert::select_maxima_for(p,b,p.descent(y));

    Ulong i = 0;
    bits::BitMap::Iterator b_end = b.end();

    for (bits::BitMap::Iterator k = b.begin(); k != b_end; ++k) {
      coxtypes::CoxNbr x = *k;
      while (e[i] < x)
	++i;
      safeSubtract(pol[i],klPol(x,z),mu_value,h+1);
      if (ERRNO) {
	Error(ERRNO,this,x,y);
	ERRNO = ERROR_WARNING;
	return;
      }
    }
  }

  return;
}

void KLContext::KLHelper::muCorrection(const coxtypes::CoxNbr& d_x, const coxtypes::CoxNbr& y,
			     const coxtypes::Generator& d_s, list::List<KLPol>& pol,
			     const Ulong& a)

/*
  This function subtracts from pol the correction terms

               P_{x,z}.mu(z,ys)q^{(l(y)-l(z))/2}

  corresponding to the z in [x,ys] s.t. mu(z,ys) != 0, zs < z, and l(ys)-l(z)
  >= 3; the correction for the coatoms in [x,ys] is dealt with in the
  function coatomCorrection.

  NOTE : as usual, we watch out for memory overflow --- it is assumed that
  CATCH_MEMORY_OVERFLOW is turned on. Overflow in the coefficients cannot
  occur during subtraction; however, we set an error if a negative coefficient
  is found, as this would be a major discovery!

*/

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
    MuRow& m = muList(ys);

    {
      coxtypes::Length ly = p.length(y);

      for (Ulong j = 0; j < m.size(); ++j) {

	coxtypes::CoxNbr z = m[j].x;

	if (p.shift(z,s) > z)
	  continue;
	if (!p.inOrder(x,z))
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

	  safeSubtract(pol[a],p_xz,m[j].mu,h);
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

void KLContext::KLHelper::prepareRow(const coxtypes::CoxNbr& y, const coxtypes::Generator& d_s)

/*
  This function prepares for the filling of row y in klList and
  muList, by making sure that all the correction terms are available.
  It is assumed that fillKLRow(ys) has already been called.

  Deals with the error if it occurs, and sets the error ERROR_WARNING;
*/

{
  const schubert::SchubertContext& p = schubert();

  coxtypes::Generator s = d_s;
  coxtypes::CoxNbr ys = p.shift(y,s);

  if (!checkMuRow(ys)) {
    if (inverse(ys) < ys) {
      readMuRow(inverse(ys));
      inverseMuRow(inverse(ys));
    }
    else
      readMuRow(ys);
  }

  const MuRow& mu_row = muList(ys);

  for (Ulong j = 0; j < mu_row.size(); ++j) {
    if (mu_row[j].mu == 0)
      continue;
    coxtypes::CoxNbr z = mu_row[j].x;
    if (p.shift(z,s) > z)
      continue;
    if (!checkKLRow(z)) {
      allocRowComputation(z);
      if (ERRNO)
	goto abort;
      fillKLRow(z);
      if (ERRNO)
	goto abort;
    }
  }

  {
    const schubert::CoatomList& c = p.hasse(ys);

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      if (p.shift(z,s) > z)
	continue;
      if (!checkKLRow(z)) {
	allocRowComputation(z);
	if (ERRNO)
	  goto abort;
	fillKLRow(z);
	if (ERRNO)
	  goto abort;
      }
    }
  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = ERROR_WARNING;
  return;
}

void KLContext::KLHelper::readMuRow(const coxtypes::CoxNbr& y)

/*
  This function fills the mu-row from the corresponding kl-row. If the
  row has not been allocated yet, it makes sure that no unnecessary
  allocations are made.
*/

{
  const schubert::SchubertContext& p = schubert();
  const KLRow& kl_row = klList(y);
  const klsupport::ExtrRow& e_row = extrList(y);

  if (!isMuAllocated(y)) { /* make row from scratch */
    MuRow mu_buf(0);
    mu_buf.setSizeValue(0);
    coxtypes::Length ly = p.length(y);

    for (Ulong j = 0; j < kl_row.size(); ++j) {

      coxtypes::CoxNbr x = e_row[j];
      coxtypes::Length lx = p.length(x);

      if ((ly-lx)%2 == 0)
	continue;

      if ((ly-lx) == 1)
	continue;

      const KLPol& pol = *kl_row[j];

      polynomials::Degree d = (ly-lx-1)/2;
      if (pol.deg() < d)
	continue;

      MuData m(x,pol[d],d);
      mu_buf.append(m);
      if (ERRNO)
	goto abort;
    }

    d_kl.d_muList[y].reset(new MuRow(mu_buf));
    if (ERRNO)
      goto abort;

    status().munodes += mu_buf.size();
    status().mucomputed += mu_buf.size();
    status().murows++;
  }
  else {

    Ulong i = 0;
    MuRow& mu_row = muList(y);

    for (Ulong j = 0; j < muList(y).size(); ++j) {
      MuData& mu = mu_row[j];
      coxtypes::CoxNbr x = mu.x;
      while(e_row[i] < x)
	++i;
      const KLPol& pol = *kl_row[i];
      coxtypes::Length d = mu.height;
      if (pol.deg() == d)
	mu.mu = pol[d];
      else {
	mu.mu = 0;
	status().muzero++;
      }
      status().mucomputed++;
    }
  }

  return;

 abort:
  Error(ERRNO);
  ERRNO = MEMORY_WARNING;
  return;
}

klsupport::KLCoeff KLContext::KLHelper::recursiveMu(const coxtypes::CoxNbr& d_x, const coxtypes::CoxNbr& y,
			       const coxtypes::Generator& d_s)

/*
  This function computes mu(x,y) using the general recursive formula for
  the descent s. In practice, this will be used only if the descent set
  of x contains the _second_ descent set of y. It is assumed that x < y
  x extremal w.r.t. y, and l(y)-l(x) odd > 1.

  Sets the error MU_FAIL, and returns the value undef_klcoeff, in case of
  failure (this can be due to memory overflow, or to coefficient over- or
  underflow.)
*/

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

  if (!p.inOrder(x,ys)) { /* value is found */
    status().mucomputed++;
    if (r == 0)
      status().muzero++;
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

    MuRow& m = muList(ys);

    for (Ulong j = 0; j < m.size(); ++j) {

      coxtypes::CoxNbr z = m[j].x;

      if (z == x)
	continue;
      if (p.shift(z,s) > z)
	continue;
      if (!p.inOrder(x,z))
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
    const schubert::CoatomList& c = p.hasse(ys);

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      coxtypes::CoxNbr zs = p.shift(z,s);
      if (zs > z)
	continue;
      if (!p.inOrder(x,z))
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

  status().mucomputed++;
  if (r == 0)
    status().muzero++;

  return r;

 abort:
  if (ERRNO != MEMORY_WARNING)
    ERRNO = MU_FAIL;
  return klsupport::undef_klcoeff;
}

void KLContext::KLHelper::secondTerm(const coxtypes::CoxNbr& y, list::List<KLPol>& pol)

/*
  This function takes care of the "second term" q.P_{x,ys} in the
  recursion formula for P_{x,y}. It is assumed that y <= inverse(y)
  and that the descent strategy is via last.

  Here all the memory allocations have been made successfully; the
  only cause of error would be an overflow condition. In that case,
  the error is treated and ERROR_WARNING is set.

  Since we want to avoid all calls to InOrder, the method here is
  to extract [e,ys], extremalize it w.r.t. the descent set of y,
  and run through it and make the correction.
*/

{
  const schubert::SchubertContext& p = schubert();
  bits::BitMap b(0);
  coxtypes::Generator s = last(y);
  coxtypes::CoxNbr ys = p.rshift(y,s);

  /* extract [e,ys] */

  p.extractClosure(b,ys);

  /* extremalize w.r.t. y */

  schubert::select_maxima_for(p,b,p.descent(y));

  /* add to appropriate terms */

  Ulong i = 0;
  bits::BitMap::Iterator b_end = b.end();
  const klsupport::ExtrRow& e = extrList(y);

  for (bits::BitMap::Iterator j = b.begin(); j != b_end; ++j) {
    coxtypes::CoxNbr x = *j;
    while(e[i] < x)
      ++i;
    safeAdd(pol[i],klPol(x,ys),1);
    if (ERRNO) {
      Error(ERRNO,this,x,y);
      ERRNO = ERROR_WARNING;
      return;
    }
  }

  return;
}


/*
  Write the polynomials from the list pol to klList(y); more precisely, it finds
  their adresses in klTree(), and writes those to klList(y). First it has to put
  the true degrees in the pol[j].

  It is assumed that y <= inverse(y).

  The only error that can occur here is memory overflow because of the
  allocation for new polynomials in klTree(). In that case, the error is
  reported, and ERROR_WARNING is set.
*/
void KLContext::KLHelper::writeKLRow
  (const coxtypes::CoxNbr& y, list::List<KLPol>& pol)
{
  KLRow& kl_row = klList(y);

  for (Ulong j = 0; j < kl_row.size(); ++j)
    if (kl_row[j]==nullptr) // dont't overwrite already stored K-L polynomials
    {
      pol[j].snap_degree();
      kl_row[j] = klTree().find(pol[j]);
      if (kl_row[j] == nullptr)
      { /* an error occurred */
	Error(ERRNO);
	ERRNO = ERROR_WARNING;
	return; // give up storing remaining  polynomials
      }
      status().klcomputed++;
    } // |for(j)| and |if(..)|
}

void KLContext::KLHelper::writeMuRow(const MuRow& row, const coxtypes::CoxNbr& y)

/*
  This function copies row to the corresponding row in the mu-list. The idea
  is to copy only the non-zero entries.
*/

{
  /* count non-zero entries */

  Ulong count = 0;

  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu != 0)
      count++;
  }

  MuRow& y_row = muList(y);
  y_row.setSize(count);
  if (ERRNO) {
    Error(ERRNO);
    ERRNO = ERROR_WARNING;
    return;
  }

  count = 0;

  for (Ulong j = 0; j < row.size(); ++j) {
    if (row[j].mu != 0) {
      new(y_row.ptr()+count) MuData(row[j].x,row[j].mu,row[j].height);
      count++;
    }
  }

  status().munodes += count;
  status().murows++;

  return;
}


/****************************************************************************

        Chapter III -- The MuFilter class

  The MuFilter class is a small functor, useful for adapting iterators
  when constructing mu-lists. It filters out elements according to length
  parity and length difference (the intention is to use it with the
  FilteredIterator adaptor.)

 ****************************************************************************/


MuFilter::MuFilter(const schubert::SchubertContext& p, const coxtypes::Length& l)
  :d_p(p), d_l(l)

{}

MuFilter::MuFilter(const schubert::SchubertContext& p, const coxtypes::CoxNbr& y)
  :d_p(p)

{
  d_l = d_p.length(y);
}

MuFilter::~MuFilter()

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

void showKLPol(FILE* file, KLContext& kl, const coxtypes::CoxNbr& d_x, const coxtypes::CoxNbr& d_y,
	       const interface::Interface& I, const coxtypes::Generator& d_s)

/*
  This function prints out the various terms appearing in the computation
  of the Kazhdan-Lusztig polynomial P_{x,y} through the standard recursion
  formula, using the generator s as descent generator.

  It is assumed that x <= y in the Bruhat order, and that s is indeed a
  descent generator for y.
*/

{
  static std::string buf;

  const schubert::SchubertContext& p = kl.schubert();

  coxtypes::CoxNbr x = d_x;
  coxtypes::CoxNbr y = d_y;
  coxtypes::Generator s = d_s;

  const KLPol& pol = kl.klPol(x,y,s);
  if (ERRNO) {
    Error (ERRNO);
    return;
  }

  unsigned long ls = io::LINESIZE;

  buf.clear();

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

  bits::Lflags f = p.descent(y);
  x = p.maximize(x,f);
  if (x > d_x) {
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

    if (!p.inOrder(x,ys)) { /* easy case */
      if (s < kl.rank()) { /* action is on the right */
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
      else { /* action is on the left */
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

    const schubert::CoatomList& c = p.hasse(ys);
    bool coatomcorrection = false;

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      if (p.shift(z,s) > z)
	continue;
      if (!p.inOrder(x,z))
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
      if (!p.inOrder(x,z))
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

void showMu(FILE* file, KLContext& kl, const coxtypes::CoxNbr& d_x, const coxtypes::CoxNbr& y,
	    const interface::Interface& I)

/*
  Maps out the computation of a mu-coefficient. See
  KLContext::KLHelper::computeMu for the algorithm.
*/

{
  static std::string buf;

  const schubert::SchubertContext& p = kl.schubert();

  coxtypes::CoxNbr x = d_x;

  klsupport::KLCoeff r = kl.mu(x,y);
  if (ERRNO)
    goto abort;

  {
    unsigned long ls = io::LINESIZE;

    buf.clear();

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

    bits::Lflags fy = p.descent(y);
    x = p.maximize(x,fy);
    if (x > d_x) {
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

    bits::Lflags f2 = p.twoDescent(y);
    bits::Lflags fx = p.descent(x);

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

void showRecursiveMu(FILE* file, KLContext& kl, const coxtypes::CoxNbr& d_x,
		     const coxtypes::CoxNbr& y, const klsupport::KLCoeff& r, const interface::Interface& I)

/*
  Maps out the computation in the case where the recursion formula has to
  be used.
*/

{
  static std::string buf;

  const schubert::SchubertContext& p = kl.schubert();
  unsigned long ls = io::LINESIZE;

  coxtypes::CoxNbr x = d_x;
  coxtypes::Generator s = kl.last(y);
  coxtypes::Length l = p.length(y) - p.length(x); /* l is odd > 1 */

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);

  if (!p.inOrder(x,ys)) { // mu(x,y) = mu(xs,ys)
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

    const schubert::CoatomList& c = p.hasse(ys);
    bool coatomcorrection = false;

    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxNbr z = c[j];
      if (p.shift(z,s) > z)
	continue;
      if (!p.inOrder(x,z))
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
     if (!p.inOrder(x,z))
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

void showSimpleMu(FILE* file, KLContext& kl, const coxtypes::CoxNbr& x,
		  const coxtypes::CoxNbr& y, const klsupport::KLCoeff& r,
		  const interface::Interface& I)

/*
  Auxiliary to ShowMu. Maps out the computation of a mu-coefficient in the
  case where there is a direct recursion. It is assumed that the computation
  proper has already been tried with success.
*/

{
  static std::string buf;

  const schubert::SchubertContext& p = kl.schubert();
  unsigned long ls = io::LINESIZE;

  coxtypes::Generator s,t;
  for (bits::Lflags f1 = p.descent(y); f1; f1 &= f1-1) {
    coxtypes::Generator u = constants::firstBit(f1);
    coxtypes::CoxNbr yu = p.shift(y,u);
    bits::Lflags fu = p.descent(yu);
    if ((p.descent(x)&fu) != fu) {
      s = u;
      t = constants::firstBit(fu & ~p.descent(x));
      break;
    }
  }

  fprintf(file,"using descent s = %lu and ascent t = %lu\n\n",
	  static_cast<Ulong>(s+1),static_cast<Ulong>(t+1));

  coxtypes::CoxNbr xs = p.shift(x,s);
  coxtypes::CoxNbr ys = p.shift(y,s);
  coxtypes::CoxNbr xt = p.shift(x,t);
  coxtypes::CoxNbr yst = p.shift(ys,t);

  /* consider four cases */

  buf.clear();

  if (p.descent(xt) & constants::lmask[s]) { /* xts < xt */

    buf.append("xs = ");
    p.append(buf,xs,I);
    buf.append("  ys = ");
    p.append(buf,ys,I);

    if (p.descent(yst) & constants::lmask[s]) { /* ysts < yst */
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
    if (p.descent(yst) & constants::lmask[s]) { /* ysts < yst */
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
      if (p.descent(xs) & constants::lmask[t]) {
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
}

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
  This function returns in h the singular locus of cl(X_y). The point is to
  do this while computing as few K-L polynomials as possible.
*/
void genericSingularities(HeckeElt& h, const coxtypes::CoxNbr& y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  bits::BitMap b(p.size());
  bits::BitMap bs(p.size());

  p.extractClosure(b,y);
  schubert::select_maxima_for(p,b,p.descent(y));

  h.clear();

  for (bits::BitMap::ReverseIterator x = b.rbegin(); x != b.rend(); ++x) {
    const KLPol& pol = kl.klPol(*x,y);
    if (ERRNO)
      return;
    if (pol.deg() > 0) { /* remove [e,x[ */
      h.emplace_back(*x,&pol);
      p.extractClosure(bs,*x);
      coxtypes::CoxNbr x1 = *x; /* *x will not be correct anymore after modification */
      b.andnot(bs);
      b.setBit(x1); /* needed to make the decrement correct */
    }
  }

  std::reverse(h.begin(),h.end());

  return;
}

bool isSingular(const HeckeElt& h)

/*
  This function answers yes if one of the polynomials in the row is distinct
  from one, no otherwise. This is equivalent to rational singularity of the
  Schubert variety, when such a geometric context is defined, and the row
  is the extremal row for an element y.

  NOTE : conjecturally, annulation of the term corresponding to the
  extremalization of the origin ensures annulation of all the others.
*/

{
  for (Ulong j = 0; j < h.size(); ++j) {
    const KLPol& pol = h[j].pol();
    if (pol.deg() != 0)
      return true;
  }

  return false;
}

bool isSingular(const KLRow& row)

/*
  This function answers yes if one of the polynomials in the row is distinct
  from one, no otherwise. This is equivalent to rational singularity of the
  Schubert variety, when such a geometric context is defined, and the row
  is the extremal row for an element y.
*/

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

   - ihBetti(h,kl,row) : puts the IH betti numbers of row in h;

 *****************************************************************************/

/*
  This function puts the IH betti numbers of the row in h, in a simple-minded
  approach. There is a serious danger of ineteger overflow here. For now, we set
  the value to undef_coxsize in case of overflow; this should be improved of
  course if it happens too often, using long integers for instance.
*/
void ihBetti(schubert::Homology& h, const coxtypes::CoxNbr& y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  bits::BitMap b(0);
  p.extractClosure(b,y);

  h.setSize(p.length(y)+1);
  h.setZero();
  bits::BitMap::Iterator b_end = b.end();

  for (bits::BitMap::Iterator x = b.begin(); x != b_end; ++x) {
    const KLPol& pol = kl.klPol(*x,y);
    coxtypes::Length d = p.length(*x);
    for (Ulong i = 0; i <= pol.deg(); ++i) {
      if (h[d+i] > coxtypes::COXSIZE_MAX - pol[i])
	h[d+i] = coxtypes::undef_coxnbr;
      else
	h[d+i] += pol[i];
    }
  }

  return;
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
void cBasis(HeckeElt& h, const coxtypes::CoxNbr& y, KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  bits::BitMap b(0);
  p.extractClosure(b,y);

  bits::BitMap::Iterator b_end = b.end();
  h.clear();

  for (bits::BitMap::Iterator x = b.begin(); x != b_end; ++x) {
    const KLPol& pol = kl.klPol(*x,y);
    h.emplace_back(*x,&pol);
  }

  return;
}


/****************************************************************************

        Chapter IX -- Utility functions

  This section defines some utility functions used in this module :

   - allocExtrRow(kl,row,y) : does the allocation of the row for y, in row;
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
  Finds |x| in |row| and returns the address of the corresponding row.
  Returns zero if x is not found. Uses binary search.
*/
MuData* find(MuRow& row, const coxtypes::CoxNbr& x)
{
  Ulong j0 = (Ulong)(-1);

  for (Ulong j1 = row.size(); j1-j0 > 1;) {
    Ulong j = j0 + (j1-j0)/2;
    if (row[j].x == x) // m was found
      return row.ptr()+j;
    if (row[j].x < x)
      j0 = j;
    else
      j1 = j;
  }

  return 0;
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

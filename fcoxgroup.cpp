/*
  This is fcoxgroup.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "fcoxgroup.h"

#include "cells.h"

#define undefined (coxtypes::ParNbr)(coxtypes::PARNBR_MAX + 1)

namespace fcoxgroup {
  using namespace cells;
};

/* local type definitions */

namespace {
  using namespace fcoxgroup;

  class Workspace {
    list::List<coxtypes::ParNbr> d_ica_arr;
    list::List<coxtypes::ParNbr> d_nfca_arr;
    list::List<coxtypes::ParNbr> d_prca_arr;
    list::List<coxtypes::ParNbr> d_rdcw_arr;
  public:
    Workspace();
    void setSize(Ulong n);
    coxtypes::ParNbr *ica_arr() {return d_ica_arr.ptr();}
    coxtypes::ParNbr *nfca_arr() {return d_nfca_arr.ptr();}
    coxtypes::ParNbr *prca_arr() {return d_prca_arr.ptr();}
    coxtypes::ParNbr *rdcw_arr() {return d_rdcw_arr.ptr();}
  };

  coxtypes::CoxSize order(FiniteCoxGroup *W);
  Workspace& workspace();
};

/****************************************************************************

  NOTE : unfinished.

  This file contains code for dealing more efficiently with finite Coxeter
  groups. There are two representations of group elements which are more
  compact than the CoxWord representation.

  The first one, which will work for any rank <= 255, represents the elements
  as arrays of rank coxtypes::ParNbr's. The computations in this representation are
  made through a cascade of small transducers. The drawback of this
  representation is that the size of the automata depends strongly on the
  choice of ordering (as opposed to the minimal root machine, which is
  completely canonical.)

  ...

  Our allocation of workspace avoids having to check the sizes at each
  operation; it is not so clear however if this really makes a difference.

 ****************************************************************************/

/****************************************************************************

      Chapter 0 -- Initialization.

 ****************************************************************************/

namespace {

Workspace::Workspace()
  :d_ica_arr(),
   d_nfca_arr(),
   d_prca_arr(),
   d_rdcw_arr()

{}

void Workspace::setSize(Ulong n)

{
  d_ica_arr.setSize(n);
  d_nfca_arr.setSize(n);
  d_prca_arr.setSize(n);
  d_rdcw_arr.setSize(n);

  return;
}

inline Workspace& workspace()

/*
  Returns its static object, which is initialized on first call.
*/

{
  static Workspace wspace;
  return wspace;
}

};

/****************************************************************************

        Chapter I -- The FiniteCoxGroup class.

  This section defines the FiniteCoxGroup class. The following functions
  are defined :

   - FiniteCoxGroup(x,l) : constructor;
   - assign(a,g) : sets a to the array form of g;
   - duflo() : returns the list of Duflo involutions;
   - inverse(a) : inverses a;
   - isFullContext() : tells if longest_elt is in context;
   - l(r,lr)Cell() : returns the partition in left (right,two-sided) cells;
   - l(r,lr)UneqCell() : returns the partition in left (right,two-sided) cells
     for unequal parameters;
   - l(r)Descent() : returns the partition in left (right) descent classes;
   - l(r)String() : returns the partition in left (right) string classes;
   - l(r)Tau() : returns the partition in left (right) generalized tau classes;
   - length() : returns the length of a;
   - normalForm(g,a) : returns the ShortLex normal form of a in g;
   - normalForm(g,h) : returns the ShortLex normal form of h in g;
   - parseModifier(P) : parses a modifier;
   - prod(a,b) : increments a by b;
   - prod(a,s) : puts in a the result of a.s;
   - prod(a,g) : puts in a the result of a.g;
   - power(a,m) : sets a to the power m;
   - rDescent(a) : returns the right descent set of the coxarr a;
   - reduced(g,a) : puts a reduced expression for a in g;

 ****************************************************************************/

namespace fcoxgroup {

/******** constructor *******************************************************/

FiniteCoxGroup::FiniteCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :CoxGroup(x,l)

/*
  Constructor for FiniteCoxGroup.
*/

{
  d_transducer = new transducer::Transducer(graph());

  workspace().setSize(l);
  for (coxtypes::Rank j = 0; j < rank(); ++j)
    transducer(j)->fill(graph());

  d_longest_coxarr = new(memory::arena()) coxtypes::ParNbr[rank()];

  /* fill longest elements */

  for (transducer::FiltrationTerm* X = transducer(); X; X = X->next())
    d_longest_coxarr[X->rank()-1] = X->size()-1;

  Ulong maxlength = length(d_longest_coxarr);

  new(&d_longest_coxword) coxtypes::CoxWord(maxlength);
  reducedArr(d_longest_coxword,d_longest_coxarr);
  d_longest_coxword.setLength(maxlength);

  d_maxlength = longest_coxword().length();

  d_order = ::order(this);
}

FiniteCoxGroup::~FiniteCoxGroup()

/*
  The only thing that the FiniteCoxGroup destructor has to do explicitly
  is to delete the transducer table, and d_longest.
*/

{
  memory::arena().free(d_longest_coxarr,rank()*sizeof(coxtypes::ParNbr));
  delete d_transducer;

  return;
}

/******** general ***********************************************************/

bool FiniteCoxGroup::isFullContext() const

/*
  Tells if the longest element is in the context. If so, it is necessarily
  the last element fo the context, and can be recognized from its left
  descent set.
*/

{
  coxtypes::CoxNbr x = schubert().size()-1;
  bits::Lflags f = ldescent(x);

  if (f == graph().supp())
    return true;
  else
    return false;
}

/******** operations with arrays ********************************************/

const coxtypes::CoxArr& FiniteCoxGroup::assign(coxtypes::CoxArr& a, const coxtypes::CoxWord& g) const

/*
  This functions returns the array-form of the element of W represented
  by the word g. It returns the result in a.
*/

{
  setZero(a);

  for(coxtypes::Length i = 0; g[i]; ++i)
    prodArr(a,g[i]-1);

  return a;
}


const coxtypes::CoxArr& FiniteCoxGroup::inverseArr(coxtypes::CoxArr& a) const

/*
  Inverse a. This is a "composite-assignment" type function, in consistency
  with our geneeral philosophy that they're the only ones really needed.

  Uses ica_arr() as workspace.
*/

{
  coxtypes::CoxArr b = workspace().ica_arr();

  assign(b,a);
  setZero(a);

  for (const transducer::FiltrationTerm* X = transducer(); X; X = X->next())
    {
      const coxtypes::CoxWord& g = X->np(b[X->rank()-1]);
      Ulong j = g.length();
      while (j) {
	j--;
	prodArr(a,g[j]-1);
      }
    }

  return a;
}


coxtypes::Length FiniteCoxGroup::length(const coxtypes::CoxArr& a) const

/*
  Returns the length of a --- overflow is not checked.
*/

{
  coxtypes::Length c = 0;

  for (const transducer::FiltrationTerm* X = transducer(); X; X = X->next())
    {
      coxtypes::ParNbr x = a[X->rank()-1];
      c += X->length(x);
    }

  return c;
}

const coxtypes::CoxArr& FiniteCoxGroup::powerArr(coxtypes::CoxArr& a, const Ulong& m) const

/*
  Raises a to the m-th power. This can be done very quickly, by squarings
  and multiplications with the original value of a (stored in b), by
  looking at the bit-pattern of m.
*/

{
  static Ulong hi_bit = (Ulong)1 << (BITS(Ulong)-1);
  static list::List<coxtypes::ParNbr> buf(0);

  if (m == 0) {
    setZero(a);
    return a;
  }

  buf.setSize(rank());
  coxtypes::CoxArr b = buf.ptr();
  Ulong p;

  assign(b,a);

  for (p = m; ~p & hi_bit; p <<= 1)  /* shift m up to high powers */
    ;

  for (Ulong j = m >> 1; j; j >>= 1)
    {
      p <<= 1;
      prodArr(a,a);  /* a = a*a */
      if (p & hi_bit)
	prodArr(a,b);  /* a = a*b */
    }

  return a;
}


int FiniteCoxGroup::prodArr(coxtypes::CoxArr& a, const coxtypes::CoxArr& b) const

/*
  Composite assignment operator : increments a by b (i.e., does a *= b).
  The algorithm goes by shifting by the successive pieces of the normal
  form of b, which are directly accessible.

  Uses prca_arr() as workspace;
*/

{
  coxtypes::CoxArr c = workspace().prca_arr();

  assign(c,b);
  int l = 0;

  for (Ulong j = 0; j < rank(); ++j)
    l += prodArr(a,transducer(rank()-1-j)->np(c[j]));

  return l;
}


int FiniteCoxGroup::prodArr(coxtypes::CoxArr& a, coxtypes::Generator s) const

/*
  Transforms the contents of a into a.s.
*/

{
  for (const transducer::FiltrationTerm* X = transducer(); X; X = X->next())
    {
      coxtypes::ParNbr x = a[X->rank()-1];
      coxtypes::ParNbr xs = X->shift(a[X->rank()-1],s);
      if (xs < undefined) {
	a[X->rank()-1] = xs;
	if (xs < x)
	  return -1;
	else
	  return 1;
      }
      s = xs - undefined - 1;
    }

  return 0; // this is unreachable
}


int FiniteCoxGroup::prodArr(coxtypes::CoxArr& a, const coxtypes::CoxWord& g) const

/*
  Shifts a by the whole string g. Returns the increase in length.
*/

{
  int l = 0;

  for (coxtypes::Length j = 0; g[j]; ++j)
    l += prodArr(a,g[j]-1);

  return l;
}


bits::Lflags FiniteCoxGroup::rDescent(const coxtypes::CoxArr& a) const

/*
  Returns the right descent set of a.

  NOTE : makes sense only when the rank is at most MEDRANK_MAX.
*/

{
  bits::Lflags f = 0;

  for (coxtypes::Generator s = 0; s < rank(); s++) /* multiply by s */
    {
      coxtypes::Generator t = s;
      for (const transducer::FiltrationTerm* X = transducer(); X; X = X->next())
	{
	  coxtypes::ParNbr x = a[X->rank()-1];
	  coxtypes::ParNbr xt = X->shift(x,t);
	  if (xt <= undefined) { /* we can decide */
	    if (xt < x)
	    f |= constants::lmask[s];
	    break;
	  }
	  t = xt - undefined - 1;
	}
    }

  return f;
}


const coxtypes::CoxWord& FiniteCoxGroup::reducedArr(coxtypes::CoxWord& g, const coxtypes::CoxArr& a) const

/*
  Returns in g a reduced expression (actually the ShortLex normal form in
  the internal numbering of the generators) of a.

  Here it is assumed that g is large enough to hold the result.
*/

{
  coxtypes::Length p = length(a);
  g[p] = '\0';

  for (const transducer::FiltrationTerm* X = transducer(); X; X = X->next())
    {
      coxtypes::ParNbr x = a[X->rank()-1];
      p -= X->length(x);
      g.setSubWord(X->np(x),p,X->length(x));
    }

  return g;
}

/******** input/output ******************************************************/

void FiniteCoxGroup::modify
  (interface::ParseInterface& P, const interface::Token& tok) const

/*
  Executes the modification indicated by tok, which is assumed to be of
  type modifier_type. It is possible that further characters may have to
  be read from str.

  In the case of a finite coxeter group, three modifies are allowed :
  *, ! and ^
*/

{
  if (interface::isLongest(tok)) {
    CoxGroup::prod(P.c,d_longest_coxword);
  }

  if (interface::isInverse(tok)) {
    CoxGroup::inverse(P.c);
  }

  if (interface::isPower(tok)) {
    Ulong m = readCoxNbr(P,ULONG_MAX);
    CoxGroup::power(P.c,m);
  }
}

bool FiniteCoxGroup::parseModifier(interface::ParseInterface& P) const

/*
  This function parses a modifier from P.str at P.offset, and acts upon
  it accordingly : in case of success, it applies the modifier to P.c,
  and advances the offset. In the case of a finite group, multiplication
  by the longest element is allowed.
*/

{
  interface::Token tok = 0;
  const interface::Interface& I = interface();

  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (not interface::isModifier(tok))
    return false;

  P.offset += p;
  modify(P,tok);

  return true;
}

/******** kazhdan-lusztig cells *********************************************/


/*
  This function returns the list of Duflo involutions in the group, in the
  order in which left cells are listed in lCell : duflo[j] is the Duflo
  involution in the j-th cell of d_lcell.

  The algorithm is as follows. We partition the involutions in the group
  according to the left cell partition. Then the Duflo involution in the
  cell is the unique involution for which l(w)-2d(w) is minimal, where
  d(w) is the degree of the Kazhdan-Lusztig polynomial P_{1,w}.

  NOTE : as for the l(r,lr)cell partitions, the list is filled upon the
  first call.
*/
const list::List<coxtypes::CoxNbr>& FiniteCoxGroup::duflo()
{
  if (d_duflo.size() == 0) { /* find duflo involutions */

    kl::KLContext& kl = *d_kl;
    const schubert::SchubertContext& p = kl.schubert();

    bits::SubSet q(0);

    /* make sure left cell partition is available */

    lCell();

    /* load involutions in q */

    q.bitMap().assign(kl.involution());
    q.readBitMap();

    /* partition involutions by left cells */

    bits::Partition pi(q.size());
    for (Ulong j = 0; j < q.size(); ++j) {
      pi[j] = d_lcell[q[j]];
    }
    pi.setClassCount(d_lcell.classCount());

    /* find Duflo involution in each cell */

    for (bits::PartitionIterator i(pi); i; ++i) {
      const list::List<Ulong>& c = i();
      if (c.size() == 1) { /* cell has single involution */
	d_duflo.append(q[c[0]]);
	continue;
      }
      coxtypes::Length m = d_maxlength;
      coxtypes::CoxNbr d = c[0];
      for (Ulong j = 0; j < c.size(); ++j) {
	coxtypes::CoxNbr x = q[c[j]]; /* current involution */
	const kl::KLPol& pol = kl.klPol(0,x);
	coxtypes::Length m1 = p.length(x) - 2*pol.deg();
	if (m1 < m) {
	  m = m1;
	  d = x;
	}
      }
      d_duflo.append(d);
    }

  }

  return d_duflo;
}

const bits::Partition& FiniteCoxGroup::lCell()

/*
  Returns the partition into left cells, making it from the right cell
  partitition.
*/

{
  if (d_lcell.classCount()) /* partition was already computed */
    return d_lcell;

  const bits::Partition& r = rCell();
  d_lcell.setSize(r.size());
  d_lcell.setClassCount(r.classCount());

  for (coxtypes::CoxNbr x = 0; x < r.size(); ++x) {
    d_lcell[x] = r(CoxGroup::inverse(x));
  }

  d_lcell.normalize();
  return d_lcell;
}

const bits::Partition& FiniteCoxGroup::lrCell()

/*
  Similar to rCell, but for two-sided cells.
*/

{
  if (d_lrcell.classCount()) /* partition was already computed */
    return d_lrcell;

  if (!isFullContext()) {
    fullContext();
    if (error::ERRNO)
      goto abort;
    kl().fillMu();
    if (error::ERRNO)
      goto abort;
  }

  if (d_lrcell.size() == 0) /* size is either zero or group order */
    cells::lrCells(d_lrcell,kl());

  return d_lrcell;

 abort:
  error::Error(error::ERRNO);
  return d_lrcell;
}

const bits::Partition& FiniteCoxGroup::lrUneqCell()

/*
  Similar to lCell, but for two-sided cells.
*/

{
  if (d_lruneqcell.classCount()) /* partition was already computed */
    return d_lruneqcell;

  if (!isFullContext()) {
    fullContext();
    if (error::ERRNO)
      goto abort;
    uneqkl().fillMu();
    if (error::ERRNO)
      goto abort;
  }

  {
    wgraph::OrientedGraph X(0);
    lrGraph(X,uneqkl());
    X.cells(d_lruneqcell);
  }

  return d_lruneqcell;

 abort:
  error::Error(error::ERRNO);
  return d_lruneqcell;
}

const bits::Partition& FiniteCoxGroup::lUneqCell()

/*
  Returns the partition in left cells for unequal parameters. The partition
  is gotten from the right one, by inversing.
*/

{
  if (d_luneqcell.classCount()) /* partition was already computed */
    return d_luneqcell;

  const bits::Partition& r = rUneqCell();
  d_luneqcell.setSize(r.size());
  d_luneqcell.setClassCount(r.classCount());

  for (coxtypes::CoxNbr x = 0; x < r.size(); ++x) {
    d_luneqcell[x] = r(CoxGroup::inverse(x));
  }

  d_luneqcell.normalize();
  return d_luneqcell;
}

const bits::Partition& FiniteCoxGroup::rCell()

/*
  This function returns the partition of the group in right cells.

  NOTE : to be on the safe side, we allow this function to respond only
  for the full group context. If the context is not full, it extends it
  first to the full group.

  NOTE : because this is a potentially very expensive operation, the
  partition is computed on request.

  NOTE : since it is not clear that the ordering in which rCells constructs
  the cells is meaningful, we normalize the partition, so that it can be
  guaranteed to always have the same meaning.
*/

{
  if (d_rcell.classCount()) /* partition was already computed */
    return d_rcell;

  if (!isFullContext()) {
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  kl().fillMu();
  if (error::ERRNO)
    goto abort;

  cells::rCells(d_rcell,kl());
  d_rcell.normalize();

  return d_rcell;

 abort:
  error::Error(error::ERRNO);
  return d_rcell;

}

const bits::Partition& FiniteCoxGroup::rUneqCell()

/*
  This function returns the partition of the group in right cells for unequal
  parameters.

  NOTE : to be on the safe side, we allow this function to respond only
  for the full group context. If the context is not full, it extends it
  first to the full group.

  NOTE : because this is a potentially very expensive operation, the
  partition is computed on request.
*/

{
  if (d_runeqcell.classCount()) /* partition was already computed */
    return d_runeqcell;

  if (!isFullContext()) {
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  d_uneqkl->fillMu();
  if (error::ERRNO)
    goto abort;

  {
    wgraph::OrientedGraph Y(0);
    rGraph(Y,uneqkl());
    Y.cells(d_runeqcell);
    d_runeqcell.normalize();
  }

  return d_runeqcell;

 abort:
  error::Error(error::ERRNO);
  return d_runeqcell;
}

const bits::Partition& FiniteCoxGroup::lDescent()

/*
  Returns the partition of the group in left descent classes, where two
  elements are equivalent iff they have the same left descent set (this is
  the non-generalized tau-invariant of Vogan.)

  It is known that this partition is coarser than the one by right cells.
*/

{
  if (d_ldescent.classCount()) /* partition was already computed */
    return d_ldescent;

  if (!isFullContext()) {
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  d_ldescent.setSize(order());

  for (coxtypes::CoxNbr x = 0; x < order(); ++x)
    d_ldescent[x] = ldescent(x);

  d_ldescent.setClassCount(1<<rank());

  return d_ldescent;

 abort:
  error::Error(error::ERRNO);
  return d_ldescent;
}

const bits::Partition& FiniteCoxGroup::rDescent()

/*
  Returns the partition of the group in right descent classes, where two
  elements are equivalent iff they have the same right descent set (this is
  the non-generalized tau-invariant of Vogan.)

  it is known that this partition is coarser than the one by left cells.
*/

{
  if (d_rdescent.classCount()) /* partition was already computed */
    return d_rdescent;

  if (!isFullContext()) {
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  d_rdescent.setSize(order());

  for (coxtypes::CoxNbr x = 0; x < order(); ++x)
    d_rdescent[x] = rdescent(x);

  d_rdescent.setClassCount(1<<rank());

  return d_rdescent;

 abort:
  error::Error(error::ERRNO);
  return d_rdescent;
}

const bits::Partition& FiniteCoxGroup::lString()

/*
  Returns the partition of the group in "left string classes" : the smallest
  subsets C with the property that for each x in C, and for each pair of
  non-commuting generators {s,t} such that L(x) contains exactly one of s,t,
  say s, and the order m(s,t) of st is finite, the whole {s,t}-string

         ... < tsx < sx < x < tx < stx ...

  (the number of elements is m-1) passing through x is contained in C (see
  the INTRO file for more details.) It is known (and easy to see) that this
  partition is finer than the one by left cells.
*/

{
  if (d_lstring.classCount()) /* partition was already computed */
    return d_lstring;

  if (!isFullContext()) { /* x is not the longest element */
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  lStringEquiv(d_lstring,schubert());

  return d_lstring;

 abort:
  error::Error(error::ERRNO);
  return d_lstring;
}

const bits::Partition& FiniteCoxGroup::rString()

/*
  The same as lString(), but on the right. Finer than the partition by right
  cells.
*/

{
  if (d_rstring.classCount()) /* partition was already computed */
    return d_rstring;

  if (!isFullContext()) { /* x is not the longest element */
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  rStringEquiv(d_rstring,schubert());

  return d_rstring;

 abort:
  error::Error(error::ERRNO);
  return d_rstring;
}

const bits::Partition& FiniteCoxGroup::lTau()

/*
  This is the left generalized-tau partition; in brief, the left-descent
  partition "stabilized" under left *-operations (see the INTRO file for
  more details.) It is known to be coarser than the partition by right cells.

  Is gotten from the corresponding right tau partition.
*/

{
  if (d_ltau.classCount()) /* partition was already computed */
    return d_ltau;

  const bits::Partition& r = rTau();
  d_ltau.setSize(r.size());
  d_ltau.setClassCount(r.classCount());

  for (coxtypes::CoxNbr x = 0; x < r.size(); ++x) {
    d_ltau[x] = r(CoxGroup::inverse(x));
  }

  d_ltau.normalize();
  return d_ltau;
}

const bits::Partition& FiniteCoxGroup::rTau()

/*
  Like lTau, but on the right. Coarser than the partition by left cells.
*/

{
  if (d_rtau.classCount()) /* partition was already computed */
    return d_rtau;

  if (!isFullContext()) { /* x is not the longest element */
    fullContext();
    if (error::ERRNO)
      goto abort;
  }

  rGeneralizedTau(d_rtau,schubert());
  d_rtau.normalize();

  return d_rtau;

 abort:
  error::Error(error::ERRNO);
  return d_rtau;
}

};

/****************************************************************************

        Chapter II -- Derived classes.

  This section contains the constructors for the derived classes of
  FiniteCoxGroup appearing in this program.

  NOTE : unfinished.

 ****************************************************************************/

namespace fcoxgroup {

FiniteBigRankCoxGroup::FiniteBigRankCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteCoxGroup(x,l)

/*
  Constructor for FiniteBigRankCoxGroup.
*/

{}

FiniteBigRankCoxGroup::~FiniteBigRankCoxGroup()

/*
  Virtual destructor for FiniteBigRankCoxGroup. Currently, nothing has to
  be done.
*/

{}

GeneralFBRCoxGroup::GeneralFBRCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteBigRankCoxGroup(x,l)

{}

GeneralFBRCoxGroup::~GeneralFBRCoxGroup()

/*
  Non-virtual destructor (leaf class). Currently, nothing has to be done.
*/

{}

FiniteMedRankCoxGroup::FiniteMedRankCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteCoxGroup(x,l)

/*
  Constructor for FiniteMedRankCoxGroup.
*/

{
  mintable().fill(graph());

  /* an error is set here in case of failure */

  return;
}

FiniteMedRankCoxGroup::~FiniteMedRankCoxGroup()

/*
  Virtual destructor for FiniteMedRankCoxGroup. The destruction of the
  mintable is the job of the CoxGroup destructor.
*/

{}

GeneralFMRCoxGroup::GeneralFMRCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteMedRankCoxGroup(x,l)

{}

GeneralFMRCoxGroup::~GeneralFMRCoxGroup()

/*
  Non-virtual destructor (leaf class). Currently, nothing has to be done.
*/

{}

FiniteSmallRankCoxGroup::FiniteSmallRankCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteMedRankCoxGroup(x,l)

/*
  Constructor for FiniteSmallRankCoxGroup.
*/

{}

FiniteSmallRankCoxGroup::~FiniteSmallRankCoxGroup()

/*
  Virtual destructor for FiniteSmallRankCoxGroup. Currently, nothing has to
  be done.
*/

{}

GeneralFSRCoxGroup::GeneralFSRCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteSmallRankCoxGroup(x,l)

{}

GeneralFSRCoxGroup::~GeneralFSRCoxGroup()

/*
  Non-virtual destructor (leaf class). Currently, nothing has to be done.
*/

{}

/********  The SmallCoxGroup class ******************************************/

SmallCoxGroup::SmallCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :FiniteSmallRankCoxGroup(x,l)

{}

SmallCoxGroup::~SmallCoxGroup()

/*
  Virtual destructor for the SmallCoxGroup class. Currently, nothing has
  to be done.
*/

{}

GeneralSCoxGroup::GeneralSCoxGroup(const type::Type& x, const coxtypes::Rank& l)
  :SmallCoxGroup(x,l)

{}

GeneralSCoxGroup::~GeneralSCoxGroup()

/*
  Non-virtual destructor (leaf class). Currently, nothing has to be done.
*/

{}

const coxtypes::CoxArr& SmallCoxGroup::assign(coxtypes::CoxArr& a, const DenseArray& d_x) const

/*
  Unpacks the DenseArray x into a.
*/

{
  const transducer::Transducer& T = d_transducer[0];
  DenseArray x = d_x;

  for (Ulong j = 0; j < rank(); ++j) {
    const transducer::FiltrationTerm& X = T.transducer(rank()-1-j)[0];
    a[j] = x%X.size();
    x /= X.size();
  }

  return a;
}

void SmallCoxGroup::assign(DenseArray& x, const coxtypes::CoxArr& a) const

/*
  Packs the array a into x.
*/

{
  x = 0;

  for (const transducer::FiltrationTerm* X = transducer(); X; X = X->next()) {
    x *= X->size();
    x += a[X->rank()-1];
  }
}

bool SmallCoxGroup::parseDenseArray(interface::ParseInterface& P) const

/*
  Tries to parse a DenseArray from P. This is a '#' character, followed
  by an integer which has to lie in the range [0,N[, where N is the size
  of the group.
*/

{
  const interface::Interface& I = interface();

  interface::Token tok = 0;
  Ulong p = I.getToken(P,tok);

  if (p == 0)
    return false;

  if (not interface::isDenseArray(tok))
    return false;

  // if we get to this point, we must read a valid integer

  P.offset += p;
  coxtypes::CoxNbr x = interface::readCoxNbr(P,d_order);

  if (x == coxtypes::undef_coxnbr) { //error
    P.offset -= p;
    error::Error(error::DENSEARRAY_OVERFLOW,d_order);
    error::ERRNO = error::PARSE_ERROR;
  }
  else { // x is valid
    coxtypes::CoxWord g(0);
    prodD(g,x);
    CoxGroup::prod(P.c,g);
  }

  return true;
}

bool SmallCoxGroup::parseGroupElement(interface::ParseInterface& P) const

/*
  This is the parseGroupElement function for the SmallCoxGroup type. In
  this class, we have one additional representation of elements, viz. the
  densearray representation. This means that an element that would be
  represented by the array [x_1, ... ,x_n] is represented by the number
  w = x_1+x_2*a_1+ ... +x_n*a_{n-1}, where a_j is the size of the j'th
  subgroup in the filtration. This will give a bijective correspondence
  between group elements and numbers in the range [0,N-1], where N is
  the size of the group.
*/

{
  Ulong r = P.offset;

  if (parseContextNumber(P)) { // next token is a context number
    if (error::ERRNO) // parse error
      return true;
    else
      goto modify;
  }

  if (parseDenseArray(P)) { // next token is a dense array
    if (error::ERRNO) // parse error
      return true;
    else
      goto modify;
  }

  // if we get to this point, we have to read a coxtypes::CoxWord

  interface().parseCoxWord(P,mintable());

  if (error::ERRNO) { // no coxtypes::CoxWord could be parsed
    if (P.offset == r) { // nothing was parsed
      error::ERRNO = 0;
      return false;
    }
    else // parse error
      return true;
  }

 modify:

  // if we get to this point, a group element was successfully read

  while (parseModifier(P)) {
    if (error::ERRNO)
      return true;
  }

  // flush the current group element

  prod(P.a[P.nestlevel],P.c);
  P.c.reset();

  if (P.offset == r) // nothing was read; c is unchanged
    return false;
  else
    return true;
}

int SmallCoxGroup::prodD(coxtypes::CoxWord& g, const DenseArray& d_x) const

/*
  Does the multiplication of g by x, by recovering the normal pieces of x.
  returns the length increase.
*/

{
  const transducer::Transducer& T = d_transducer[0];

  DenseArray x = d_x;
  int l = 0;

  for (Ulong j = 0; j < rank(); ++j) {
    const transducer::FiltrationTerm& X = T.transducer(rank()-1-j)[0];
    coxtypes::ParNbr c = x%X.size();
    l += CoxGroup::prod(g,X.np(c));
    x /= X.size();
  }

  return l;
}

int SmallCoxGroup::prodD(DenseArray& x, const coxtypes::CoxWord& g) const

/*
  Does the multiplication of x by g.
*/

{
  static list::List<coxtypes::ParNbr> al(0);

  al.setSize(rank());
  coxtypes::CoxArr a = al.ptr();
  assign(a,x);
  int l = prodArr(a,g);
  assign(x,a);

  return l;
}

};

/****************************************************************************

        Chapter III -- Auxiliary functions.

  This section contains some auxiliary functions for the construction of
  finite Coxeter groups. The following functions are defined :

   - fillLongest(W) : fills in the longest element;
   - order(W) : returns the order of the group;

 ****************************************************************************/

namespace {


/*
  Initializes the following constants :

    - W->longest_coxarr : array form of the longest element in W;
    - W->longest_coxword : string form of the longest element in W;

*/


coxtypes::CoxSize order(FiniteCoxGroup *W)

/*
  This function fills in the order of W. It sets order to the order
  of W if the order is <= COXSIZE_MAX, sets it to 0 (for the time being)
  otherwise.

  Assumes that the order of subgroup has been filled in already.
*/

{
  coxtypes::CoxSize order = 1;

  for (const transducer::FiltrationTerm* X = W->transducer(); X; X = X->next())
    {
      if (X->size() > coxtypes::COXSIZE_MAX/order) /* overflow */
	return 0;
      order *= X->size();
    }

  return order;
}

};

/****************************************************************************

      Chapter IV -- Types and sizes.

  This section regroups some auxiliary functions for type recognition
  and size computations.

  The functions provided are :

  - isFiniteType(W) : recognizes a finite group --- defines the notion of
    a finite group in this program.
  - maxSmallRank : the maximum rank for a SmallCoxGroup on the current
    machine;

 ****************************************************************************/

bool fcoxgroup::isFiniteType(coxgroup::CoxGroup *W)

/*
  Recognizes the type of a finite group. Non-irreducible types are
  allowed; they are words in the irreducible types. This function
  defines the class of groups that will be treated as finite groups
  in this program; the i/o setup is flexible enough that there is
  no reason that a finite group should be entered otherwise.
*/

{
  return isFiniteType(W->type());
}


coxtypes::Rank fcoxgroup::maxSmallRank(const type::Type& x)

/*
  Returns the smallest rank for which a CoxNbr holds the given element.
  It is assumed that x is one of the finite types A-I.
*/

{
  coxtypes::Rank l;
  unsigned long c;

  switch(x[0])
    {
    case 'A':
      c = 1;
      for (l = 1; l < coxtypes::RANK_MAX; l++) {
	if (c > coxtypes::COXNBR_MAX/(l+1))  /* l is too big */
	  return l-1;
	c *= (l+1);
      }
      return l;
    case 'B':
    case 'C':
      c = 2;
      for (l = 2; l < coxtypes::RANK_MAX; l++) {
	if (c > coxtypes::COXNBR_MAX/2*l)  /* l is too big */
	  return l-1;
	c *= 2*l;
      }
      return l;
    case 'D':
      c = 4;
      for (l = 3; l < coxtypes::RANK_MAX; l++) {
	if (c > coxtypes::COXNBR_MAX/2*l)  /* l is too big */
	  return l-1;
	c *= 2*l;
      }
      return l;
	return l;
    case 'E':
      if (coxtypes::COXNBR_MAX < 2903040)
	return 6;
      else if (coxtypes::COXNBR_MAX < 696729600)
	return 7;
      else
	return 8;
    case 'F':
      return 4;
    case 'G':
      return 2;
    case 'H':
      return 4;
    case 'I':
      return 2;
    default: // unreachable
      return 0;
    };
}

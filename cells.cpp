/*
  This is cells.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "cells.h"
#include "wgraph.h"

#include "stack.h"

namespace cells {
  using namespace klsupport;
  using namespace stack;
};

/****************************************************************************

  This module contains code for the computation of Kazhdan-Lustig cells for
  the group (say in the finite case; in the infinite case the best we can
  hope for is to get pieces of cells.)

  The main tool for this seems to be the systematic use of *-operations;
  this is particularly effective in the simply laced cases, but very useful
  also in the other ones. These operations are recalled and implemented
  in the Schubert module.

  Let us recall the results that we use here, all contained in the original
  Inventiones paper of Kazhdan and Lusztig. Each left cell is stable under
  left *-operations (insofar as they are defined); in other terms for each
  x such that *x (w.r.t. some simple edge {s,t}) is defined, *x is left-
  equivalent to x. So left cells are unions of left star-orbits (in type A,
  it is even so that the left cells are the left star-orbits.) Furthermore,
  the domain of every right star-operation is a union of right equi-descent
  classes, and therefore a union of left cells; and the operation takes left
  cells to left cells isomorphically as W-graphs. In particular also, each
  two-sided cell is a union of two-sided star orbits.

  So our first goal should be to try to get at the functions defining the
  partitions in left, right and two-sided cells, computing as few
  mu-coefficients as possible; then we will want to determine the full
  W-graph structure on the cells, classify them up to isomorphism, and
  have functions such as "get the (left, right, 2-sided) cell of an element."

  More sophisticated data for cells (distinguished involutions, a-functions
  ...) will have to wait a little bit more.

 ****************************************************************************/

namespace {
  using namespace cells;

  typedef list::List<coxtypes::CoxNbr> CoxList;
};

/****************************************************************************

        Chapter I -- Partitions

  This section contains functions that define partitions of (subsets of)
  the context, useful for cell-computations. The following functions are
  defined :

    - lCells(pi,p), rCells(pi,kl), lrCells(pi,kl) : partition of the context
      into left (right, two-sided) kazhdan-lusztig cells;
    - lDescentPartition(pi,p), rDescentPartititon(pi,p) : partition of p
      according to left (right) descent sets;
    - lStringEquiv(pi,p), rStringEquiv(pi,p) : partition of p according
      to weak Bruhat equivalences;
    - lGeneralizedTau(pi,p), rGeneralizedTau(pi,p) : left (right) descent
      partition, stabilized under star operations;

 ****************************************************************************/

namespace cells {


/*
  This function puts in pi the partition of p into left cells --- in the case
  of an incomplete context, they will be the cells defined by the links in
  the graph.

  The idea will be to minimize K-L computations. We proceed as follows. First,
  we determine the generalized-tau partition of the context. Then, we look
  at the star-orbits among the tau-classes, and decompose one representative
  of each into cells; we propagate the cells using star-operations.
*/

void lCells(bits::Partition& pi, kl::KLContext& kl)
{
  static bits::SubSet q(0);
  static bits::SubSet a(0);
  static bits::Partition qcells(0);
  static list::List<Ulong> cell_count(0);
  static list::List<Ulong> qcell_count(0);
  static wgraph::OrientedGraph P(0);
  static Fifo<Ulong> orbit;

  schubert::SchubertContext& p = kl.schubert();
  q.setBitMapSize(p.size());
  a.setBitMapSize(p.size());
  a.reset();
  cell_count.setSize(0);

  rGeneralizedTau(pi,p);

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x) {

    /* a holds the elements already processed */

    if (a.isMember(x))
      continue;

    /* put the next generalized-tau class in q */

    q.reset();
    pi.writeClass(q.bitMap(),pi(x));
    q.readBitMap();

    /* put cell-partition of q in qcells */

    wgraph::WGraph X = W_graph<'l'>(q,kl);
    X.graph().cells(qcells,&P);

    /* the fifo-list orbit is used to traverse the *-orbit of the first
       element of the current generalized-tau class */

    orbit.push(a.size());
    qcell_count.setSize(0);

    /* get class counts and mark off cells in q */

    for (bits::PartitionIterator i(qcells); i; ++i) {
      const bits::Set& c = i();
      qcell_count.append(c.size());
      cell_count.append(c.size());
      for (Ulong j = 0; j < c.size(); ++j)
	a.add(q[c[j]]);
    }

    /* propagate cells with star-operations; the idea is that star operations
       act on the level of generalized-tau classes, so each element in a given
       generalized-tau class accepts the same language. */

    while (orbit.size()) {

      Ulong c = orbit.pop();
      coxtypes::CoxNbr z = a[c];

      for (coxtypes::StarOp j = 0; j < p.nStarOps(); ++j) {

	coxtypes::CoxNbr zj = p.star(z,j);

	if (zj == coxtypes::undef_coxnbr)
	  continue;
	if (a.isMember(zj))
	  continue;

	/* mark off orbit */

	orbit.push(a.size());

	for (Ulong i = 0; i < q.size(); ++i) {
	  coxtypes::CoxNbr y = a[c+i];
	  coxtypes::CoxNbr yj = p.star(y,j);
	  a.add(yj);
	}

	for (Ulong i = 0; i < qcell_count.size(); ++i) {
	  cell_count.append(qcell_count[i]);
	}

      }
    }
  }

  /* write down the partition */

  Ulong c = 0;

  for (Ulong j = 0; j < cell_count.size(); ++j) {
    for (Ulong i = 0; i < cell_count[j]; ++i) {
      pi[a[c+i]] = j;
    }
    c += cell_count[j];
  }

  pi.setClassCount(cell_count.size());

  return;
}


/*
  Same as lCells, but does the partition into right cells.
*/
void rCells(bits::Partition& pi, kl::KLContext& kl)
{
  static bits::SubSet q(0);
  static bits::SubSet a(0);
  static bits::Partition qcells(0);
  static list::List<Ulong> cell_count(0);
  static list::List<Ulong> qcell_count(0);
  static wgraph::OrientedGraph P(0);
  static Fifo<Ulong> orbit;
  static bits::Permutation v(0);

  schubert::SchubertContext& p = kl.schubert();
  q.setBitMapSize(p.size());
  a.setBitMapSize(p.size());
  a.reset();
  cell_count.setSize(0);

  lGeneralizedTau(pi,p);

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x) {

    if (a.isMember(x))
      continue;

    /* get the next generalized-tau class */

    q.reset();
    pi.writeClass(q.bitMap(),pi(x));
    q.readBitMap();

    /* find cells in class */

    wgraph::WGraph X = W_graph<'r'>(q,kl);
    X.graph().cells(qcells,&P);

    /* the fifo-list orbit is used to traverse the *-orbit of the first
       element of the current generalized-tau class */

    orbit.push(a.size());
    qcell_count.setSize(0);

    /* get class counts and mark off cells in q */

    for (bits::PartitionIterator i(qcells); i; ++i) {
      const bits::Set& c = i();
      qcell_count.append(c.size());
      cell_count.append(c.size());
      for (Ulong j = 0; j < c.size(); ++j)
	a.add(q[c[j]]);
    }

    /* propagate cells with star-operations */

    while (orbit.size()) {

      Ulong c = orbit.pop();
      coxtypes::CoxNbr z = a[c];

      for (coxtypes::StarOp j = p.nStarOps(); j < 2*p.nStarOps(); ++j) {

	coxtypes::CoxNbr zj = p.star(z,j);

	if (zj == coxtypes::undef_coxnbr)
	  continue;
	if (a.isMember(zj))
	  continue;

	/* mark off orbit */

	orbit.push(a.size());

	for (Ulong i = 0; i < q.size(); ++i) {
	  coxtypes::CoxNbr y = a[c+i];
	  coxtypes::CoxNbr yj = p.star(y,j);
	  a.add(yj);
	}

	for (Ulong i = 0; i < qcell_count.size(); ++i) {
	  cell_count.append(qcell_count[i]);
	}

      }
    }
  }

  /* write down the partition */

  Ulong c = 0;

  for (Ulong j = 0; j < cell_count.size(); ++j) {
    for (Ulong i = 0; i < cell_count[j]; ++i) {
      pi[a[c+i]] = j;
    }
    c += cell_count[j];
  }

  pi.setClassCount(cell_count.size());

  return;
}


/*
  This function computes the two-sided cells in the context. There are
  certainly better ways to do this, but I'm afraid I don't know enough
  to do it other than by filling in all the mu's ...
*/
void lrCells(bits::Partition& pi, kl::KLContext& kl)
{
  kl.fillMu();

  wgraph::WGraph X = cells::W_graph<'b'>(kl);
  X.graph().cells(pi);
}


/*
  This function writes in pi the partition of p according to the left
  descent sets.
*/
void lDescentPartition(bits::Partition& pi, const schubert::SchubertContext& p)
{
  static list::List<GenSet> d(0); /* holds the appearing descent sets */

  pi.setSize(p.size());
  d.setSize(0);

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
    insert(d,p.ldescent(x));

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
    pi[x] = find(d,p.ldescent(x));

  pi.setClassCount(d.size());

  return;
}

void lStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p)

/*
  This function writes in pi the partition of p according to the (left)
  weak Bruhat links which are equivalences for the W-graph. In other words,
  x is equivalent to sx if sx > x and the left descent set of x is not
  contained in the left descent set of sx; this means that there is a t,
  not commuting with s, in the left descent set of x s.t. x and sx are
  in the same left chain for {s,t}.
*/

{
  static bits::BitMap b(0);
  static Fifo<coxtypes::CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(p.size());
  Ulong count = 0;

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x) {
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[x] = count;
    orbit.push(x);
    while (orbit.size()) {
      coxtypes::CoxNbr z = orbit.pop();
      for (coxtypes::Generator s = 0; s < p.rank(); ++s) {
	coxtypes::CoxNbr sz = p.lshift(z,s);
	if (b.getBit(sz))
	  continue;
	GenSet fz = p.ldescent(z);
	GenSet fsz = p.ldescent(sz);
	GenSet f = fz & fsz;
	if ((f == fz) || (f == fsz)) /* inclusion */
	  continue;
	b.setBit(sz);
	pi[sz] = count;
	orbit.push(sz);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void lStringEquiv(bits::Partition& pi, const bits::SubSet& q, const schubert::SchubertContext& p)

/*
  Does the partition of the subset q into left string classes. It is assumed
  that q is stable under the equivalence relation.
*/

{
  static bits::BitMap b(0);
  static Fifo<coxtypes::CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(q.size());
  Ulong count = 0;

  for (Ulong j = 0; j < q.size(); ++j) {
    const coxtypes::CoxNbr x = q[j];
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[j] = count;
    orbit.push(x);
    while (orbit.size()) {
      coxtypes::CoxNbr z = orbit.pop();
      for (coxtypes::Generator s = 0; s < p.rank(); ++s) {
	coxtypes::CoxNbr sz = p.lshift(z,s);
	if (b.getBit(sz))
	  continue;
	GenSet fz = p.ldescent(z);
	GenSet fsz = p.ldescent(sz);
	GenSet f = fz & fsz;
	if ((f == fz) || (f == fsz)) /* inclusion */
	  continue;
	if (!q.isMember(sz)) { // q is not stable! this shouldn't happen
	  error::ERRNO = error::ERROR_WARNING;
	  return;
	}
	b.setBit(sz);
	orbit.push(sz);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void rDescentPartition(bits::Partition& pi, const schubert::SchubertContext& p)

/*
  This function writes in pi the partition of p according to the right
  descent sets.
*/

{
  static list::List<GenSet> d(0); /* holds the appearing descent sets */

  pi.setSize(p.size());
  d.setSize(0);

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
    insert(d,p.rdescent(x));

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
    pi[x] = find(d,p.rdescent(x));

  pi.setClassCount(d.size());

  return;
}

void rStringEquiv(bits::Partition& pi, const schubert::SchubertContext& p)

/*
  Same as lStringEquiv, but on the other side.
*/

{
  static bits::BitMap b(0);
  static Fifo<coxtypes::CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(p.size());
  Ulong count = 0;

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x) {
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[x] = count;
    orbit.push(x);
    while (orbit.size()) {
      coxtypes::CoxNbr z = orbit.pop();
      for (coxtypes::Generator s = 0; s < p.rank(); ++s) {
	coxtypes::CoxNbr zs = p.rshift(z,s);
	if (b.getBit(zs))
	  continue;
	GenSet fz = p.rdescent(z);
	GenSet fzs = p.rdescent(zs);
	GenSet f = fz & fzs;
	if ((f == fz) || (f == fzs)) /* inclusion */
	  continue;
	b.setBit(zs);
	pi[zs] = count;
	orbit.push(zs);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}

void rStringEquiv(bits::Partition& pi, const bits::SubSet& q, const schubert::SchubertContext& p)

/*
  Same as lStringEquiv, but on the other side.
*/

{
  static bits::BitMap b(0);
  static Fifo<coxtypes::CoxNbr> orbit;

  b.setSize(p.size());
  b.reset();

  pi.setSize(q.size());
  Ulong count = 0;

  for (Ulong j = 0; j < q.size(); ++j) {
    const coxtypes::CoxNbr x = q[j];
    if (b.getBit(x))
      continue;
    b.setBit(x);
    pi[j] = count;
    orbit.push(x);
    while (orbit.size()) {
      coxtypes::CoxNbr z = orbit.pop();
      for (coxtypes::Generator s = 0; s < p.rank(); ++s) {
	coxtypes::CoxNbr zs = p.rshift(z,s);
	if (b.getBit(zs))
	  continue;
	GenSet fz = p.rdescent(z);
	GenSet fzs = p.rdescent(zs);
	GenSet f = fz & fzs;
	if ((f == fz) || (f == fzs)) /* inclusion */
	  continue;
	if (!q.isMember(zs)) { // q is not stable! this shouldn't happen
	  error::ERRNO = error::ERROR_WARNING;
	  return;
	}
	b.setBit(zs);
	orbit.push(zs);
      }
    }
    count++;
  }

  pi.setClassCount(count);

  return;
}


/*
  This is the most delicate of the partition functions. It is the maximal
  refinement of the right descent partition under right star operations.
  In other words, two elements x and y are in the same class for this
  partition, if for each right star-word a (i.e. a sequence of right
  star-operations), x*a and y*a have the same right descent set.

  The algorithm is very much like the minimization algorithm for a finite
  state automaton.

  NOTE : this could probably be simplified with a PartitionIterator; be
  wary though of modifications in pi during the loop.
*/
void rGeneralizedTau(bits::Partition& pi, schubert::SchubertContext& p)
{
  static bits::Permutation v(0);
  static list::List<Ulong> b(0);
  static list::List<Ulong> cc(0); // sizes of parts of the partition |pi|
  static list::List<Ulong> a(0);

  /* initialize pi with partition into right descent sets */

  Ulong prev;
  rDescentPartition(pi,p);
  v.setSize(pi.size());

  do {
    prev = pi.classCount();

    /* refine */

    for (Ulong r = 0; r < p.nStarOps(); ++r) {

      pi.sortI(v);   // set |v| to inverse standardization of partition values
      Ulong count = pi.classCount();
      cc.setSize(count);
      cc.setZero();

      for (Ulong j = 0; j < pi.size(); ++j)
	cc[pi[j]]++; // count each class of |pi|

      Ulong i = 0;

      for (Ulong c = 0; c < pi.classCount(); ++c) {

	coxtypes::CoxNbr x = v[i]; /* first element in class */

	if (p.star(x,r) == coxtypes::undef_coxnbr)
	  goto next_class;

	/* find possibilities for v[.]*r */

	b.setSize(0);

	for (Ulong j = 0; j < cc[c]; ++j) {
	  assert(pi[v[i]]==pi[v[i+j]]); // we traverse a class of |pi|
	  auto star = p.star(v[i+j],r);
	  assert(star != coxtypes::undef_coxnbr); // same descent set, same stars
	  Ulong cr = pi[star];
	  insert(b,cr); // add |cr| to list of class values if new
	}

	if (b.size() > 1) { /* there is a refinement */
	  a.setSize(cc[c]);
	  for (Ulong j = 0; j < a.size(); ++j)
	    a[j] = find(b,pi[p.star(v[i+j],r)]);
	  for (Ulong j = 0; j < cc[c]; ++j) {
	    if (a[j] > 0)
	      pi[v[i+j]] = count+a[j]-1; // make all but one into a new class
	  }
	  count += b.size()-1;
	}

      next_class:
	i += cc[c];
	continue;
      }

      pi.setClassCount(count);

    }

  } while (prev < pi.classCount());

  return;
} // |rGeneralizedTau|

// Like |rGeneralizedTau|, but on the left.
void lGeneralizedTau(bits::Partition& pi, schubert::SchubertContext& p)
{
  static bits::Permutation v(0);
  static list::List<Ulong> b(0);
  static list::List<Ulong> cc(0);
  static list::List<Ulong> a(0);

  /* initialize pi with partition into right descent sets */

  Ulong prev;
  lDescentPartition(pi,p);
  v.setSize(pi.size());

  do {
    prev = pi.classCount();

    /* refine */

    for (Ulong r = p.nStarOps(); r < 2*p.nStarOps(); ++r)
    {

      pi.sortI(v);   /* sort partition */
      Ulong count = pi.classCount();
      cc.setSize(count);
      cc.setZero();

      for (Ulong j = 0; j < pi.size(); ++j)
	cc[pi[j]]++; // count each class of |pi|

      Ulong i = 0;

      for (Ulong c = 0; c < pi.classCount(); ++c) {

	coxtypes::CoxNbr x = v[i]; /* first element in class */

	if (p.star(x,r) == coxtypes::undef_coxnbr)
	  goto next_class;

	/* find possibilities for v[.]*r */

	b.setSize(0);

	for (Ulong j = 0; j < cc[c]; ++j) {
	  assert(pi[v[i]]==pi[v[i+j]]); // we traverse a class of |pi|
	  auto star = p.star(v[i+j],r);
	  assert(star != coxtypes::undef_coxnbr); // same descent set, same stars
	  Ulong cr = pi[star];
	  insert(b,cr); // add |cr| to list of class values if new
	}

	if (b.size() > 1) { /* there is a refinement */
	  a.setSize(cc[c]);
	  for (Ulong j = 0; j < a.size(); ++j)
	    a[j] = find(b,pi[p.star(v[i+j],r)]);
	  for (Ulong j = 0; j < cc[c]; ++j) {
	    if (a[j] > 0)
	      pi[v[i+j]] = count+a[j]-1; // make all but one into a new class
	  }
	  count += b.size()-1;
	}

      next_class:
	i += cc[c];
	continue;
      }

      pi.setClassCount(count);

    }

  } while (prev < pi.classCount());

  return;
}


};

/*****************************************************************************

        Chapter II -- W-graph construction

  This section defines functions for the construction of W-graphs :

   - graph<l/r/b>(kl) : the graph part only, no descent sets;
   - W_graph<l/r/b>(kl) : constructs a W-graph directly from the K-L data;
   - W_graph<l/r/b>(q,kl) : the same, restricted to a subset;

 *****************************************************************************/

namespace cells {

template<char side>
  wgraph::OrientedGraph graph(kl::KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  auto desc =  [&p](coxtypes::CoxNbr x)
    { return side=='r' ? p.rdescent(x) : side=='l' ? p.ldescent(x) : p.descent(x); };

  wgraph::OrientedGraph X(kl.size());

  for (coxtypes::CoxNbr y = 0; y < kl.size(); ++y) {
    const kl::MuRow& mu = kl.muList(y);
    for (Ulong j = 0; j < mu.size(); ++j) {
      if (mu[j].mu != 0) {
	coxtypes::CoxNbr x = mu[j].x;
	if (desc(x) != desc(y)) /* make an edge from x to y */
	  X.edge(x).append(y);
      }
    }
  }

  for (coxtypes::CoxNbr y = 0; y < kl.size(); ++y) {
    const schubert::CoxNbrList& c = p.hasse(y);
    for (Ulong j = 0; j < c.size(); ++j) {
      if ((desc(c[j])&desc(y)) != desc(c[j]))
	X.edge(c[j]).append(y);
      if ((desc(c[j])&desc(y)) != desc(y))
	X.edge(y).append(c[j]);
    }
  }

  return X;
}

template wgraph::OrientedGraph graph<'l'>(kl::KLContext& kl);
template wgraph::OrientedGraph graph<'r'>(kl::KLContext& kl);
template wgraph::OrientedGraph graph<'b'>(kl::KLContext& kl);

/*
  Construct a W-graph directly from the K-L data. In other
  words, we construct a graph with vertex set the elements of p; for each
  x < y s.t. mu(x,y) != 0, and L(x) != L(y), we set an edge from x to y
  if L(y) \subset L(x), from y to x if L(x) \subset L(y); the coefficient
  of this edge will be mu(x,y) in both cases.

  Also, to each vertex is associated the descent set L(x).

  Assumes that the mu-table has been filled.

  NOTE : this should be changed when there will no longer be a mu-table
  in the current sense.
*/
template<char side> // one of 'l', 'r', 'b'
    wgraph::WGraph W_graph(kl::KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  wgraph::WGraph X(kl.size());
  X.graph() = cells::graph<'l'>(kl);

  // fill in coefficients

  for (coxtypes::CoxNbr y = 0; y < kl.size(); ++y)
  {
    coxtypes::Length ly = p.length(y);
    wgraph::CoeffList& c = X.coeffList(y);
    const wgraph::EdgeList& e = X.edge(y);
    c.setSize(e.size());
    for (Ulong i = 0; i < c.size(); ++i) // traverse both |e| and |c|
    {
      coxtypes::CoxNbr x = e[i];
      coxtypes::Length lx = p.length(x);
      if (lx < ly or lx == ly+1)
	c[i] = 1;
      else
	c[i] = kl.mu(y,x);
    }
  }

  // fill in descent sets

  for (coxtypes::CoxNbr y = 0; y < kl.size(); ++y)
    X.descent(y) =
      side=='l' ? p.ldescent(y) : side=='r'? p.rdescent(y) : p.descent(y);

  return X;
}

template wgraph::WGraph W_graph<'l'>(kl::KLContext& kl);
template wgraph::WGraph W_graph<'r'>(kl::KLContext& kl);
template wgraph::WGraph W_graph<'b'>(kl::KLContext& kl);


/*
  Construct the right/left/two-sided W-graph for the subset q. It is
  assumed that q is a union of left cells (typically, q might be a right
  descent class, or one of the classes provided by GeneralizedTau).

  The difference with the full lWGraph, is that we do _not_ assume that
  the mu-coefficients have already been computed; we compute them as
  needed.

  It is assumed that q is sorted in increasing order.
*/
template<char side>
  wgraph::WGraph W_graph(const bits::SubSet& q, kl::KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();
  auto desc =  [&p](coxtypes::CoxNbr x)
    { return side=='r' ? p.rdescent(x) : side=='l' ? p.ldescent(x) : p.descent(x); };

  wgraph::WGraph X(q.size());

  wgraph::OrientedGraph& Y = X.graph();
  Y.reset();

  for (Ulong j = 0; j < q.size(); ++j)
  {

    coxtypes::CoxNbr y = q[j];
    coxtypes::Length ly = p.length(y);

    // set descent set

    X.descent(j) = desc(y);

    bitmap::BitMap b = p.closure(y);

    for (Ulong pos = 0; pos < q.size(); ++pos)
      if (b.is_member(q[pos]))
      {
	coxtypes::CoxNbr x = q[pos];
	coxtypes::Length dl = ly-p.length(x);

	if (dl%2 == 0)
	  continue;
	if (dl == 1)
	{ /* found a hasse edge */
	  if ((desc(x)&desc(y)) != desc(x))
	  {
	    Y.edge(pos).append(j);
	    X.coeffList(pos).append(1);
	  }
	  if ((desc(x)&desc(y)) != desc(y))
	  {
	    Y.edge(j).append(pos);
	    X.coeffList(j).append(1);
	  }
	  continue;
	}

	klsupport::KLCoeff mu = kl.mu(x,y);

	if (mu != 0 and desc(x) != desc(y))
	{
	  Y.edge(pos).append(j);
	  X.coeffList(pos).append(mu);
	}
      } // |for(pos)|
  } // |for(j)|
  return X;
}

template wgraph::WGraph W_graph<'l'>(const bits::SubSet& q, kl::KLContext& kl);
template wgraph::WGraph W_graph<'r'>(const bits::SubSet& q, kl::KLContext& kl);
template wgraph::WGraph W_graph<'b'>(const bits::SubSet& q, kl::KLContext& kl);




};

/*****************************************************************************

        Chapter III -- Graph construction for unequal parameters.

  In the case of unequal parameters, there is no W-graph associated to the
  situation. However, there is still an oriented graph, defined exactly
  as in the case of equal parameters, in terms of the mu-coefficients. Here
  in fact mu^s_{x,y} is defined only when ys > y, xs < x; hence the condition
  that R(x) and R(y) are distinct is automatically fulfilled. So the edges
  from y point to the ws s.t. wy > w, and to the z in muTable(s,y) s.t.
  mu^s_{z,y} != 0.

 *****************************************************************************/

namespace cells {

/*
  Return the left/right/two-sided graph corresponding to the edges in the context.
  It is assumed that the mu-table has been filled.
*/
template<char side>
  wgraph::OrientedGraph graph(uneqkl::KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();
  wgraph::OrientedGraph X(kl.size());
  const Lflags S = constants::lt_mask[kl.rank()];

  if (side=='b')
    X = graph<'r'>(kl); // start with right edges
  else // reset X
    for (coxtypes::CoxNbr y = 0; y < X.size(); ++y) {
      wgraph::EdgeList& e = X.edge(y);
      e.setSize(0);
    }

  // make edges (on the left when two-sided)
  for (coxtypes::CoxNbr y = 0; y < X.size(); ++y) {
    coxtypes::CoxNbr yi = side=='r' ? y : kl.inverse(y);
    for (GenSet f = ~p.rdescent(y) & S; f; f &= (f-1)) {
      coxtypes::Generator s = constants::firstBit(f);
      const uneqkl::MuRow& muRow = kl.muList(s,y);
      for (Ulong j = 0; j < muRow.size(); ++j) {
	wgraph::EdgeList& e =
	  X.edge(side=='r' ? muRow[j].x : kl.inverse(muRow[j].x));
	if (side=='b')
	  insert(e,wgraph::Edge(yi));
	else
	  e.append(yi);
      }
      wgraph::Vertex sy  = kl.inverse(p.shift(y,s));
      wgraph::EdgeList& e = X.edge(sy);
      if (side=='b')
	insert(e,wgraph::Edge(yi));
      else
	e.append(yi);
    }
   }

  if (side=='l') // then we must sort edgelists
    for (coxtypes::CoxNbr y = 0; y < X.size(); ++y) {
      wgraph::EdgeList& e = X.edge(y);
      e.sort();
    }

  return X;
}

template wgraph::OrientedGraph graph<'l'>(uneqkl::KLContext& kl);
template wgraph::OrientedGraph graph<'r'>(uneqkl::KLContext& kl);
template wgraph::OrientedGraph graph<'b'>(uneqkl::KLContext& kl);


};
/*****************************************************************************

        Chapter IV -- Utilities

  This section defines some utility functions :

   - checkClasses(pi,p) : checks the classes of a refined partition;

 *****************************************************************************/

namespace cells {

coxtypes::CoxNbr checkClasses
  (const bits::Partition& pi, const schubert::SchubertContext& p)

/*
  This function checks if the classes of a refined partition are stable
  under weak equivalence, as they should.
*/

{
  static bits::Permutation v(0);
  static bits::Partition pi_q(0);
  static bits::SubSet q(0);

  q.setBitMapSize(p.size());

  v.setSize(pi.size());
  pi.sortI(v);

  Ulong i = 0;

  for (Ulong j = 0; j < pi.classCount(); ++j) {
    q.reset();
    for (; pi(v[i]) == j; ++i) {
      q.add(v[i]);
    }
    lStringEquiv(pi_q,q,p);
    if (error::ERRNO) {
      printf("error in class #%lu\n",j);
      return q[0];
    }
  }

  return 0;
}

};

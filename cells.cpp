/*
  This is cells.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "cells.h"

#include <algorithm>
#include "wgraph.h"

#include "stack.h"

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

template<char side> // one of 'l', 'r'
  bits::Partition descent_partition(const schubert::SchubertContext& p)
{
  auto desc =  [&p](coxtypes::CoxNbr x) -> GenSet
  { return p.descent_set<side>(x); };

  return bits::Partition(p.size(),desc);
} // |descent_partition|

template bits::Partition descent_partition<'l'>
  (const schubert::SchubertContext& p);
template bits::Partition descent_partition<'r'>
  (const schubert::SchubertContext& p);

/*
  This is the most delicate of the partition functions. It is the maximal
  refinement of the right or left descent partition under right respectivelu
  left star operations. In other words, two elements $x$ and $y$ are in the same
  class for this partition, if for each composition $f$ of right
  star-operations, the lements $f(x)$ and $f(y)$ have the same right descent set.

  The algorithm is very much like the construction of a finite state automaton
  for a regular language, where a coarse partition of the set of words is
  successively refined as applying the extension by certain letters detects that
  equivalence classes need to be split up, and this is repeated up to the point
  where no further refinement is obtained. A similar algorithm also occurs when
  in the presence of a collection of possibly recursive type definitions one
  wants to decide (structural) equivalence of type expressions: by definition
  they are equivalent if no sequence of type operations can detect a difference.
*/
template<char side> // one of 'l', 'r'
  bits::Partition generalized_tau(schubert::SchubertContext& p)
{
  /* initialize pi with partition into right descent sets */

  bits::Partition pi = descent_partition<side>(p);

  bool go_on; // stupid C++ rule requires a declaration outside |do..while();|
  do
  {
    bits::Permutation v = pi.inverse_standardization();
    containers::vector<unsigned> class_size(pi.classCount(),0);
    for (Ulong j = 0; j < pi.size(); ++j)
      ++class_size[pi[j]]; // count each class of |pi|

    Ulong i = 0;
    containers::vector<bits::Partition> refinements;
    refinements.reserve(pi.classCount());

    for (Ulong c=0; c<pi.classCount(); i += class_size[c], ++c) // do class |c|
    {
      const auto* star = p.star_base<side>(v[i]); // for first element in class

      bitmap::BitMap star_set(p.nStarOps());
      for (Ulong r = 0; r < p.nStarOps(); ++r)
	star_set.set_to(r,p.in_context(star[r]));

      auto star_view = // what the values of |pi| look like one- *-op away
	[&,i] (unsigned j) -> containers::vector<Ulong>
        { containers::vector<Ulong> result;
	  for (Ulong r : star_set)
	    result.push_back(pi(p.star_base<side>(v[i+j])[r]));
	  return result;
	};
      refinements.emplace_back(class_size[c],star_view);
    } // |for(c)| loop over classes for |pi|
    go_on = pi.refine(refinements);
  } while(go_on); // whether to do another iteration

  return pi;
} // |generalized_tau|

template bits::Partition generalized_tau<'l'> (schubert::SchubertContext& p);
template bits::Partition generalized_tau<'r'> (schubert::SchubertContext& p);

/*
  Return the partition of |p| into left cells --- in the case of an incomplete
  context, they will be the cells defined by the links in the graph.

  The idea will be to minimize K-L computations. We proceed as follows. First,
  we determine the generalized-tau partition of the context. Then, we look
  at the star-orbits among the tau-classes, and decompose one representative
  of each into cells; we propagate the cells using star-operations.
*/
template<char side> // one of 'l', 'r',
  bits::Partition cells(kl::KLContext& kl)
{
  schubert::SchubertContext& p = kl.schubert();

  bits::SubSet seen(p.size());
  containers::vector<Ulong> cell_sizes; // sizes of cells on |seen| list

  constexpr char opposite_side = 'l' + 'r' - side;
  bits::Partition pi = generalized_tau<opposite_side>(p);

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
  {
    if (seen.is_member(x))
      continue;

    bits::SubSet gen_tau_class = pi.class_of(x);

    // compute Wgraph for the class, then find its strong components (cells)
    bits::Partition qcells =
      W_graph<side>(gen_tau_class,kl).graph().cells(); // no induced graph needed

    /* the fifo-list orbit is used to traverse the *-orbit of the first
       element of the current generalized-tau class */
    containers::queue<Ulong> orbit { seen.size() }; // start with next-to-add
    containers::vector<Ulong> comp_sizes;
    comp_sizes.reserve(qcells.classCount());

    // get cell sizes and record their members; they fill |x|'s |gen_tau_class|
    for (bits::PartitionIterator pit(qcells); pit; ++pit)
    {
      const bits::Set& comp = pit(); // a strongly connected component (cell)
      comp_sizes.push_back(comp.size());
      // mark all elements of these cells as seen
      for (coxtypes::CoxNbr y : comp)
	seen.add(gen_tau_class[y]);
    }

    // record cell structure by copying sizes-list |comp_sizes| to |cell_sizes|
    cell_sizes.insert(cell_sizes.end(), comp_sizes.begin(),comp_sizes.end());

    /* from the generalized tau class of |x|, partitioned into cells, we can
       cheaply deduce other, isomorphic, packets of cells by repeatedly applying
       opposite-side star operations
     */
    while (not orbit.empty())
    {
      Ulong c = orbit.front(); // index of start of cell packet on |seen|
      // the packet runs from |seen[c]| up to |seen[c+gen_tau_class.size()]|
      // and is partitioned into consecutive cells of sizes |comp_sizes|
      orbit.pop();

      coxtypes::CoxNbr z = seen[c]; // the head of the packet
      const auto* z_star = p.star_base<opposite_side>(z);

      for (coxtypes::StarOp j = 0; j < p.nStarOps(); ++j)
      {
	coxtypes::CoxNbr zj = z_star[j];

	if (zj == coxtypes::undef_coxnbr)
	  continue;
	if (seen.is_member(zj))
	  continue;

	GenSet d = p.descent_set<opposite_side>(zj);

	// mark start of new cell packet that will be added to |seen| below
	orbit.push(seen.size()); // index of next-to-add packet leader

	// add image of |gen_tau_class| obtained by applying star operation |j|
	for (Ulong i = 0; i < gen_tau_class.size(); ++i)
	{
	  coxtypes::CoxNbr y = seen[c+i];
	  coxtypes::CoxNbr yj = p.star_base<opposite_side>(y)[j];
	  assert(p.descent_set<opposite_side>(yj)==d);
	  seen.add(yj);
	}

	// record new cells by copying sizes-list |comp_sizes| to |cell_sizes|
	cell_sizes.insert(cell_sizes.end(),
			  comp_sizes.begin(),comp_sizes.end());
      } // |for(j)|, different star operations

    } // |while(not orbit.empty())|
  } // |for(x)||

  return bits::Partition(&seen[0],seen.size(),cell_sizes);
} // |cells|

template bits::Partition cells<'l'> (kl::KLContext& kl);
template bits::Partition cells<'r'> (kl::KLContext& kl);


/*
  This function computes the two-sided cells in the context. There are
  certainly better ways to do this, but I'm afraid I don't know enough
  to do it other than by filling in all the mu's ... [Fokko]
*/
template<> bits::Partition cells<'b'>(kl::KLContext& kl)
{
  kl.fillMu();
  return cells::W_graph<'b'>(kl).graph().cells(); // strong components partition
}

/*
  Return the partition of |p| according to the left/right weak Bruhat links
  which are equivalences for the W-graph. In other words, $x$ is equivalent to
  $sx$ if $sx > x$ and the left descent set of $x$ is not contained in the left
  descent set of $sx$; this means that there is a $t$, not commuting with $s$,
  in the left descent set of $x$ such that $x$ and $sx$ are in the same left
  chain for $\{s,t\}$.
*/
template<char side> // one of 'l', 'r'
  bits::Partition string_equiv(const schubert::SchubertContext& p)
{
  containers::vector<Ulong> classify(p.size(),Ulong(-1));
  Ulong count = 0; // value that will caracterize next orbit

  for (coxtypes::CoxNbr start = 0; start<p.size(); ++start)
    if (classify[start]>=count)
    { // now |start| will start a new orbit
      containers::queue<coxtypes::CoxNbr> orbit { start };
      classify[start] = count;
      do
      {
	coxtypes::CoxNbr x = orbit.front();
	orbit.pop();
	assert(classify[x]==count); // older orbits should not be encountered

	for (coxtypes::Generator s = 0; s < p.rank(); ++s) // visit neighbours
	{
	  coxtypes::CoxNbr y = side=='l' ? p.lshift(x,s) : p.rshift(x,s);
	  if (classify[y]<=count) // either an old element or already in orbit
	    continue; // don't even try to see if this is in same orbit
	  const GenSet dx = p.descent_set<side>(x), dy = p.descent_set<side>(y);
	  if ((dx & dy)!=dx and (dx & dy)!=dy) // inclusion neither way
	  {
	    classify[y] = count; // we flag |y| as part of the current orbit
	    orbit.push(y); // and record it for inspection of its neighbours
	  }
	} // |for(s)||
      } while (not orbit.empty());
      ++count;
    } // |for(x) if (classify[start]>=count)|
  assert(std::all_of(classify.begin(),classify.end(),
		     [count,&classify]
		     (coxtypes::CoxNbr x){ return classify[x]<count;}));

  return bits::Partition(std::move(classify),count);
} // |string_equiv|

template bits::Partition string_equiv<'l'>(const schubert::SchubertContext& p);
template bits::Partition string_equiv<'r'>(const schubert::SchubertContext& p);

/*
  Do the partition of the subset q into left string classes. It is assumed
  that q is stable under the equivalence relation. (This function is unused.)
*/
template<char side> // one of 'l', 'r'
  bits::Partition string_equiv
    (const bits::SubSet& q, const schubert::SchubertContext& p)
{
  containers::vector<Ulong> classify(q.size(),Ulong(-1));
  auto q_index = [&q](coxtypes::CoxNbr x) -> Ulong
  { return std::lower_bound(q.begin(),q.end(),x)-q.begin(); };

  Ulong count = 0;

  for (Ulong i=0; i<q.size(); ++i)
    if (classify[i]>=count)
    { // now |q[i]| will start a new orbit
      classify[i] = count; // this is the first orbit element
      containers::queue<coxtypes::CoxNbr> orbit
	{ static_cast<coxtypes::CoxNbr>(q[i]) };
      do
      {
	coxtypes::CoxNbr x = orbit.front();
	orbit.pop();
	assert(classify[q_index(x)]==count); // older orbits should not occur

	for (coxtypes::Generator s = 0; s < p.rank(); ++s) // visit neighbours
	{
	  coxtypes::CoxNbr y = side=='l' ? p.lshift(x,s) : p.rshift(x,s);
	  const GenSet dx = p.descent_set<side>(x), dy = p.descent_set<side>(y);
	  const GenSet common = dx & dy;
	  if (common != dx and common != dy) // inclusion neither way
          {
	    if (not q.is_member(y))
	    { // we found that |q| is not stable! this shouldn't happen
	      error::ERRNO = error::ERROR_WARNING;
	      return bits::Partition(0);
	    }
	    classify[q_index(y)] = count; // mark |y| as in the current orbit
	    orbit.push(y); // and record it for inspection of its neighbours
	  }
	} // |for(s)||
      } while (not orbit.empty());
      ++count; // pass to a next orbit
    } // |for(x)| |if(not seen)|

  assert(std::all_of(classify.begin(),classify.end(),
		     [count,&classify]
		     (coxtypes::CoxNbr x){ return classify[x]<count;}));
  return bits::Partition(std::move(classify),count);
} // |string_equiv|

template bits::Partition string_equiv<'l'>
  (const bits::SubSet& q, const schubert::SchubertContext& p);
template bits::Partition string_equiv<'r'>
  (const bits::SubSet& q, const schubert::SchubertContext& p);


/*****************************************************************************

        Chapter II -- W-graph construction

  This section defines functions for the construction of W-graphs :

   - graph<l/r/b>(kl) : the graph part only, no descent sets;
   - W_graph<l/r/b>(kl) : constructs a W-graph directly from the K-L data;
   - W_graph<l/r/b>(q,kl) : the same, restricted to a subset;

 *****************************************************************************/

template<char side> // one of 'l', 'r', 'b'
  wgraph::OrientedGraph graph(kl::KLContext& kl)
{
  const schubert::SchubertContext& p = kl.schubert();

  wgraph::OrientedGraph X(kl.size());

  for (coxtypes::CoxNbr y = 0; y < kl.size(); ++y) {
    const kl::MuRow& mu = kl.muList(y);
    for (Ulong j = 0; j < mu.size(); ++j) {
      if (mu[j].mu != 0) {
	coxtypes::CoxNbr x = mu[j].x;
	if (p.descent_set<side>(x) != p.descent_set<side>(y))
	  X.edge(x).append(y); // add an edge from |x| to |y|
      }
    }
  }

  for (coxtypes::CoxNbr y = 0; y < kl.size(); ++y) {
    const schubert::CoxNbrList& c = p.hasse(y);
    for (Ulong j = 0; j < c.size(); ++j)
    { const Lflags dx = p.descent_set<side>(c[j]), dy=p.descent_set<side>(y);
      if ((dx&dy) != dx)
	X.edge(c[j]).append(y);
      if ((dx&dy) != dy)
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
  X.graph() = cells::graph<side>(kl);

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
  X.setSize(q.size()); // this is still necessary

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

/*
  Return the left/right/two-sided graph corresponding to the edges in the
  context. It is assumed that the mu-table has been filled.
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

/*****************************************************************************

        Chapter IV -- Utilities

  This section defines some utility functions :

   - checkClasses(pi,p) : checks the classes of a refined partition;

 *****************************************************************************/


/*
  Check if the classes of a refined partition are stable
  under weak equivalence, as they should.
*/
coxtypes::CoxNbr checkClasses
  (const bits::Partition& pi, const schubert::SchubertContext& p)
{
  static bits::Partition pi_q(0);
  static bits::SubSet q(0);

  q.setBitMapSize(p.size());

  bits::Permutation v = pi.inverse_standardization();

  Ulong i = 0;

  for (Ulong j = 0; j < pi.classCount(); ++j) {
    q.reset();
    for (; pi(v[i]) == j; ++i) {
      q.add(v[i]);
    }
    pi_q=string_equiv<'l'>(q,p);
    if (error::ERRNO) {
      printf("error in class #%lu\n",j);
      return q[0];
    }
  }

  return 0;
} // |checkClasses|

};

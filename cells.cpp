/*
  This is cells.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "cells.h"
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
    { return side=='r' ? p.rdescent(x) : p.ldescent(x); };

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
  in the presence of a collcetion of possibly recursive type definitions one
  wants to decide (structural) equivalence of type expressions: by definition
  they are equivalent if no sequence of type operations can detect a difference.

  NOTE : this could probably be simplified with a PartitionIterator; be
  wary though of modifications in pi during the loop.
*/
template<char side> // one of 'l', 'r'
  bits::Partition generalized_tau(schubert::SchubertContext& p)
{
  /* initialize pi with partition into right descent sets */

  bits::Partition pi = descent_partition<side>(p);

  Ulong prev_count; // hold previous count to measure progress; must be outside loop
  do
  {
    prev_count = pi.classCount(); // base value for |do...while(...)| loop body

    for (Ulong r = 0; r < p.nStarOps(); ++r) // try to refine |pi| by *-op |r|
    {
      bits::Permutation v = pi.inverse_standardization();
      containers::vector<unsigned> class_size(pi.classCount(),0);

      for (Ulong j = 0; j < pi.size(); ++j)
	++class_size[pi[j]]; // count each class of |pi|

      Ulong i = 0;
      const auto n_classes = pi.classCount(); // fix this since loop increases it
      for (Ulong c=0; c<n_classes; i += class_size[c], ++c) // handle class |c|
      {
	const coxtypes::CoxNbr x = v[i]; // first element in class
	const coxtypes::CoxNbr x_star = p.star_base<side>(x)[r];
	if (x_star == coxtypes::undef_coxnbr)
	  continue; // if this one fails, the whold class will (?)

	containers::multimap<Ulong,Ulong> star_values
	  { std::make_pair(pi[x_star],x) };

	for (Ulong j = 1; j < class_size[c]; ++j)
	{
	  const coxtypes::CoxNbr xx = v[i+j];
	  assert(pi[x]==pi[xx]); // we traverse the |pi|-class of |x|
	  const coxtypes::CoxNbr xx_star = p.star_base<side>(xx)[r];
	  assert(p.in_context(xx_star)); // same descent set, same stars
	  star_values.insert(std::make_pair(pi[xx_star],xx));
	}

	for (auto it = star_values.upper_bound(star_values.begin()->first);
	     it != star_values.end(); // increasing |it| happens in loop body
	     pi.incr_class_count() ) // traverse any new subclasses to split off
	  for (auto key = it->first;
	       it != star_values.end() and it->first==key;
	       ++it)
	    pi[it->second] = pi.classCount(); // give new classifier value to |xx|

      } // |for(c)| (class for |pi| before loop)

    } // |for(r)| (star operation)
  } while(pi.classCount() > prev_count); // whether to do another iteration

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

  auto op_desc = [&p](unsigned x) -> GenSet
    { return side=='l' ? p.rdescent(x) : p.ldescent(x) ; };

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
  {
    if (seen.is_member(x))
      continue;

    bits::SubSet gen_tau_class = pi.class_of(x);

    // compute Wgraph for the class, then find its strong components (cells)
    bits::Partition qcells(0); // will represent the strong components
    wgraph::WGraph X = W_graph<side>(gen_tau_class,kl);
    X.graph().cells(qcells,nullptr); // get classes |qcells|, no induced graph

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

	GenSet d = op_desc(zj);

	// mark start of new cell packet that will be added to |seen| below
	orbit.push(seen.size()); // index of next-to-add packet leader

	// add image of |gen_tau_class| obtained by applying star operation |j|
	for (Ulong i = 0; i < gen_tau_class.size(); ++i)
	{
	  coxtypes::CoxNbr y = seen[c+i];
	  coxtypes::CoxNbr yj = p.star_base<opposite_side>(y)[j];
	  assert(op_desc(yj)==d);
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

  bits::Partition pi;
  wgraph::WGraph X = cells::W_graph<'b'>(kl);
  X.graph().cells(pi); // partition graph into strong components
  return pi;
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
  auto desc =  [&p](coxtypes::CoxNbr x) -> GenSet
    { return side=='r' ? p.rdescent(x) : p.ldescent(x); };

  bitmap::BitMap seen(p.size());

  bits::Partition pi(p.size());
  Ulong count = 0;

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x)
    if (not seen.is_member(x))
    {
      seen.insert(x);
      pi[x] = count;
      containers::queue<coxtypes::CoxNbr> orbit { x };
      do
      {
	coxtypes::CoxNbr z = orbit.front();
	orbit.pop();

	for (coxtypes::Generator s = 0; s < p.rank(); ++s)
	{
	  coxtypes::CoxNbr sz = side=='l' ? p.lshift(z,s) : p.rshift(z,s);
	  if (seen.is_member(sz))
	    continue;
	  GenSet fz = desc(z);
	  GenSet fsz = desc(sz);
	  GenSet f = fz & fsz;
	  if (f != fz and f != fsz) // inclusion neither way
          {
	    seen.insert(sz);
	    pi[sz] = pi.classCount();
	    orbit.push(sz);
	  }
	} // |for(s)||
      } while (not orbit.empty());
      pi.incr_class_count();
    } // |for(x) if (not seen.is_member(x))|

  return pi;
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
  auto desc =  [&p](coxtypes::CoxNbr x) -> GenSet
    { return side=='r' ? p.rdescent(x) : p.ldescent(x); };

  bitmap::BitMap seen(p.size());

  bits::Partition pi(q.size());
  Ulong count = 0;

  for (Ulong j = 0; j < q.size(); ++j)
  {
    const coxtypes::CoxNbr x = q[j];
    if (not seen.is_member(x))
    {
      seen.insert(x);
      pi[j] = count;
      containers::queue<coxtypes::CoxNbr> orbit { x };
      do
      {
	coxtypes::CoxNbr z = orbit.front();
	orbit.pop();
	for (coxtypes::Generator s = 0; s < p.rank(); ++s)
	{
	  coxtypes::CoxNbr sz = side=='l' ? p.lshift(z,s) : p.rshift(z,s);
	  if (seen.is_member(sz))
	    continue;
	  GenSet fz = desc(z);
	  GenSet fsz = desc(sz);
	  GenSet f = fz & fsz;
	  if (f != fz and f != fsz) // inclusion neither way
          {
	    if (not q.is_member(sz))
	    { // q is not stable! this shouldn't happen
	      error::ERRNO = error::ERROR_WARNING;
	      return pi;
	    }
	    seen.insert(sz);
	    pi[sz] = pi.classCount();
	    orbit.push(sz);
	  }
	} // |for(s)||
      } while (not orbit.empty());
      pi.incr_class_count();
      count++;
    } // |if(not seen)|
  } // |for(x)|

  return pi;
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

  auto desc =  [&p](coxtypes::CoxNbr x)
    { return side=='r' ? p.rdescent(x) : side=='l' ? p.ldescent(x)
      : p.descent(x); };

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

/*
  This is posets.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "posets.h"
#include "sl_list.h"

/****************************************************************************

  This file contains some code for the analysis of posets --- this would
  be a subject for a program per se! Currently the only "general" poset
  that appears in the program is the primitive ideal spectrum, i.e. the
  set of cells in a W-graph, considered as an ordered set. We just include
  some functions that are required to normalize and output this poset.

 ****************************************************************************/

namespace posets {

namespace {

  PosetElt first_sink_in
    (const wgraph::OrientedGraph& G, const bitmap::BitMap& b);

}

/****************************************************************************

        Chapter I -- The Poset class

  The Poset class is just a straightforward implementation of a poset; the
  idea is that we will never have to deal with poseets that have more than
  a few thousand elements; therefore we can afford the luxury to carry the
  full incidence graph, as a list of BitMaps. This allows for very fast
  comparison testing, and efficient functions to analyze the poset.

  The following functions are provided :

    - constructors and destructors :

      - Poset();
      - Poset(const Ulong&);
      - Poset(const wgraph::OrientedGraph&);
      - ~Poset();

    - manipulators :

    - accessors :

      - maxima_within(D) : the maximal elements of D in a;
      - hasseDiagram(H) : writes the Hasse diagram in H;
      - isTriangular() : checks if the poset is triangular;
      - size(); (inlined)

 ****************************************************************************/


Poset::Poset()

{}


/*
  Constructs a Poset structure capable of accomodating a Poset of size n.
*/
Poset::Poset(const Ulong& n)
  : d_closure(n,bitmap::BitMap(n)) // start with |n| empty bitmaps
{}


/*
  Construct the poset defined by the graph G, assumed to be acyclic; i.e.,
  the underlying set is the vertex set of G, and x <= y iff there is an
  oriented path in G from y to x (we assume that G describes dominance
  relations, as a matter of convention.)
*/
Poset::Poset(const wgraph::OrientedGraph& G)
  : d_closure(G.size(),bitmap::BitMap(G.size()))
{
  bitmap::BitMap candidates(size());
  candidates.fill();

  while (not candidates.empty())
  {
    PosetElt x = first_sink_in(G,candidates);
    d_closure[x].insert(x);
    for (auto y : G.edge(x))
    {
      assert(not candidates.is_member(y));
      d_closure[x] |= d_closure[y];
    }
    candidates.remove(x); // |candidates.andnot(d_closure[x])| wouldn't do more
  }
}

/******** manipulators ******************************************************/

/******** accessors *********************************************************/


/*
  Return the maximal elements of |D|. Assumes that the poset is in triangular form.

  The algorithm is as follows. The largest element |elt| in |D| is certainly
  maximal. Then remove cl(z) from D, and iterate until reaching the empty set.
*/
containers::vector<Ulong> Poset::maxima_within(bitmap::BitMap D) const
{
  containers::sl_list<Ulong> result;
  Ulong elt = D.capacity();
  while (D.back_up(elt))
  {
    result.push_front(elt);
    D.andnot(d_closure[elt]);
  }

  return result.to_vector();
}


/*
  Check whether the poset is enumerated in a way compatible with the ordering.
  The condition is that x <= y in the poset implies x <= y as numbers. If not,
  it is always possible to permute the poset in order to get such an
  enumeration.
*/
bool Poset::isTriangular() const
{
  for (PosetElt x = 0; x < size()-1; ++x)
  { Ulong n = size();
    d_closure[x].back_up(n); // for side effect; it always returns |true|
    if (n>x)
      return false;
  }
  return true;
}

void Poset::hasseDiagram(wgraph::OrientedGraph& H) const
{
  H.setSize(size());

  for (PosetElt x = 0; x < size(); ++x)
    H.edge(x) = maxima_within(bitmap::BitMap(d_closure[x]).remove(x));
}

/******** input/output ******************************************************/


/*****************************************************************************

        Chapter II -- Auxiliary functions.

  This chapter defines some auxiliary functions used in this module :

    - first_sink_in(G,b) : return the first minimal |x| in subset |b| of |G|

 *****************************************************************************/

namespace {

PosetElt first_sink_in(const wgraph::OrientedGraph& G, const bitmap::BitMap& b)
{
  auto outside_b = [&b](Ulong y) { return not b.is_member(y); };
  for (PosetElt x : b)
    if (std::all_of(G.edge(x).begin(),G.edge(x).end(),outside_b))
      return x;

  return G.size(); // not found indication
}


}; // |namespace|
}; // |namespace posets|

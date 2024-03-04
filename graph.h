/*
  This is graph.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/* type definitions */

#ifndef GRAPH_H  /* guarantee single inclusion */
#define GRAPH_H

#include "globals.h"
#include "list.h"
#include "bits.h"
#include "coxtypes.h"
#include "memory.h"
#include "type.h"

namespace graph {

/* type declarations */

  class CoxGraph;
  using CoxEntry = unsigned short;
  using CoxMatrix = containers::vector<CoxEntry>;


/* constants */

  const Ulong SBITMAP_MAX =
    coxtypes::RANK_MAX/CHAR_BIT
    + (bool)(coxtypes::RANK_MAX % CHAR_BIT);
  /* a CoxNbr should hold at least 2 COXENTRY_MAX elements */
  static const CoxEntry COXENTRY_MAX = 32763;
  static const CoxEntry undef_coxentry = USHRT_MAX;
  static const CoxEntry infty = 0;


/******** function declarations **********************************************/

  containers::vector<GenSet> conjugacy_classes(const CoxGraph& G);
  bool isAffine(CoxGraph& G, GenSet I);
  bool isConnected(CoxGraph& G, GenSet I);
  bool isCrystallographic(CoxGraph& G, GenSet I);
  bool isFinite(CoxGraph& G, GenSet I);
  bool isLoop(CoxGraph& G, GenSet I);
  bool isSimplyLaced(CoxGraph& G, GenSet I);
  bool isTree(CoxGraph& G, GenSet I);
  coxtypes::CoxSize order(CoxGraph& G, GenSet I);
  coxtypes::ParSize quotOrder(CoxGraph& G, GenSet I, GenSet J);
  coxtypes::Generator *standardEnumeration(CoxGraph& G, GenSet I);
  const type::Type& type(CoxGraph& G, GenSet I);
};


/* type definitions */

class graph::CoxGraph
{
 private:
  type::Type d_type;
  coxtypes::Rank d_rank;
  CoxMatrix d_matrix; // flattened, by row
  GenSet d_S; // all generators
  containers::vector<GenSet> d_star; // neighbourhoods of the vertices
  containers::vector<GenSet> d_finite_edges; // list of non-infinite edges
 public:
/* constructors and destructors */
  CoxGraph(const type::Type& x, const coxtypes::Rank& l);
  ~CoxGraph();
/* accessors */
  GenSet component(GenSet I,coxtypes::Generator s) const;
  GenSet extremities(GenSet I) const;
  CoxEntry M(coxtypes::Generator s, coxtypes::Generator t) const
    { return(d_matrix[s*d_rank + t]); } // this defines the matrix layout
  GenSet nodes(GenSet I) const;
  coxtypes::Rank rank() const { return d_rank; }
  GenSet supp() const { return d_S; }
  GenSet star(coxtypes::Generator s) const { return(d_star[s]); }
  GenSet star(GenSet I, coxtypes::Generator s) const
    { return(d_star[s]&I); }
  const containers::vector<GenSet>& finite_edges() const
    { return d_finite_edges; }
  const type::Type& type() const { return d_type; }

}; // |class graph::CoxGraph|

#endif

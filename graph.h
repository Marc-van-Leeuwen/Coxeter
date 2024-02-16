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

  containers::vector<bits::Lflags> conjugacy_classes(const CoxGraph& G);
  bool isAffine(CoxGraph& G, bits::Lflags I);
  bool isConnected(CoxGraph& G, bits::Lflags I);
  bool isCrystallographic(CoxGraph& G, bits::Lflags I);
  bool isFinite(CoxGraph& G, bits::Lflags I);
  bool isLoop(CoxGraph& G, bits::Lflags I);
  bool isSimplyLaced(CoxGraph& G, bits::Lflags I);
  bool isTree(CoxGraph& G, bits::Lflags I);
  coxtypes::CoxSize order(CoxGraph& G, bits::Lflags I);
  coxtypes::ParSize quotOrder(CoxGraph& G, bits::Lflags I, bits::Lflags J);
  coxtypes::Generator *standardEnumeration(CoxGraph& G, bits::Lflags I);
  const type::Type& type(CoxGraph& G, bits::Lflags I);
};


/* type definitions */

class graph::CoxGraph
{
 private:
  type::Type d_type;
  coxtypes::Rank d_rank;
  CoxMatrix d_matrix; // flattened, by row
  bits::Lflags d_S; // all generators
  containers::vector<bits::Lflags> d_star; // neighbourhoods of the vertices
  containers::vector<bits::Lflags> d_finite_edges; // list of non-infinite edges
 public:
/* constructors and destructors */
  CoxGraph(const type::Type& x, const coxtypes::Rank& l);
  ~CoxGraph();
/* accessors */
  bits::Lflags component(bits::Lflags I,coxtypes::Generator s) const;
  bits::Lflags extremities(bits::Lflags I) const;
  CoxEntry M(coxtypes::Generator s, coxtypes::Generator t) const
    { return(d_matrix[s*d_rank + t]); } // this defines the matrix layout
  bits::Lflags nodes(bits::Lflags I) const;
  coxtypes::Rank rank() const { return d_rank; }
  bits::Lflags supp() const { return d_S; }
  bits::Lflags star(coxtypes::Generator s) const { return(d_star[s]); }
  bits::Lflags star(bits::Lflags I, coxtypes::Generator s) const
    { return(d_star[s]&I); }
  const containers::vector<bits::Lflags>& finite_edges() const
    { return d_finite_edges; }
  const type::Type& type() const { return d_type; }

}; // |class graph::CoxGraph|

#endif

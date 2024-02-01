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
  typedef unsigned short CoxEntry;
  typedef list::List<CoxEntry> CoxMatrix;


/* constants */

  const Ulong SBITMAP_MAX =
    coxtypes::RANK_MAX/CHAR_BIT
    + (bool)(coxtypes::RANK_MAX % CHAR_BIT);
  /* a CoxNbr should hold at least 2 COXENTRY_MAX elements */
  static const CoxEntry COXENTRY_MAX = 32763;
  static const CoxEntry undef_coxentry = USHRT_MAX;
  static const CoxEntry infty = 0;

/******** function declarations **********************************************/



  void getConjugacyClasses(list::List<bits::Lflags>& cl, const CoxGraph& G);
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
  CoxMatrix d_matrix;
  bits::Lflags d_S;
  list::List<bits::Lflags> d_star;
  list::List<bits::Lflags> d_starOps;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(CoxGraph));}
  CoxGraph(const type::Type& x, const coxtypes::Rank& l);
  ~CoxGraph();
/* accessors */
  bits::Lflags component(bits::Lflags I,coxtypes::Generator s) const;
  bits::Lflags extremities(bits::Lflags I) const;
  CoxEntry M(coxtypes::Generator s, coxtypes::Generator t) const;                  /* inlined */
  bits::Lflags nodes(bits::Lflags I) const;
  coxtypes::Rank rank() const;                                           /* inlined */
  bits::Lflags supp() const;                                         /* inlined */
  bits::Lflags star(coxtypes::Generator s) const;                              /* inlined */
  bits::Lflags star(bits::Lflags I, coxtypes::Generator s) const;                    /* inlined */
  const list::List<bits::Lflags>& starOps() const;                         /* inlined */
  const type::Type& type() const;                                    /* inlined */
};

/******** inline definitions **********************************************/

namespace graph {

  inline CoxEntry CoxGraph::M(coxtypes::Generator s, coxtypes::Generator t) const
    {return(d_matrix[s*d_rank + t]);}
  inline coxtypes::Rank CoxGraph::rank() const {return d_rank;}
  inline bits::Lflags CoxGraph::supp() const {return d_S;}
  inline bits::Lflags CoxGraph::star(coxtypes::Generator s) const {return(d_star[s]);}
  inline bits::Lflags CoxGraph::star(bits::Lflags I, coxtypes::Generator s) const
    {return(d_star[s]&I);}
  inline const list::List<bits::Lflags>& CoxGraph::starOps() const {return d_starOps;}
  inline const type::Type& CoxGraph::type() const {return d_type;}

};

#endif

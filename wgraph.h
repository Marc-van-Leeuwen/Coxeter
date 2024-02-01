/*
  This is wgraph.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef WGRAPH_H  /* guard against multiple inclusions */
#define WGRAPH_H

#include "globals.h"
#include "list.h"

/******** type declarations *************************************************/

namespace wgraph {
  class WGraph;
  class OrientedGraph;

  typedef Ulong Vertex;
  typedef unsigned short Coeff;
  typedef list::List<Coeff> CoeffList;
  typedef Vertex Edge;
  typedef list::List<Edge> EdgeList;
};

/******** type definitions **************************************************/

#include "bits.h"
#include "interface.h"

class wgraph::OrientedGraph {
 private:
  list::List<EdgeList> d_edge;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(OrientedGraph));}
  OrientedGraph(const Ulong &n):d_edge(n) {};
  ~OrientedGraph();
/* accessors */
  void cells(bits::Partition& pi, OrientedGraph* P = 0) const;
  const EdgeList& edge(const Vertex& x) const;                 /* inlined */
  Vertex firstMinimal(const bits::BitMap& b) const;
  void levelPartition(bits::Partition& pi) const;
  void print(FILE* file) const;
  Ulong size() const;                                          /* inlined */
/* modifiers */
  EdgeList& edge(const Vertex& x);                             /* inlined */
  void permute(const bits::Permutation& a);
  void reset();
  void setSize(const Ulong& n);                                /* inlined */
};

class wgraph::WGraph {
 private:
  OrientedGraph* d_graph;
  list::List<CoeffList> d_coeff;
  list::List<bits::Lflags> d_descent;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(WGraph));}
  WGraph(const Ulong &n);
  ~WGraph();
/* accessors */
  const CoeffList& coeffList(const Vertex& x) const;             /* inlined */
  const bits::Lflags& descent(const Vertex& x) const;            /* inlined */
  const EdgeList& edge(const Vertex& x) const;                   /* inlined */
  const OrientedGraph& graph() const;                            /* inlined */
  Ulong size() const;                                            /* inlined */
/* modifiers */
  CoeffList& coeffList(const Vertex& x);                         /* inlined */
  bits::Lflags& descent(const Vertex& x);                        /* inlined */
  EdgeList& edge(const Vertex& x);                               /* inlined */
  OrientedGraph& graph();                                        /* inlined */
  void reset();
  void setSize(const Ulong& n);
/* input/output */
  void print(FILE* file, const interface::Interface& I) const;
};

namespace wgraph {

  inline const CoeffList& WGraph::coeffList(const Vertex& x) const
    {return d_coeff[x];}
  inline const bits::Lflags& WGraph::descent(const Vertex& x) const
    {return d_descent[x];}
  inline const EdgeList& WGraph::edge(const Vertex& x) const
    {return d_graph->edge(x);}
  inline const OrientedGraph& WGraph::graph() const {return *d_graph;}
  inline CoeffList& WGraph::coeffList(const Vertex& x) {return d_coeff[x];}
  inline EdgeList& WGraph::edge(const Vertex& x) {return d_graph->edge(x);}
  inline Ulong WGraph::size() const {return d_graph->size();}
  inline OrientedGraph& WGraph::graph() {return *d_graph;}

  inline bits::Lflags& WGraph::descent(const Vertex& x) {return d_descent[x];}
  inline const EdgeList& OrientedGraph::edge(const Vertex& x) const
    {return d_edge[x];}
  inline Ulong OrientedGraph::size() const {return d_edge.size();}
  inline EdgeList& OrientedGraph::edge(const Vertex& x) {return d_edge[x];}
  inline void OrientedGraph::setSize(const Ulong& n) {d_edge.setSize(n);}

};

#endif

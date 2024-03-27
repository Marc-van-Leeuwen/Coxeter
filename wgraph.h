/*
  This is wgraph.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef WGRAPH_H  /* guard against multiple inclusions */
#define WGRAPH_H

#include "globals.h"
#include "containers.h"
#include "list.h"

/******** type declarations *************************************************/

namespace wgraph {
  class WGraph;
  class OrientedGraph;

  typedef Ulong Vertex;
  typedef unsigned short Coeff;
  typedef list::List<Coeff> CoeffList;
  typedef Vertex Edge;
  typedef containers::vector<Edge> EdgeList;
};

/******** type definitions **************************************************/

#include "bits.h"
#include "interface.h"

class wgraph::OrientedGraph
{
 private:
  containers::vector<EdgeList> d_edge;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(OrientedGraph));}
  OrientedGraph(const Ulong &n):d_edge(n) {}
/* accessors */
  bits::Partition cells(OrientedGraph* P = nullptr) const;
  const EdgeList& edge(const Vertex& x) const { return d_edge[x]; }
  bits::Partition level_partition() const;
  void print(FILE* file) const;
  Ulong size() const {return d_edge.size(); }
/* modifiers */
  EdgeList& edge(const Vertex& x) { return d_edge[x]; }
  template<typename I> void add_source_vertex(I begin,I end)
  { d_edge.push_back(EdgeList(begin,end)); }
  void permute(const bits::Permutation& a);
  void reset();
  void setSize(const Ulong& n) { d_edge.resize(n); } // no new edges yet
}; // |class wgraph::OrientedGraph|

class wgraph::WGraph
{
 private:
  OrientedGraph* d_graph;
  list::List<CoeffList> d_coeff;
  list::List<GenSet> d_descent;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(WGraph));}
  WGraph(const Ulong &n);
  ~WGraph();
/* accessors */
  const CoeffList& coeffList(const Vertex& x) const;             /* inlined */
  const GenSet& descent(const Vertex& x) const;                  /* inlined */
  const EdgeList& edge(const Vertex& x) const;                   /* inlined */
  const OrientedGraph& graph() const { return *d_graph; }
  Ulong size() const;                                            /* inlined */
/* modifiers */
  CoeffList& coeffList(const Vertex& x);                         /* inlined */
  GenSet& descent(const Vertex& x);                              /* inlined */
  EdgeList& edge(const Vertex& x);                               /* inlined */
  OrientedGraph& graph() { return *d_graph; }
  void reset();
  void setSize(const Ulong& n);
/* input/output */
  void print(FILE* file, const interface::Interface& I) const;
}; // |class wgraph::WGraph|

namespace wgraph {

  inline const CoeffList& WGraph::coeffList(const Vertex& x) const
    {return d_coeff[x];}
  inline const GenSet& WGraph::descent(const Vertex& x) const
    {return d_descent[x];}
  inline const EdgeList& WGraph::edge(const Vertex& x) const
    {return d_graph->edge(x);}
  inline CoeffList& WGraph::coeffList(const Vertex& x) {return d_coeff[x];}
  inline EdgeList& WGraph::edge(const Vertex& x) {return d_graph->edge(x);}
  inline Ulong WGraph::size() const {return d_graph->size();}

  inline GenSet& WGraph::descent(const Vertex& x) {return d_descent[x];}

};

#endif

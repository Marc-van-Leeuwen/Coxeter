/*
  This is wgraph.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "wgraph.h"

#include <algorithm>

#include "sl_list.h" // to complete |containers::queue|

/****************************************************************************

  This module contains code for the implementation of W-graphs. Recall that
  for a given Coxeter group W, a W-graph is the datum of a graph, together
  with a labelling of the vertices by subsets of the generating set S of W,
  and a labelling of the (oriented) edges by integers mu(x,y); these data
  have to be such that certain explicit formulae define a representation of
  the Hecke algebra of W on the free Z[q^{1/2},q^{-1/2}]-module generated
  by the vertices of the graph.

  Of course the main source of such graphs (the only one in this program)
  is from (subsets of) the group W itself, where the subsets of S are descent
  sets, and the mu(x,y) come from the k-l polynomials.

  We are interested in the following problems :

    - decompose a given graph in cells;
    - compare graphs up to isomorphism;
    - find the cells in a graph up to isomorphism;
    - check if a graph is a W-graph;

 ****************************************************************************/

namespace wgraph {


/****************************************************************************

        Chapter I -- The WGraph class

  Recall that a W-graph is an oriented graph, together with the datum of a
  subset of the generating set S for each vertex, and a coefficient mu(x,y)
  for each edge, so that certain formulae (abstracted from the formulae
  that give the action of the standard generators of the Hecke algebra on
  the k-l basis) yield a representation of the Hecke algebra. In particular,
  such graphs (left,right,and two-sided) are defined on the full group W,
  and their restriction to left (right,two-sided) cells is again a W-graph.

  In our implementation, we have separated the oriented graph proper from
  the subsets and coefficients; this is reasonable as some important
  computations, such as the determination of the cells, only use the graph
  structure, and moreover are most naturally implemented (and will actually
  be used later on) in the setting of abstract graphs.

  The following functions are defined :

    - constructors and destructors :

      - WGraph(n) : allocates a W-graph with n vertices, no edges;
      - ~WGraph() : destructor;

    - accessors :

      - edge(x) : returns a reference to the list of edges originating from x
                  (inlined);
      - graph() : returns a constant reference to the graph (inlined);
      - size() : returns the number of vertices (inlined);

    - modifiers :

      - graph() : returns a non-constant reference to the graph (inlined);
      - reset() : resets the structure;
      - setSize(const Ulong &n) : sets the size to n;

    - input/output :

      - print(file) : prints the graph on a file in ascii format;

 ****************************************************************************/


WGraph::WGraph(const Ulong& n):d_coeff(n),d_descent(n)

/*
  Constructor for the WGraph class.
*/

{
  d_graph = new OrientedGraph(n);
}

WGraph::~WGraph()

/*
   The only non-automatic part is the deletion of d_graph.
*/

{
  delete d_graph;
}


/*
  This function resets the structure to an empty graph of the same size.
*/
void WGraph::reset()
{
  d_graph->reset();
  d_coeff.setZero();
  d_descent.setZero();
}


/*
  Sets the sizes of the data structures so that the graph can accomodate
  n vertices.
*/
void WGraph::setSize(const Ulong& n)
{
  d_graph->setSize(n);
  d_coeff.setSize(n);
  d_descent.setSize(n);
}


/*
  Prints the graph on a file in ascii format.
*/
void WGraph::print(FILE* file, const interface::Interface& I) const
{
  const OrientedGraph& Y = *d_graph;

  int d = io::digits(size()-1,10);

  /* count number of edges */

  Ulong count = 0;

  for (Vertex x = 0; x < size(); ++x) {
    const EdgeList& e = Y.edge(x);
    count += e.size();
  }

  // find alignement

  std::string str;
  Lflags f = constants::lt_mask[I.rank()];
  interface::append(str,f,I);
  Ulong descent_maxwidth = str.length();

  fprintf(file,"%lu vertices, %lu edges\n\n",size(),count);

  for (Vertex x = 0; x < size(); ++x) {
    fprintf(file,"%*lu : ",d,x);
    str.clear();
    interface::append(str,descent(x),I);
    io::pad(str,descent_maxwidth );
    io::print(file,str);
    fprintf(file," ");
    const EdgeList e = Y.edge(x);
    const CoeffList c = coeffList(x);
    for (Ulong j = 0; j < e.size(); ++j) {
      fprintf(file,"%lu(%lu)",e[j],static_cast<Ulong>(c[j]));
      if (j+1 < e.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"\n");
  }
}

/****************************************************************************

        Chapter II -- The OrientedGraph class.

  An oriented graph is just a set of vertices, and for each vertex, a list
  of vertices representing the edges originating from that vertex.

  The following functions are defined :

    - constructors and destructors :

      - OrientedGraph(n) : initializes a graph with n vertices and no edges
        (inlined);
      - ~OrientedGraph() : destructor;

    - accessors :

      - edge(x): returns a constant reference to the list of edges originating
        from x (inlined);
      - size(): returns the number of vertices (inlined);
      - normalPermuation(a) : gives the permutation to a normalized form;
      - print(): prints out the graph;
      - cells(): compute the partition corresponding to strong components
        of the preorder relation defined by the cell edges (as going to smaller)

    - modifiers :

      - edge(x): returns a non-copnstant reference to the list of edges
        originating from x (inlined);
      - permute(a): permutes the graph according to a;
      - reset(): resets the structure;
      - setSize(n): sets the size of the data to accomodate n vertices
        (inlined);

 ****************************************************************************/


/*
  The |cells| function is a subtle algorithm, due to Tarjan. Fokko's original
  implementation was copied into Atlas (and adapted to its evolving support
  structure). It is however incorrect, although some effort is needed to find a
  graph for which it fails, and apparently none was ever encountered in
  practive. It was therefore rewritten by Marc, which version was then (much)
  later imported back into Coxeter by Marc. We given below both a polished (but
  still incorrect) version of the algorithm that was in Coxeter, with an
  explanation of why it is incorrect, and then a version imported from Atlas
  (with a few improvements).
*/


#if 0 // go to #else: use new implementation of Tarjan's algorithm down below

namespace {

  void get_class(const OrientedGraph& X, const Vertex& y, bitmap::BitMap& b,
		 bits::Partition& pi, OrientedGraph* P);
}; // |namespace|

/*
  Define a preorder relation on the vertices by setting $x\leq y$ if there is an
  oriented path from $y$ to $x$. This function returns the partition into
  equivalence classes (strong components) of this preorder. We use the Tarjan
  algorithm, explained in one of the Knuth books, but which I learned from Bill
  Casselman. [dixit Fokko; I have increasingly adapted/added what follows MvL]

  The vertices for which the partition function is already defined will be said
  to be settled; a bitmap |removed| indicated those elements. Some other
  attributes are attributed to each vertex via local arrays. Each major round of
  computation starts with a vertex that has not been seen in previous rounds,
  say x0. In this round we shall visit all vertices visible from x0 by forward
  edges. We do so by a depth-first traversal that avoids re-traversing vertices
  that are already seen during this traversal, while keeping a stack of vertices
  along the first path found leading to the current node. Data is propagated
  back up during the search, in such a manner that when we have completed the
  traversal of descendants of a node |x|, we can decide whether any of its
  predecessors on the stack can be reached from it. If this is not the case, |x|
  heads a strong component, which can be find as everything reachable from |x|
  that is not already removed; we find and remove those vertices in a separate
  traversal in |get_class|. In any case we then continue the depth-first
  traversal by resuming the pass over the children of the parent of |x|.

  The extra data consists of some values that would have been held in local
  variables of the function performing the depth first search if that were done
  recursively, represented by the stack |active| in our version below, in which
  the recustion has been transformed into iteration; in addition there is one
  numeric attribute |min| that is maintained for all vertices, regardless if they
  have an |active| record or not (but for vertices that become |removed|, this
  value will no longer be looked at). The "local" data on the stack are: the
  number of the vertex for which this record was created, and a pointer to its
  list of edges together with a counter of our progress in traversing this list;
  these allow a straightforward implemention of depth-first traversal.

  The |min| field of a vertex serves to recognize that it has been encountered,
  as its initial high value is replaced by the depth of the stack at which its
  record is placed when this happens. However it then gets lowered if some edge
  leads to an already encountered vertex with a lower |min| value, and upon
  completion of the traversal of edges of |x|, its final |min| value gets passed
  back up to its parent, where the minimum with the value stored for the parent
  is agin taken. Just before returning to the parent, the value |min[x]| is
  compared to its initial (stack depth) value; if it is unchanged, then |x|
  heads a strong component which is split off and marked in |removed| before
  returning to the parent (also we skip comparison with the parent |min| value
  in this case, as it will certainly be less than |min[x]|).

  While it is tempting to think of |min[x]| as pointing to the oldest ancestor
  of |x| reachable from |x| (Fokko said as much), this is not true; for one, the
  depth recorded might come from a vertex that no longer |active| and therefore
  not an ancestor of |x|, and also that vertex might have lowered its |min|
  value after it was recorded and passed up the chain. As a consequence, the
  decision to split off a component might be premature; here is an exemple. Tak
  a graph on 7 vertices numbered 0 to 6, with outgoing edges: 0->{4}, 1->{2},
  2->{1,4}, 3->{2,5,6}, 4->{3}, 5->{1} 6->{}. The algoithm pushes in order the
  vertices 0, 4, 3, 2, 1 with their min values increasing from 0 to 4; then 1
  sees 2 and sets min[1]=3 (no component), returning to 2 it sees 4, setting
  min[2]=1 (no component), returning to 3 sees a fresh vertex 5, pushing it at
  level 3 so min[5]=3; then vertex 5 sees 1 with (also) min[1]==3 so min[5]==3
  remains, and upon completion a (supposedly singleton) component headed by 5 is
  falsely detected. In reality it heads a cycle 5->1->2->4->3->5, which should
  not be detected as component before the sink 6 is removed as singleton
  component, while at this point that vertex has not even been visited.
*/
bits::Partition OrientedGraph::cells(OrientedGraph* P) const // by incorrect Tarjan
{
  struct stack_record { Vertex v; const EdgeList* elist; Ulong ecount; };
  containers::vector<stack_record> active; // like a |stack|, but we need |size()|

  struct vertex_record { Ulong rank; Ulong min; };

  // area of minimal stack depth of descendendent of vertex encountered
  containers::vector<Ulong> min(size(),size()); // length and value both |size()|

  bits::Partition pi(size());

  bitmap::BitMap removed(size());

  for (Vertex x0 = 0; x0 < size(); ++x0)
  {
    if (removed.is_member(x0)) // |x0| was dealt with previously
      continue;

    assert(active.empty());
    active.push_back(stack_record{0,&edge(x0),0}); // initial stack record

    while(not active.empty())
    {
    continue_while:
      auto& cur = active.back();
      Vertex x = cur.v;

      for (const EdgeList& e = *cur.elist; cur.ecount < e.size(); ++cur.ecount)
      {
	Vertex y = e[cur.ecount];
	if (removed.is_member(y))
	  continue;
	if (min[y] == size()) // whether |y| is new
	{ // push a new vertex record
	  min[y] = active.size(); // level at which next record will be stored
	  active.push_back(stack_record{y,&edge(y),0});
	  goto continue_while; // C++ does not let us |continue| a non-inner loop
	}
	else if (min[x] > min[y])
	  min[x] = min[y]; // take into account low |min| value encountered
      } // |for(cur.ecount)|
      // at this point we have exhausted the edges of |x|
      active.pop_back(); // we still have |x|
      if (min[x] == active.size())
	get_class(*this,x,removed,pi,P); // record and remove a class
      else if (min[x] < min[active.back().v]) /* if t=1, previous case holds */
	min[active.back().v] = min[x]; // propagate low value to parent
    } // |while(t>0)|

  }

  return pi;
} // |OrientedGraph::cells|


namespace {


/*
  This is an auxiliary to (strong component) |OrientedGraph::cells|.

  After the element y has been identified as minimal among the elements not
  already marked in b, this function marks off the equivalence class of y;
  these are just the elements visible from y and not already marked in b.
  The class is written as a new class into |pi|, and induced edges of the
  quotient graph on strong components are added to |P|, if it is non-null.

  This function modifies |pi| by direct access to its components; that should be
  avoided, but this function will go away when Tarjan gets inported from Atlas.
*/
void get_class(const OrientedGraph& X, const Vertex& x0, bitmap::BitMap& recorded,
	       bits::Partition& pi, OrientedGraph* P)
{
  const Ulong cc = pi.classCount(); // prepare for adding a new class to |pi|
  pi[x0] = cc;

  containers::queue<Ulong> q { x0 };
  recorded.insert(x0);
  std::set<Ulong> dest_classes;
  do // |not q.empty()|
  {
    Vertex x = q.front();
    q.pop();
    for (Vertex y : X.edge(x))
    {
      if (recorded.is_member(y)) // take note of edges to earlier components
      {
	if (P!=nullptr and pi[y] < cc) // then possibly add a new destination class
	  dest_classes.insert(pi[y]);
      }
      else // |y| not in |recorded|
      {
	q.push(y);
	recorded.insert(y);
	pi[y] = cc;
      }
    } // |for(y)|
  } while (not q.empty());

  pi.incr_class_count();
  if (P!=nullptr)
    P->add_source_vertex(dest_classes.begin(),dest_classes.end());
} // |get_class|

}; // |namespace|


#else // use new implementation
/*
  Define a preorder relation on the vertices by setting x <= y iff there is an
  oriented path from x to y. This function puts in pi the partition function
  corresponding to the equivalence classes of this preorder. We use the
  Tarjan algorithm, explained in one of the Knuth books, but which I learned
  from Bill Casselman. [Fokko; I have increasingly adapted/added what follows MvL]

  The vertices for which the partition function is already defined will be said
  to be settled. That status and some other attributes will be attributed to
  each vertex via local arrays. Start each major round with a vertex that has
  not been seen in previous rounds, say x0. In this round we shall visit all
  vertices visible from x0 by forward edges. We do so by a depth-first
  traversal, keeping a stack of vertices along the first path found leading to
  the current node. Data is propagated back up during the search, in such a
  manner that when we have completed the traversal of descendants of a node |x|,
  we can decide whether any of its predecessors on the stack can be reached from
  it. If this is the case we back up to the parent of |x| but leave |x| in a
  temporarily sidelined part of the stack: when growing (for other descendants
  of the parent of |x| for instance) the new items will be added the same place
  the stack top was before (somewhere above |x|). The idea is that |x| and any
  vertices left from the traversal of its descendants will be scooped up when a
  predecessor of |x| later decides that it is the head of a strong component: at
  that time the vertices of that component will for an upper part of the stack,
  which is then amputated as a whole. This action is precisely what is done in
  the case we left aside in our discussion, namely when the data gathered allow
  concluding that none of the predecessors of |x| can be reached from it.

  This extra data consists of some values that would have been held in local
  variables of the function performing the depth first search if that were done
  recursively, but which will be held in vector of records for active nodes in
  the iterative version, and one numeric attribute that is added to all vertices
  (and which therefore can be inspected at any time, whether or not the vertex
  if the currently visited one). The latter is simply a sequence number
  attributed when the vertex is first visited, with a special value for vertices
  not yet visited, and other (higher) values attributed to their vertices as
  each strong component gets split off (we can use different values for each
  component to represent the partition into strong components). The other
  "local" data are: the number of the vertex for which this record was created,
  a reference to the record of its parent so we can continue our search after
  having exhausted its edges, an indication of our progress in traversing its
  edges (so we can continue after having exhausted an edge) and finally a number
  |min| that is initially the sequence number of the vertex, but which can get
  lowered at the completion of the traversal of a child to the least (oldest)
  sequence number seen during that traversal when it encounters already active
  vertices (in which case the search does not make that vertex current, but its
  vertex number is propagated back to |min| values of predecessors as long as it
  is lower than the previously recorded value).

  It is tempting to think of |min| for a vertex as the minimal sequence number
  reachable from the vertex at which the value is stored [Fokko said as much],
  but this is not true: we only record the vertex sequence number for an already
  active vertex itself, which does not take into account anything reachable from
  that vertex, and even if we did record the current minimum for that vertex, it
  might not yet reflect all vertices reachable from it. However it is true that
  a vertex heads a component if and only if the |min| value when completing the
  traversal of all edges is still at its initial value, and this it precisely
  what |min| is used for. That this decision is justified, as well as the
  consequent action of making into a component all still active vertices from
  the current vertex onward, is based on two invariants. (1) all active vertices
  newer than the current vertex $c$ can be reached from $c$, and (2) from any
  (still) active vertex, $c$ can be reached. The former is clear from the
  depth-first search description, the second is because if $a$ is not an
  ancestor of $c$, then its traversal is complete, and since it was not removed,
  some younger vertex $a'$ can be reached from it. Then if $a'$ is still no
  ancestor of $c$, the same applies to $a'$, and so forth until reaching an
  ancestor of $c$. Applying the two statements at the point where the traversal
  of $c$ is complete, we see that $c$ heads a strong component iff the traversal
  saw no vertices younger than $c$, and if so its consists of all still active
  vertices not younger than $c$.
*/

bits::Partition OrientedGraph::cells(OrientedGraph* gr) const // Tarjan's algorithm
{
  using seqno = Vertex; // sequence number in depth-first traversal
  using work_addr = unsigned; // reference to an active vertex by stack location

  struct info
  {
    Vertex v;    // identification of vertex in the graph
    work_addr parent; // location of parent in depth-first traversal
    unsigned next_edge; // index in edge list of next edge to consider
    seqno min;        // minimal rank of vertex reachable as indicated above

    // constructor
    info (const Vertex x, work_addr p,
	  containers::vector<seqno>& rank, seqno& counter)
      : v(x), parent(p), next_edge(0), min(rank[x]=counter++) { }
  };

  containers::vector<seqno> rank(size(),0);

  const work_addr nil= size(); // impossible index into |active|
  const seqno limit = size()+1; // ranks from here up are classified vertices


  containers::vector<info> active;
  Ulong comp_count=0; // incremented after each component is split off

  // if (gr!=nullptr) gr->resize(0); // start induced graph with a clean slate

  /* the following loop guarantees that all strong components will be found
     regardless of where in the graph we start, and whether or not it is
     connected
  */
  for (Vertex x0 = 0; x0 <size(); ++x0)
    if (rank[x0]<limit)
    {
      seqno count=1;
      active.push_back(info(x0,nil,rank,count)); // x0 has no parent
      work_addr cur_pos=0; // point current position to x0

      while(cur_pos!=nil)
      {
	info& x = active[cur_pos];
	const EdgeList& edges = edge(x.v);

	if (x.next_edge < edges.size())
	{
	  Vertex y = edges[x.next_edge++];
	  if (rank[y]==0) // y is a fresh vertex
	  {
	    work_addr y_pos = active.size();
	    active.push_back(info(y,cur_pos,rank,count));
	    cur_pos = y_pos; // and |continue| the |while| loop from there
	  }
	  else // |y| was seen before (cross edge), or |y| is settled
	  {  // if |y| is settled nothing will happen
	    if (rank[y] < x.min) // then record that we can reach y
	      x.min = rank[y];
	    // now |continue| the while loop without changing |cur_pos|
	  }
	}
	else // we have exhausted the edges of |x.v|, so it matures
	  if (x.min == rank[x.v]) // no older vertex reachable from |x|
	  { // split off strong component
	    const auto here = active.begin()+cur_pos;
	    for (auto it = here; it!=active.end(); ++it)
	      rank[it->v] = limit+comp_count; // settle vertex, class |comp_count|
	    if (gr!=nullptr)
	    {
	      std::set<Vertex> out_edges;
	      for (auto it = here; it!=active.end(); ++it)
		for (Vertex v : edge(it->v))
		  if (rank[v] != limit+comp_count) // avoid induced self-edges
		    out_edges.insert(rank[v]-limit);
	      gr -> add_source_vertex(out_edges.begin(),out_edges.end());
	    }
	    active.erase(here,active.end()); // amputate no longer |active| part
	    ++comp_count;

	    cur_pos = x.parent; // we shall resume the traversal there
	  }
	  else // |x| matures but does not head a new strong component
	  { // note that |x| cannot be |x0|, so active[x.parent] exists
	    cur_pos = x.parent; // we shall continue the traversal there
	    if (x.min < active[cur_pos].min) // then update parent info
	      active[cur_pos].min=x.min; // what |x| sees, its parent sees
	  }
      }  // while(cur_pos!=nil)

      assert(active.empty());

    } //for (x0) if (rank[x0]<infinity)

  for (auto& entry : rank)
  { assert(entry>=limit); // node was incorporated
    entry -= limit; // this will be the actual classifier value
  }

  bits::Partition pi(std::move(rank),comp_count); //
  return pi;
}


#endif

/*
  Assuming the graph has no oriented cycles, this function writes in pi the
  partition of the vertices according to their level, where sinks have level
  0, then sinks in the remaining poset have level one, etc.

  NOTE : the implementation is simple-minded : we traverse the graph as many
  times as there are levels. [This method is never used.]
*/
bits::Partition OrientedGraph::level_partition() const
{
  containers::vector<Ulong> classify(size(),Ulong(-1));
  Ulong count = 0; // count elements encountered
  Ulong cur_level=0;

  do // run at least once, since |size()>0|
  {
    for (bits::SetElt x = 0; x < size(); ++x)
    {
      if (classify[x]<cur_level)
	continue;
      const EdgeList e = d_edge[x];
      for (Ulong j = 0; j < e.size(); ++j)
	if (classify[e[j]]>=cur_level)
	  goto next_x; // |x| not a sink when removing lower levels

      // if we get here, |x| is at the current level
      classify[x] = cur_level;
      ++count; // keep track of number of elements marked, for termination
    next_x:
      continue;
    } // |for(x)|
    ++cur_level; // this level is complete
  } while (count < size());

  return bits::Partition(std::move(classify),cur_level);
} // |levelPartition|


/*
  Permute the graph according to the permutation a, according
  to the usual rule : the edges of a(x) should be the image under a of the
  edge set of x.

  As usual, permuting values is easy : it is enough to apply a to the
  elements in the various edgelists. Permuting ranges is trickier, because
  it involves a^-1.

  It is assumed of course that |a| holds a permutation of size |size()|.
*/
void OrientedGraph::permute(const bits::Permutation& a)
{
  static bits::BitMap b(0);
  static EdgeList e_buf(0);

  /* permute values */

  for (bits::SetElt x = 0; x < size(); ++x) {
    EdgeList& e = d_edge[x];
    for (Ulong j = 0; j < e.size(); ++j) {
      e[j] = a[e[j]];
    }
  }

  /* permute ranges */

  b.setSize(size());
  b.reset();

  for (bits::SetElt x = 0; x < size(); ++x) {
    if (b.getBit(x))
      continue;
    if (a[x] == x) { /* fixed point */
      b.setBit(x);
      continue;
    }
    for (bits::SetElt y = a[x]; y != x; y = a[y])
      std::swap(d_edge[x],d_edge[y]);
    b.setBit(x);
  }
} // |OrientedGraph::permute|


// Print out the graph on |file|
void OrientedGraph::print(FILE* file) const
{
  fprintf(file,"size : %lu\n\n",size());

  int d = io::digits(size(),10);

  for (Vertex x = 0; x < size(); ++x) {
    const EdgeList& e = edge(x);
    fprintf(file,"%*lu : ",d,x);
    for (Ulong j = 0; j < e.size(); ++j) {
      fprintf(file,"%*lu",d,e[j]);
      if (j < e.size()-1) { /* there is more to come */
	fprintf(file,",");
      }
    }
    fprintf(file,"\n");
  }

  fprintf(file,"\n");
} // |print|


// Reset the structure to hold an edge-less graph of the same size.
void OrientedGraph::reset()
{
  for (Ulong j = 0; j < size(); ++j) {
    d_edge[j].clear();
  }
}


/****************************************************************************

        Chapter III -- Auxiliaries

  This chapter contains some auxiliary functions for the main functions in
  this module. The following functions are defined :

    - get_class(X,y,b) : gets the class of y in X, using the bitmap b;

 ****************************************************************************/
}; // |namespace wgraph|

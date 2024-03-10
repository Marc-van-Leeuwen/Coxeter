/*
  This is schubert.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "schubert.h"

#include "error.h"

namespace schubert {
  using namespace error;
};

/****************************************************************************

  This file contains the definition of the general functions pertaining
  to Schubert contexts, and the definition of the SchubertContext
  class. Other implementations of Schubert contexts will not be found in
  other files, maybe they were meant to be wriiten

  A Schubert context contains a description of a finite downwards closed subset
  Q of the Coxeter group (initially the singleton {e}), and should be capable of
  providing data and services related to the Bruhat ordering and to the
  (partially defined) action of the generators on the Q. The elements of Q are
  assumed to be enumerated by the integers in the range [0,size()[, in an
  ordering compatible with the Bruhat ordering.

  Specifically, the following functions are required :

  - managing the context :

    - append(g,x) : appends to the CoxWord g a reduced expression of the
      element #x in Q;
    - contextNumber(g) : returns the number x in [0,size()[ corresponding
      to the CoxWord g, or undef_coxnbr if there is no such x;

  - descent sets :

    - descent(x) : returns the two-sided descent set of x, flagged in a
      Lflags;
    - ldescent(x), rdescent(x) : same for left (right) descent sets;
    - firstDescent(x), firstLDescent(); firstRDescent(x) : returns the
      first set bit in the corresponding descent sets;
    - downset(s) : a bitmap of the set of x s.t. xs < x;
    - isDescent(x,s) : tells if s is a descent for x;
    - maximize(x,f), minimize(x,f) : extremalizes x w.r.t. the action of the
      generators flagged by f;

  - Bruhat ordering :

    - extendSubSet(q,s) : given q holding a decreasing subset of Q, and
      s s.t. q.s is contained in Q, puts q.s in q;
    - extractClosure(q,x) : puts into q the interval [e,x];
    - hasse(x) : returns the coatom list of x;
    - inOrder(x,y) : tells if x <= y;

  - action of generators :

    - shift(x,s) : returns the shift xs (left shift by s-rank() if s > rank());
    - lshift(x,s), rshift(x,s) : left or right shifts;

  - length :

    - length(x) : returns the length of element #x;
    - maxlength() : largest length of element in context;
    - parity(x) : a bitmap of the set of z s.t. length(z) = length(x) mod 2;

  - sizes :

    - rank() : returns the rank of the underlying group;
    - setSize(n) : this should really be a private function; it can be used
      safely only in order to _reduce_ the size of the context (typically,
      when the Schubert extension succeeded but the kl extension failed,
      and we wish to revert;)
    - size() : returns the current size of the context;

  The SchubertContext class contains lengths, shifts, coatoms and
  downsets as tables, so the above functions are really simply table accesses.
  The context extension and interval extraction is done as described in my
  paper "Computing Kazhdan-Lusztig polynomials for arbitrary Coxeter groups",
  Experiment. Math. 11 (2002), no. 3, pp. 371-381.

 ****************************************************************************/

namespace {

  using namespace schubert;
  using namespace bits;

  const char* undef_str = "undefined";

};


class SchubertContext::ContextExtension
{
private:
  SchubertContext& d_schubert;
  Ulong d_size;
  coxtypes::CoxNbr* d_star;
public:
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
  {return memory::arena().free(ptr,sizeof(ContextExtension));}
  ContextExtension(SchubertContext& p, const Ulong& c);
  ~ContextExtension();
  Ulong size() {return d_size;}
}; // |class ContextExtension|

/****************************************************************************

        Chapter I -- The SchubertContext class.

  This section defines the functions in the AbstractSchubertContext class :

  - constructors and destructors :

    - SchubertContext(G);
    - ~SchubertContext

  - accessors :

    - append(g,x) : appends the normal form of x to g;
    - closure(x) : returns a bitmap of the interval [e,x];
    - contextNumber(g) : returns the number of g in the context;
    - descent(x), ldescent(x), rdescent(x) : descent sets; (inlined)
    - downset(s) : bitmap of {x, xs < x}; (inlined)
    - firstDescent(x), firstLDescent(x), firstRDescent(x) : smallest
      generator taking x down; (inlined)
    - hasse(x) : the list of coatoms of x; (inlined)
    - inOrder(x,y) : tells if x <= y;
    - length (x) : the length of x; (inlined)
    - maximize(x,f) : maximizes x w.r.t. the action of the generators in f;
    - maxlength() : maximal length in context; (inlined)
    - minimize(x,f) : minimizes x w.r.t. the action of the generators in f;
    - normalForm(g,x,order) : returns the normal form of x for order;
    - parity(x) : bitmap of {z, l(z)=l(x) mod 2}; (inlined)
    - rank() : the rank of the group; (inlined)
    - shift(x,s), lshift(x,s), rshift(x,s) : shifts x by s; (inlined)
    - size() : the size of the context; (inlined)
    - twoDescent(y) : returns the "super-descent" set of y;

  - modifiers :

    - extendContext(g) : extends the context to accomodate g;
    - extendSubSet(s) : extends the subset to accomodate multiplication by s;
    - permute(q) : applies the permutation q to the context;
    - revertSize(n) : reverts to a previous size;
    - setSize(n) : resets the size;
    - subset() : gives access to the subset; (inlined)


 ****************************************************************************/

namespace schubert {

/******** constructors ******************************************************/


/*
  Constructor for the SchubertContext class. The data are initialized for
  the one-element context {e} : size is one, and for the single element
  0, length is 0, descent sets are empty, coatom set is empty, shifts
  are all undefined.
*/
SchubertContext::SchubertContext(const graph::CoxGraph& G)
  : d_graph(G)
  , d_rank(G.rank())
  , d_maxlength(0)
  , d_size(1)
  , d_length(1,0)
  , d_hasse(1,CoxNbrList())
  , d_descent(1,Lflags())
  , d_shift(1,2*d_rank,coxtypes::undef_coxnbr)
  , d_star(1)
  , d_downset(2*d_rank,bitmap::BitMap(1))
  , d_parity{bitmap::BitMap(1),bitmap::BitMap(1)}
  , d_history()
{
  d_star.setSizeValue(1);

  d_star[0] = new(memory::arena()) coxtypes::CoxNbr[2*nStarOps()];
  for (coxtypes::StarOp j = 0; j < 2*nStarOps(); ++j)
    d_star[0][j] = coxtypes::undef_coxnbr;

  d_parity[0].insert(0);
}


/*
  Destructing a SchubertContext turns out to be a little bit tricky, because of
  the way memory has been allocated for the various structures, in the
  successive extensions. In particular, in d_star, only *some* pointers have
  been allocated by new.

  This problem will be much easier to handle once the memory allocations
  go through a private arena; it will then be enough to simply free the
  arena. For now, we introduce a monitoring through a stack of
  ContextExtensions.
*/
SchubertContext::~SchubertContext()
{
  /* reverse history */

  memory::arena().free(d_star[0],2*nStarOps()*sizeof(coxtypes::CoxNbr));
}

/******** accessors ********************************************************/


/*
  Appends to |g| the ShortLex normal form of x. The normal form is
  easily obtained using the left descent sets.

  NOTE : it is the programmer's responsibilty when using this function, to
  guarantee that the result is reduced. Otherwise, use "prod".
*/

coxtypes::CoxWord& SchubertContext::append
  (coxtypes::CoxWord& g, coxtypes::CoxNbr d_x) const
{
  coxtypes::CoxNbr x = d_x;

  while (x>0) {
    coxtypes::Generator s = firstLDescent(x);
    g.append(s+1);
    x = lshift(x,s);
  }

  return g;
}

coxtypes::CoxNbr SchubertContext::contextNumber(const coxtypes::CoxWord& g) const

/*
  This functions returns the number corresponding to g in the current
  context; returns coxtypes::undef_coxnbr if g is not in the context.
*/

{
  coxtypes::CoxNbr x = 0;

  for (Ulong j = 0; j < g.length(); ++j) {
    coxtypes::Generator s = g[j]-1;
    x = rshift(x,s);
    if (x == coxtypes::undef_coxnbr)
      break;
  }

  return x;
}


/*
  Put into |b| the subset $[e,x]$ of |p|. It is assumed that |b|
  is capable of holding a subset of |p|.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.

*/
void SchubertContext::extractClosure
  (bits::BitMap& b, coxtypes::CoxNbr x) const
{
  bits::SubSet q(d_size);
  q.reset();
  q.add(0);

  for (coxtypes::CoxNbr x1 = x; x1;) {
    coxtypes::Generator s = firstLDescent(x1);
    extendSubSet(q,s);
    x1 = lshift(x1,s); // remove leftmost generator
  }

  b = q.bitMap();

  return;
}

void SchubertContext::spread_subset
  (bitmap::BitMap& q, coxtypes::Generator s) const
{
  auto q0 = q; // fix original subset to loop over
  for (coxtypes::CoxNbr x : q0)
    if (not isDescent(x,s)) // unnecessary, if |q| is downwards closed
      q.insert(shift(x,s)); // could be left or right shift (upwards)
}

containers::sl_list<coxtypes::Generator>
  SchubertContext::word (coxtypes::CoxNbr x) const
{
  containers::sl_list<coxtypes::Generator> result;
  while (x>0)
  {
    coxtypes::Generator s = firstLDescent(x);
    result.push_back(s);
    x = lshift(x,s);
  }
  return result;
}

bitmap::BitMap SchubertContext::closure(coxtypes::CoxNbr x) const
{
  assert(x<size());
  bitmap::BitMap result(size()); // full size probably needed by callers
  CoxNbrList elements; // copy for faster iteration
  elements.reserve(x+1); // ensure, vitally, no reallocation will occur
  result.insert(0);
  elements.push_back(0);
  for (auto s : word(x))
  { const auto high_level = elements.end();
    for (auto it = elements.begin(); it!=high_level; ++it)
    {
      const coxtypes::CoxNbr cand = rshift(*it,s);
      if (not result.is_member(cand)) // actually only excludes old elements
      {
	result.insert(cand); // may or may not be new
	elements.push_back(cand);
      }
    }
  }

  return result;
}


/*
  Checks if x <= y in the Bruhat ordering, using the well-known recursive
  algorithm : if y = 0, then x must be 0; otherwise, take s s.t. ys < y;
  then if xs < x, x <= y iff xs < ys; if xs > x, x <= y iff x <= ys.

  NOTE : this function is not intended for heavy use!
*/
bool SchubertContext::inOrder(coxtypes::CoxNbr x, coxtypes::CoxNbr y) const
{
  if (x == 0)
    return true;
  if (x == y)
    return true;
  if (x > y)
    return false;

  coxtypes::Generator s = firstDescent(y); // in practice a right descent

  coxtypes::CoxNbr xs = shift(x,s);
  coxtypes::CoxNbr ys = shift(y,s);

  if (xs < x)
    return inOrder(xs,ys);
  else /* xs > x */
    return inOrder(x,ys);
}


/*
  This function maximizes x w.r.t. the flags in f. The return value is
  coxtypes::undef_coxnbr if the extremalization takes us outside the context. It
  is assumed that f is a valid set of flags, i.e., contained in S \coprod S.
*/
coxtypes::CoxNbr SchubertContext::maximize
  (coxtypes::CoxNbr x, const Lflags& f)  const
{
  coxtypes::CoxNbr x1 = x;
  Lflags g = f & ~d_descent[x1];

  while (g!=0)
    {
      coxtypes::Generator s = constants::firstBit(g);
      x1 = shift(x1,s);
      if (x1 == coxtypes::undef_coxnbr)
	break;
      g = f & ~d_descent[x1];
    }

  return x1;
}


/*
  This function minimizes x w.r.t. the flags in f. Here the return value
  is always defined. It is assumed that f is a valid set of flags, i.e.,
  contained in S \coprod S.
*/
coxtypes::CoxNbr SchubertContext::minimize
  (coxtypes::CoxNbr x, const Lflags& f) const
{
  coxtypes::CoxNbr x1 = x;
  Lflags g = f & d_descent[x1];

  while (g!=0)
    {
      coxtypes::Generator s = constants::firstBit(g);
      x1 = shift(x1,s);
      g = f & d_descent[x1];
    }

  return x1;
}


/*
  This function returns the normal form of x for the given ordering of the
  generators. The |order| argument is typically |d_out| of the interface; so
  order[j] is the external number of the internal generator #j.

  NOTE : this function is more expensive than append; its main intention
  is for use in i/o functions.
*/
coxtypes::CoxWord& SchubertContext::normalForm
  (coxtypes::CoxWord& g, coxtypes::CoxNbr d_x,
   const bits::Permutation& order) const
{
  g.reset();
  coxtypes::CoxNbr x = d_x;

  while (x>0) {
    coxtypes::Generator s = minDescent(ldescent(x),order);
    g.append(s+1);
    x = lshift(x,s);
  }

  return g;
}


/*
  Returns the "super-descent" set of x; this is the union of the descent
  set of x, and of the descent sets of the xs, where s runs through the
  descent set of x.
*/
Lflags SchubertContext::twoDescent(coxtypes::CoxNbr x) const
{
  Lflags f = descent(x);

  for (Lflags f1 = f; f1; f1 &= f1-1)
  {
    coxtypes::Generator s = constants::firstBit(f1);
    coxtypes::CoxNbr xs = shift(x,s);
    f |= descent(xs);
  }

  return f;
}

/******** modifiers ********************************************************/


/*
  This function extends the context to the smallest one containing the
  exixting one and the given g, i.e., the new context is the union of
  the old context and [e,g]. Apart from some previously undefined shifts
  becoming defined, this doesn't induce _any_ modification in the data
  for the old context; the numbers of the new elements come in at the top.

  Sets the error |EXTENSION_FAIL| in case of failure.

  The outline of the function is as follows. First, we determine the largest
  initial subword h of g which is alreaady in the context, and construct the
  interval [e,h] as a subset of the context. Then, for each remaining generator
  s in g, we add the elements in [e,hs] not already in the context, and we
  update everything.
*/
coxtypes::CoxNbr SchubertContext::extendContext
  (const coxtypes::CoxWord& g)
{
  coxtypes::CoxNbr y = 0;
  bits::SubSet q(d_size);
  q.reset();
  q.add(0);

  Ulong j = 0;

  CATCH_MEMORY_OVERFLOW = true;

  for (; j < g.length(); ++j) {
    coxtypes::Generator s = g[j]-1;
    if (rshift(y,s) == coxtypes::undef_coxnbr)
      break;
    extendSubSet(q,s);
    if (ERRNO)
      goto error_handling;
    y = rshift(y,s);
  }

  for (; j < g.length(); ++j) {
    coxtypes::Generator s = g[j]-1;
    fullExtension(q,s);
    if (ERRNO)
      goto error_handling;
    if (j >= d_maxlength)
      d_maxlength = j+1; // keep the high water mark
    y = rshift(y,s);
  }

  CATCH_MEMORY_OVERFLOW = false;

  return y;

 error_handling:
  Error(ERRNO);
  ERRNO = EXTENSION_FAIL;
  return(coxtypes::undef_coxnbr);
}



/*
  Given a subset q of p holding a decreasing subset, and a geneator s s.t.
  q.s. is contained in the context, this function puts in q
  the set q.s ( here s can be either a right or a left shift.)

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.
*/
void SchubertContext::extendSubSet
  (bits::SubSet& q, coxtypes::Generator s) const
{
  Ulong a = q.size();

  for (Ulong j = 0; j < a; ++j) { /* run through q */
    coxtypes::CoxNbr x = (coxtypes::CoxNbr)q[j];
    coxtypes::CoxNbr xs = shift(x,s);
    if (xs < x)
      continue;
    if (q.isMember(xs))
      continue;
    /* if we get here a new element is found */
    q.add(xs);
    if (ERRNO)
      return;
  }

  return;
}


/*
  Apply the permutation |a| to the context, meaning that for each |i| the
  Coxeter element |i| will henceforth have number |a[i]|. We have explained in
  kl.cpp how this should be done. The objects to be modified are the following :

   - d_length : a table with range in the context;
   - d_hasse : each row is a list with values in the context; the table itself
     has range in the context; in addition the permuted rows should be sorted;
   - d_descent : a table with range in the context;
   - d_shift : a table with range in the context; each row has values in the
     context, or coxtypes::undef_coxnbr;
   - d_downset : a table of bitmaps ranging over the context;
   - d_parity : a pair of bitmaps ranging over the context;
*/
void SchubertContext::permute(const bits::Permutation& a)
{

  /* permute values */

  for (CoxNbrList& c : d_hasse)
    for (auto& elt : c)
      elt = a[elt];

  for (coxtypes::CoxNbr x = 0; x < d_size; ++x)
    for (coxtypes::Generator s = 0; s < 2*d_rank; ++s)
    { auto& xs = d_shift.entry(x,s);
      if (xs != coxtypes::undef_coxnbr)
	xs = a[xs];
    }

  /* permute the ranges */
  bitmap::BitMap seen(a.size());

  for (coxtypes::CoxNbr x = 0; x < this->size(); ++x)
  {
    if (seen.is_member(x))
      continue;
    if (a[x] == x) {
      seen.insert(x); // not really useful, we should never come here again
      continue;
    }

    for (coxtypes::CoxNbr y = a[x]; y != x; y = a[y])
    {
      std::swap(d_length[x],d_length[y]);
      std::swap(d_hasse[x],d_hasse[y]);
      d_shift.swap_rows(x,y);

      /* back up values for y */

      Lflags descent_buf = d_descent[y];

      /* put values for x in y */

      d_descent[y] = d_descent[x];

      /* store backup values in x */

      d_descent[x] = descent_buf;

      // swap bits at |x| and |y| in all downsets
      for (coxtypes::Generator s = 0; s < 2*this->rank(); ++s)
      {
	bool t = d_downset[s].is_member(y);
	d_downset[s].set_to(y,d_downset[s].is_member(x));
	d_downset[s].set_to(x,t);
      }

      /* modify parity bitmaps */

      for (unsigned i : {0,1})
      {
	bool t = d_parity[i].is_member(y);
	d_parity[i].set_to(y,d_parity[i].is_member(x));
	d_parity[i].set_to(x,t);
      }

      seen.insert(y); // mark |y| as seen|
    }
    seen.insert(x); // mark |x| as seen|
  }

}


/*
  Revert the size of the context to some previous value. It is very important
  that |n| is indeed a previous value of the context size (i.e. that it can be
  found through the extension history), and that no permutation has taken place
  between that previous size and the current size. This function is intended to
  be used immediately after the extension of the rest of the context failed, so
  that everything returns to where it was before. Because of the way things are
  allocated, we are actually able to return most of the allocated memory (only
  the allocations for the basic lists will remain as they were.)
*/
void SchubertContext::revertSize(const Ulong& n)
{
  Ulong m = size();

  while (m > n)
  {
    m -= d_history.top().size();
    d_history.pop();
  }
}


/*
  Resizes the various data structures to accomodate a context of size n.
  This means that the Lists d_length, d_hasse, d_descent and d_shift
  are resized to size n, and that memory is allocated for the new shift
  tables; we cannot do this for coatom lists, since they are variable in
  size.

  It is assumed that n is greater than the current size.

  Sets the error MEMORY_WARNING in case of overflow, if CATCH_MEMORY_OVERFLOW
  had been set.
*/
void SchubertContext::increase_size(const Ulong& n)
{
  Ulong prev_size = size();

  try {
    d_history.emplace(*this,n-size());
  }
  catch(...)
  {
    revertSize(prev_size);
  }
}

/******** input/output ****************************************************/

std::string& SchubertContext::append(std::string& str, coxtypes::CoxNbr x)
  const

{
  if (x == coxtypes::undef_coxnbr)
    str.append(undef_str);
  else
    coxtypes::append(str,x);

  return str;
}

std::string& SchubertContext::append
  (std::string& str, coxtypes::CoxNbr x,
   const interface::Interface& I) const

{
  if (x == coxtypes::undef_coxnbr)
    return str.append(undef_str);
  else {
    coxtypes::CoxWord g(0);
    normalForm(g,x,I.order());
    return I.append(str,g);
  }
}

void SchubertContext::print(FILE* file, coxtypes::CoxNbr x) const

{
  if (x == coxtypes::undef_coxnbr)
    fprintf(file,"%s",undef_str);
  else
    fprintf(file,"%lu",static_cast<Ulong>(x));

  return;
}

void SchubertContext::print(FILE* file, coxtypes::CoxNbr x,
				    const interface::Interface& I) const

{
  if (x == coxtypes::undef_coxnbr)
    fprintf(file,"%s",undef_str);
  else {
    coxtypes::CoxWord g(0);
    normalForm(g,x,I.order());
    I.print(file,g);
  }

  return;
}

/******** private member functions ******************************************

  The following functions are defined as private member functions. The
  main reason for this is that this gives access to the representation,
  so that we can inline the access (to the shift table for instance) that
  would otherwise go to a virtual function call.

   - fillCoatoms(first,s) : fills in the coatom lists of new elements of
     the extension by s;
   - fillDihedralShifts(x,s) : fills in the shifts in the case where
     x is dihedral;
   - fillShifts(first,s) : fills in the shift tables of new elements of
     the extension by s;
   - fillStar(first) : fills in the star tables of new elements;
   - fullExtension(q,s) : fills in the extension obtained by adding
     the elements xs, x in q;

*****************************************************************************/


/*
  This auxiliary fills the coatom lists of the new elements after extending for
  a shift by |s|. It is assumed that the context has been resized appropriately,
  thatd new elements start at |first|, and that lengths and shifts by |s| have
  already been filled in.
*/
void SchubertContext::fillCoatoms(const Ulong& first,  coxtypes::Generator s)
{

  for (coxtypes::CoxNbr x = first; x < d_size; ++x)
  {
    coxtypes::CoxNbr xs = shift(x,s);
    assert(xs<x);

    containers::sl_list<coxtypes::CoxNbr> c = {xs}; // start with this singleton

    for (coxtypes::CoxNbr z : d_hasse[xs])
    {
      coxtypes::CoxNbr zs = shift(z,s);
      assert(zs!=coxtypes::undef_coxnbr);
      if (zs > z) /* z moves up */
	c.push_back(zs);
    }

    assert(d_hasse.size()==x);
    d_hasse.push_back(c.to_vector()); // convert list tightly to vector
  }
}


/*
  Fill in the shifts for x in the dihedral case. It is assumed that the shift by
  s is already filled in, and that length(x) is > 1. We have denoted on the
  right the action of s, on the left the action on the side different from s.
  This works even if in fact the action of s is on the left.
*/

void SchubertContext::fillDihedralShifts
  (coxtypes::CoxNbr x, coxtypes::Generator s)
{
  coxtypes::CoxNbr xs = shift(x,s);

  /* find the other generator involved, on the same side as s */

  coxtypes::Generator s1, t, t1;
  graph::CoxEntry m;

  if (s < d_rank) { /* action is on the right */
    t = firstRDescent(xs);
    s1 = s + d_rank;
    t1 = t + d_rank;
    m = d_graph.M(s,t);
  }
  else { /* action is on the left */
    s1 = s - d_rank;
    t1 = firstLDescent(xs);
    t = t1 + d_rank;
    m = d_graph.M(s1,t1);
  }

  const CoxNbrList& c = d_hasse[x];
  coxtypes::CoxNbr z; /* the other coatom of x */

  if (c[0] == xs)
    z = c[1];
  else
    z = c[0];

  if (d_length[x] == m) { /* descents for s,t on both sides */
    d_descent[x] |=
      constants::eq_mask[t] | constants::eq_mask[s1] | constants::eq_mask[t1];
    d_downset[t].insert(x);
    d_downset[s1].insert(x);
    d_downset[t1].insert(x);
    d_shift.entry(x,t) = z;
    d_shift.entry(z,t) = x;
    if (m % 2) { /* xs = tx; xt = sx */
      d_shift.entry( x,s1) = z;
      d_shift.entry( z,s1) = x;
      d_shift.entry( x,t1) = xs;
      d_shift.entry(xs,t1) = x;
    }
    else { /* xs = sx; xt = tx */
      d_shift.entry( x,s1) = xs;
      d_shift.entry(xs,s1) = x;
      d_shift.entry( x,t1) = z;
      d_shift.entry( z,t1) = x;
    }
  }
  else { /* descent on one side only */
    if (d_length[x] % 2) { /* xs and sx */
      d_shift.entry(x,s1) = z;
      d_shift.entry(z,s1) = x;
      d_descent[x] |= constants::eq_mask[s1];
      d_downset[s1].insert(x);
    }
    else { /* xs and tx */
      d_shift.entry(x,t1) = z;
      d_shift.entry(z,t1) = x;
      d_descent[x] |= constants::eq_mask[t1];
      d_downset[t1].insert(x);
    }
  }
}


/*
  Fill in the shift tables of the new elements in p. It is assumed that first is
  the first new element, that the coatom tables, lengths and shifts by s have
  already been filled in. We use the algorithm deduced form Dyer's theorem,
  alluded to in the introduction.
*/

void SchubertContext::fillShifts(coxtypes::CoxNbr first, coxtypes::Generator s)
{
  coxtypes::CoxNbr x = first;

  /* check if something happens in length one; if there is a new element
   of length one, it is unique and equal to s */

  if (d_length[x] == 1) { /* x = s */
    coxtypes::Generator t;
    if (s < d_rank) /* s acts on the right */
      t = s + d_rank;
    else /* s acts on the left */
      t = s - d_rank;
    d_shift.entry(0,t) = x;
    d_shift.entry(x,t) = 0;
    d_descent[x] |= constants::eq_mask[t];
    d_downset[t].insert(x);
    ++x;
  }

  for (; x < d_size; ++x) {
    const CoxNbrList& c = d_hasse[x];

    if (c.size() == 2) { /* dihedral case */
      fillDihedralShifts(x,s);
      continue;
    }

    for (coxtypes::Generator t = 0; t < 2*d_rank; ++t) { /* examine shift by t */
      if (t == s)
	continue;
      bool firstplus = true;
      coxtypes::CoxNbr z = coxtypes::undef_coxnbr;
      for (Ulong j = 0; j < c.size(); ++j) {
	if (!(constants::eq_mask[t] & d_descent[c[j]])) { /* coatom has ascent */
	  if (firstplus) { /* it's the first time */
	    firstplus = false;
	    z = c[j]; // z is the coatom that goes up
	  }
	  else {
	    goto nextt;
	  }
	}
      }
      /* if we reach this point there was exactly one ascent */
      d_shift.entry(x,t) = z;
      d_shift.entry(z,t) = x;
      d_descent[x] |= constants::eq_mask[t];
      d_downset[t].insert(x);
    nextt:
      continue;
    }
  }

  return;
}


/*
  This function fills in the star operations for the new elements. Each star
  operation is a partially defined involution. The tables have already been
  initially set to |coxtypes::undef_coxnbr|; we fill in the operation in pairs,
  using the element that goes down.

  Recall that a star operation is associated to each edge {s,t} in the
  Coxeter graph such that m(s,t) < infty. The domain of the left star operation
  is the set of elements for which ldescent() intersects {s,t} in one element
  exactly (in other words, the elements that are neither minimal nor maximal
  in the left coset under the dihedral subgroup generated by s and t.)

  NOTE : a value coxtypes::undef_coxnbr means that either the element is not in
  the domain of the star operation, or that the star operation takes us out of
  context; hence an undefined value may become defined after extension; but this
  will always happen for elements paired up with a new element, so we only have
  to go through the new ones.
*/
void SchubertContext::fillStar(coxtypes::CoxNbr first)
{
  const containers::vector<GenSet>& ops = d_graph.finite_edges();

  for (coxtypes::CoxNbr x = first; x < d_size; ++x) {

    GenSet fx = rdescent(x);
    for (coxtypes::StarOp j = 0; j < nStarOps(); ++j) {

      // determine if x is in domain for right-star right descent set |fx|
      GenSet f = fx & ops[j];
      if ((f == 0) || (f == ops[j])) // must have singleton-intersect |ops[j]|
	continue;

      coxtypes::CoxNbr x_min = minimize(x,ops[j]);
      coxtypes::Length d = d_length[x] - d_length[x_min];
      coxtypes::Generator s =
	constants::firstBit(f); /* the _only_ bit in f, actually */
      coxtypes::Generator t = constants::firstBit(ops[j] & ~f);
      graph::CoxEntry m = d_graph.M(s,t);

      if (2*d < m) /* star is either coxtypes::undef_coxnbr or increasing */
	continue;

      /* if we get here we fill in a pair in d_star */

      if (2*d == m)
	d_star[x][j] = x;
      else {
	coxtypes::CoxNbr x1 = x;
	while ((d_length[x1] - d_length[x_min]) > (m - d)) {
	  GenSet f1 = rdescent(x1) & ops[j];
	  coxtypes::Generator s1 = constants::firstBit(f1);
	  x1 = shift(x1,s1);
	}
	d_star[x][j] = x1;
	d_star[x1][j] = x;
      }
    }

    fx = ldescent(x);
    for (coxtypes::StarOp j = 0; j < nStarOps(); ++j) {
      /* determine if x is in left domain */
      GenSet f = fx & ops[j];
      if ((f == 0) || (f == ops[j]))
	continue;
      Lflags lops = static_cast<Lflags>(ops[j]) << d_rank;
      coxtypes::CoxNbr x_min = minimize(x,lops);
      coxtypes::Length d = d_length[x] - d_length[x_min];
      coxtypes::Generator s =
	constants::firstBit(f); /* the _only_ bit in f, actually */
      coxtypes::Generator t = constants::firstBit(ops[j] & ~f);
      graph::CoxEntry m = d_graph.M(s,t);

      if (2*d < m) /* star is either coxtypes::undef_coxnbr or increasing */
	continue;

      /* if we get here we fill in a pair in d_star */

      if (2*d == m)
	d_star[x][j+nStarOps()] = x;
      else {
	coxtypes::CoxNbr x1 = x;
	while ((d_length[x1] - d_length[x_min]) > (m - d)) {
	  GenSet f1 = ldescent(x1) & ops[j];
	  coxtypes::Generator s1 = constants::firstBit(f1);
	  x1 = lshift(x1,s1);
	}
	d_star[x][j+nStarOps()] = x1;
	d_star[x1][j+nStarOps()] = x;
      }
    }
  }

  return;
}


/*
  Given a context p, a subset q of p holding [e,y], and a generator s s.t.
  y.s is not contained in p, this function extends p to hold y, and puts in
  q the interval [e,y.s] (here s can be either a right or a left shift.)

  Sets the following errors :

    - LENGHT_OVERFLOW if the new element has length greater than LENGTH_MAX
      (presumably this could have been checked before.)
    - COXNBR_OVERFLOW if the size of the extension would be greater than
      COXNBR_MAX;

  A more delicate problem is the handling of memory overflow. It has to be
  assumed that fullExtension labours under the constraint that
  CATCH_MEMORY_OVERFLOW is set; i.e., we don't want to exit brutally if
  we get a memory overflow, losing all previous computations.

  If an overflow error occurs, it is guaranteed that the context stays in
  its original form (except for sizes of varlists.)
*/
void SchubertContext::fullExtension(bits::SubSet& q, coxtypes::Generator s)
{
  /* check length overflow */

  coxtypes::CoxNbr y = q[q.size()-1]; /* largest element in q */

  if (d_length[y] == coxtypes::LENGTH_MAX) { /* overflow */
    ERRNO = LENGTH_OVERFLOW;
    return;
  }

  /* determine the size of the extension */

  coxtypes::CoxNbr c = 0; // counts number of new elements

  for (Ulong j = 0; j < q.size(); ++j) { /* run through q */
    if (shift(q[j],s) == coxtypes::undef_coxnbr)
      ++c;
  }

  /* check for size overflow */

  if (c > coxtypes::COXNBR_MAX - d_size) { // overflow predicited
    ERRNO = COXNBR_OVERFLOW;
    return;
  }

  /* resize context */

  coxtypes::CoxNbr prev_size = d_size;
  increase_size(d_size+c);
  q.setBitMapSize(d_size+c);

  if (ERRNO) /* memory overflow */
    goto revert;

  /* fill in lengths and shifts by s */

  { coxtypes::CoxNbr xs = prev_size; /* first new element */

    for (Ulong j = 0; j < q.size(); ++j)
    { coxtypes::CoxNbr x = q[j];
      if (shift(x,s) == coxtypes::undef_coxnbr)
      {
	d_shift.entry( x,s) = xs;
	d_shift.entry(xs,s) = x;
	d_length.push_back(d_length[x] + 1);
	d_parity[d_length[xs]%2].insert(xs);
	d_descent.push_back(constants::eq_mask[s]);
	d_downset[s].insert(xs);
	xs++;
      }
    }

    /* fill in the new elements */

    fillCoatoms(prev_size,s);
    fillShifts(prev_size,s);
    fillStar(prev_size);

    /* update q */

    extendSubSet(q,s);

    if (ERRNO)
      goto revert;
  }

  return;

 revert:
  revertSize(prev_size);
  return;
}

};

/****************************************************************************

        Chapter II -- The ContextExtension class.

  The ContextExtension class is provided to manage the resizings of the
  context, so that we can keep track of the pointers that are allocated
  to new memory.

  The following functions are provided :

    - ContextExtension(p,c) : builds the extension;
    - ~ContextExtension();

  NOTE : with reasonably managed arenas, this could probably be dropped;
  memory allocation in small chunks could be almost as fast and efficient
  as in big ones.

 ****************************************************************************/

namespace schubert {


/*
  This constructor also manages the resizing of the SchubertContext p from its
  current |size| to |size+c|. It is called through |d_history.emplace| from
  |SchubertContext::increase_size| (which on its turn is called by
  |SchubertContext::fullExtension|), so the object constructed ends up at the
  top of the |d_history| stack.
*/
SchubertContext::ContextExtension::ContextExtension
  (SchubertContext& p, const Ulong& c)
  : d_schubert(p)
  , d_size(c)
{
  if (c == 0)
    return;

  Ulong n = p.size()+c; // the new size for the context

  p.d_length.reserve(n);
  p.d_hasse.reserve(n);
  p.d_descent.reserve(n);
  p.d_shift.grow(n,coxtypes::undef_coxnbr); // extend with entirely undefined rows
  p.d_star.setSize(n);
  if (ERRNO)
    goto revert;

  /* make room for shift tables and star tables */

  d_star = new(memory::arena()) coxtypes::CoxNbr[2*p.nStarOps()*c];
  if (ERRNO)
    goto revert;
  memset(d_star,0xFF,2*p.nStarOps()*c*sizeof(coxtypes::CoxNbr));
  p.d_star[p.d_size] = d_star;
  for (Ulong j = p.d_size+1; j < n; ++j)
    p.d_star[j] = p.d_star[j-1] + 2*p.nStarOps();

  for (Ulong j = 0; j < 2ul*p.rank(); ++j)
    p.d_downset[j].set_capacity(n);

  p.d_parity[0].set_capacity(n);
  p.d_parity[1].set_capacity(n);

  if (ERRNO)
    goto revert;

  p.d_size = n;

  return;

 revert:
  p.d_shift.shrink(p.d_size);
  for (Ulong j = 0; j < 2*static_cast<Ulong>(p.rank()); ++j) {
    p.d_downset[j].set_capacity(p.d_size);
  }
  p.d_parity[0].set_capacity(p.d_size);
  p.d_parity[1].set_capacity(p.d_size);
  return;
}


/*
  Destruction of a context extension.

  NOTE : this is currently unfinished, and usable only for the destruction
  of a whole context. It should resize the lists downwards, and put
  coxtypes::undef_coxnbr values where appropriate (this can be determined by running
  through the deleted elements.) Also, it could take care of freeing
  the superfluous coatom lists.
*/
SchubertContext::ContextExtension::~ContextExtension()
{
  SchubertContext& p = d_schubert;
  Ulong prev_size = p.d_size-d_size;

  /* the pointers d_shift  and d_star were allocated previously */

  memory::arena().free(d_star,2*p.nStarOps()*d_size*sizeof(coxtypes::CoxNbr));

  p.d_size = prev_size;
} // |SchubertContext::ContextExtension::~ContextExtension|

};

/****************************************************************************

        Chapter II -- The ClosureIterator class.

  The closureIterator class is an iterator designed to loop over the context,
  providing at each step a SubSet holding the closure of element #y. This
  sort of loop will be needed in constructing tables of kl polynomials and
  mu-coefficients; it will be much more efficient to update the subset rather
  than reconstruct it from scratch at each stage (in fact, before this was
  introduced, closure extraction was the dominant function for mu-tables.)

  The following functions are defined :

    - constructors and destructors :

      - ClosureIterator(n) : initializes a ClosureIterator of size n;
      - ~ClosureIterator() : calls standard destructors on components;

    - iterator operators :

      - operator bool() : validity check (inlined);
      - operator()() : returns the current closure (inlined);
      - operator++() : increments the structure;

****************************************************************************/

namespace schubert {

ClosureIterator::ClosureIterator(const SchubertContext& p)
  : d_schubert(p)
  , state()
  , elements{0}
  , d_visited(p.size())
{
  d_visited.insert(0);
  state.push(node{p.closure(0), 1, coxtypes::CoxNbr(0), p.rascent(0)});
}


/*
  This function increments the iterator. In this case, this means updating
  |d_closure| to hold the closure of the next element.

  The important thing here is to choose carefully the order of traversal of
  p. We wish to do this in such a way that the subset can be managed as a
  stack. What we do is traverse in the lexicographical order of normal forms
  (so it will be Lex rather than ShortLex).

  The control structure that manages the traversal is a coxtypes::CoxWord,
  representing the normal form of the current element. The current element is
  initially zero. On update, the next element is the first extension of the
  current element within the context if there is such; otherwise the next
  extension of the parent of the current element, and so on. In order to
  facilitate things, we keep track of the number of elements visited.
*/
void ClosureIterator::operator++()
{
  const SchubertContext& p = d_schubert;
  elements.reserve(p.size()); // ensure, vitally, no reallocation will occur

  // find right ascents of |d_current|, possible extensions of the current word
  while(true) // loop exists by |return|, whether past the end or not
  { // |not state.empty()| is loop invariant
    auto& asc = state.top().asc;
    while (asc!=0) // there are candidates right ascents for |current()| left
    {
      coxtypes::Generator s = constants::firstBit(asc);
      asc &= asc-1; // clear the bit for |s|, whether or not it is productive
      coxtypes::CoxNbr x = p.shift(current(),s);
      if (x == coxtypes::undef_coxnbr) // don't leave the current context
	continue;
      if (d_visited.is_member(x)) // nor visit an element already visited
	continue;

      // now |x| is new and obtained by extending |current()| on the right by |s|
      d_visited.insert(x);
      const auto& prev_closure = closure(); // read this from old |state|
      state.push(node{prev_closure,0,x,p.rascent(x)}); // closure starts a copy
      auto& new_closure = state.top().closure;
      const auto end = elements.end(); // must fix current high water mark
      for (auto it = elements.begin(); it<end; ++it) // iterate over old part
      { auto ys = p.rshift(*it,s);
	if (not prev_closure.is_member(ys)) // avoid elements that are not new
	  { // all others will be added as new, and distinct
	    elements.push_back(ys); // no realloc, |it| and |end| still valid
	    new_closure.insert(ys);
	  }
      }
      state.top().closure_size = elements.size();
      return; // we're done incrementing, and still valid
    } // |for(s)|

    // if we get here, there are no extensions of the current word
    state.pop();
    if (state.empty())
      return; // our iterator is past the end
    elements.resize(state.top().closure_size); // drop no longer present elements
  } // |while(true)|
}

};

/****************************************************************************

        Chapter III -- Bruhat order and descent

  This section contains the definitions for the functions dealing with the
  Bruhat order and descent sets:

    - extractInvolutions(p,b) : extracts involutions;
    - shortLexOrder(p,x,y) : checks if x <= y in ShortLex order;

 ****************************************************************************/

namespace schubert {


/*
  This function extracts from b the involutions contained in it. First
  we check if L(x) = R(x); if yes, we check if x.x = e.
*/
void extractInvolutions(const SchubertContext& p, bits::BitMap& b)
{
  bits::BitMap::Iterator last = b.end();

  for (bits::BitMap::Iterator i = b.begin(); i != last; ++i) {
    coxtypes::CoxNbr x = *i;
    if (p.rdescent(x) != p.ldescent(x))
      goto not_involution;
    /* check if x.x = e */
    {
      coxtypes::CoxNbr xl = x;
      coxtypes::CoxNbr xr = x;
      while (xl) {
	coxtypes::Generator s = p.firstRDescent(xl);
	xl = p.rshift(xl,s);
	xr = p.lshift(xr,s);
	if (p.rdescent(xl) != p.ldescent(xr))
	  goto not_involution;
      }
    }
    /* if we get here, we have an involution */
    continue;
  not_involution:
    b.clearBit(x);
    continue;
  }
}


/*
  This function extracts from b the maximal elements w.r.t. f, by
  intersecting with the appropriate downsets.
*/
void select_maxima_for
  (const SchubertContext& p, bits::BitMap& b, const Lflags& f)
{
  Lflags f1 = f;

  while(f1) {
    coxtypes::Generator s = constants::firstBit(f1);
    b &= p.downset(s); // retain elements for which |s| is a right descent
    f1 &= f1-1;
  }
}

void select_maxima_for
  (const SchubertContext& p, bitmap::BitMap& b, Lflags f)
{
  while(f!=0)
  {
    coxtypes::Generator s = constants::firstBit(f); // either on right or left
    b &= p.down_set(s); // retain elements for which |s| is a descent
    f &= f-1;
  }
}


bool shortLexOrder(const SchubertContext& p, coxtypes::CoxNbr d_x,
		   coxtypes::CoxNbr d_y, const bits::Permutation& order)

/*
  This function checks if x <= y in the ShortLex order of the normal forms
  w.r.t. the given order (as usual, this means that we compare order[] to
  compare generators.) In other words, the result is true iff either length(x)
  < length(y), or lengths are equal and firstterm(x) < firstterm(y), or
  firstterms are equal (to s), and sx <= sy in ShortLex order.
*/

{
  if (d_x == d_y)
    return true;
  if (p.length(d_x) < p.length(d_y))
    return true;
  if (p.length(d_x) > p.length(d_y))
    return false;

  coxtypes::CoxNbr x = d_x;
  coxtypes::CoxNbr y = d_y;

  coxtypes::Generator s_x = p.firstLDescent(x,order);
  coxtypes::Generator s_y = p.firstLDescent(y,order);

  while (s_x == s_y) {
    x = p.lshift(x,s_x);
    y = p.lshift(y,s_y);
    s_x = p.firstLDescent(x,order);
    s_y = p.firstLDescent(y,order);
  }

  if (order[s_x] < order[s_y])
    return true;
  if (order[s_x] > order[s_y])
    return false;

  return false; // unreachable
}

};

/****************************************************************************

        Chapter IV -- Input/output

  This section defines the input/output functions declared in schubert.h :

   - print(file,p) : prints the context on a file;
   - printBitMap(file,b,p,I) : prints a BitMap;
   - printPartition(file,pi,p,I) : prints a partition;
   - printPartition(file,pi,b,p,I) : prints a partition restricted to a BitMap;

 ****************************************************************************/

namespace schubert {

void print(FILE* file, const SchubertContext& p)

/*
  This function prints out the contents of the Schubert context.
*/

{
  fprintf(file,"size : %lu  maxlength : %lu",static_cast<Ulong>(p.size()),
	  static_cast<Ulong>(p.maxlength()));
  fprintf(file,"\n\n");

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x) {
    fprintf(file,"%4lu : ",static_cast<Ulong>(x));

    for (coxtypes::Generator s = 0; s < p.rank(); ++s) {
      if (p.rshift(x,s) == coxtypes::undef_coxnbr)
	fprintf(file,"%4s","*");
      else
	fprintf(file,"%4lu",static_cast<Ulong>(p.rshift(x,s)));
    }
    fprintf(file,";");

    for (coxtypes::Generator s = 0; s < p.rank(); ++s) {
      if (p.lshift(x,s) == coxtypes::undef_coxnbr)
	fprintf(file,"%4s","*");
      else
	fprintf(file,"%4lu",static_cast<Ulong>(p.lshift(x,s)));
    }
    fprintf(file,";");

    fprintf(file,"  {");
    const CoxNbrList& c = p.hasse(x);
    for (Ulong j = 0; j < c.size(); ++j) {
      fprintf(file,"%lu",static_cast<Ulong>(c[j]));
      if (j+1 < c.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"}");

    fprintf(file,"  R:(");
    for (GenSet f = p.rdescent(x); f;) {
      fprintf(file,"%lu",static_cast<Ulong>(constants::firstBit(f)+1));
      f &= f-1;
      if (f) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,")");

    fprintf(file,"  L:(");
    for (GenSet f = p.ldescent(x); f;) {
      fprintf(file,"%lu",static_cast<Ulong>(constants::firstBit(f)+1));
      f &= f-1;
      if (f) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,")");

    fprintf(file,"\n");
  }

  fprintf(file,"\nStar operations :\n\n");

  for (coxtypes::CoxNbr x = 0; x < p.size(); ++x) {
    fprintf(file,"%4lu : ",static_cast<Ulong>(x));
    for (Ulong r = 0; r < 2*p.nStarOps(); ++r) {
      if (p.star(x,r) == coxtypes::undef_coxnbr)
	fprintf(file,"%5s","*");
      else
	fprintf(file,"%5lu",static_cast<Ulong>(p.star(x,r)));
    }
    fprintf(file,"\n");
  }

  fprintf(file,"\n");

  return;
}

void printBitMap(FILE* file, const bits::BitMap& b, const SchubertContext& p,
		    const interface::Interface& I)

/*
  This function prints the elements of the bitmap (assumed to hold a subset
  of the context) on the file.
*/

{
  bool first = true;

  fprintf(file,"{");

  for (bits::BitMap::Iterator i = b.begin(); i != b.end(); ++i) {
    if (first)
      first = false;
    else
      fprintf(file,",");
    coxtypes::CoxWord g(0);
    p.append(g,*i);
    I.print(file,g);
  }

  fprintf(file,"}");
}

void printPartition(FILE* file, const bits::Partition& pi, const SchubertContext& p,
		    const interface::Interface& I)

/*
  This function prints the partition pi, assumed to hold a partition of the
  context p.
*/

{
  Ulong count = 0;

  for (bits::PartitionIterator i(pi); i; ++i) {
    const Set& c = i();
    fprintf(file,"%lu(%lu):{",count,c.size());
    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxWord g(0);
      p.append(g,c[j]);
      I.print(file,g);
      if (j+1 < c.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"}\n");
    ++count;
  }

  return;
}

void printPartition(FILE* file, const bits::Partition& pi, const bits::BitMap& b,
		    const SchubertContext& p, const interface::Interface& I)

/*
  Prints the partition pi restricted to the subset flagged by b.
*/

{
  list::List<Ulong> q(b.begin(),b.end()); // replaces readBitMap
  bits::Partition pi_b(b.begin(),b.end(),pi);

  Ulong count = 0;

  for (bits::PartitionIterator i(pi_b); i; ++i) {
    const Set& c = i();
    fprintf(file,"%lu(%lu):{",count,c.size());
    for (Ulong j = 0; j < c.size(); ++j) {
      coxtypes::CoxWord g(0);
      p.append(g,q[c[j]]);
      I.print(file,g);
      if (j+1 < c.size()) /* there is more to come */
	fprintf(file,",");
    }
    fprintf(file,"}\n");
    ++count;
  }

  return;
}

};

/****************************************************************************

        Chapter V -- Utilities

  This section defines some utility functions :

   - betti(h,y,p) : puts in h the betti numbers of [e,y];
   - min(c,nfc) : extracts the minimal element from c;
   - extractMaximals(p,c) : extracts from c the list of its maximal elements;
   - minDescent(f,order) : finds the smallest element in f w.r.t. order;
   - readBitMap(c,b) : reads b into c;
   - setBitMap(b,c) : reads c into b;
   - sum(h) : returns the sum of the terms in h;

 ****************************************************************************/

namespace schubert {

void betti(Homology& h, coxtypes::CoxNbr y, const SchubertContext& p)

/*
  This function puts the ordinary betti numbers of the row in h, in a
  simple-minded approach. No overflow is possible here.

  It is assumed that row is a row in kllist.
*/

{
  bits::BitMap b(0);
  p.extractClosure(b,y);

  h.setSize(p.length(y)+1);
  h.setZero();

  bits::BitMap::Iterator b_end = b.end();

  for (bits::BitMap::Iterator x = b.begin(); x != b_end; ++x) {
    h[p.length(*x)]++;
  }

  return;
}

Ulong min(const Set& c, NFCompare& nfc)

/*
  This function extracts the minimal element form c. It is defined for
  Set instead of list::List<coxtypes::CoxNbr> so as to be able to apply it directly to
  partition classes.
*/

{
  if (c.size() == 0)
    return coxtypes::undef_coxnbr;

  Ulong m = c[0];

  for (Ulong j = 1; j < c.size(); ++j) {
    if (!nfc(m,c[j]))
      m = c[j];
  }

  return m;
}

void extractMaximals(const SchubertContext& p, list::List<coxtypes::CoxNbr>& c)

/*
  This function erases from c all elements that are not maximal elements for
  the Bruhat order among the entries in c.

  It is assumed that c is sorted in an ordering compatible with the Bruhat
  order; so if we start from the top, we will always encounter a maximal
  element before any lower one.
*/

{
  Ulong extr_count = 0;

  for (Ulong j = c.size(); j;) {
    --j;
    for (Ulong i = c.size()-extr_count; i < c.size(); ++i) {
      if (p.inOrder(c[j],c[i])) /* forget j */
	goto nextj;
    }
    extr_count++;
    c[c.size()-extr_count] = c[j];
  nextj:
    continue;
  }

  c.setData(c.ptr()+c.size()-extr_count,0,extr_count);
  c.setSize(extr_count);

  return;
}

void extractMaximals(const SchubertContext& p, list::List<coxtypes::CoxNbr>& c,
		     list::List<Ulong>& a)

/*
  Like the previous one, but puts the indices in c of the maximal elements
  in the list a.
*/

{
  list::List<coxtypes::CoxNbr> e(0);
  a.setSize(0);

  for (Ulong j = c.size(); j;) {
    --j;
    for (Ulong i = 0; i < e.size(); ++i) {
      if (p.inOrder(c[j],e[i])) /* forget j */
	goto nextj;
    }
    a.append(j);
    e.append(c[j]);
  nextj:
    continue;
  }

  a.reverse();

  return;
}


/*
  Returns the set bit position in f for which order is smallest. In practice,
  order is the external numbering of the generators; so this gives the
  internal number of the descent generator with smallest external number.

  NOTE : the return value is BITS(Ulong) if f is empty.
*/
Ulong minDescent(const GenSet& d_f, const bits::Permutation& order)
{
  GenSet f = d_f;
  Ulong m = constants::firstBit(f);
  f &= f-1;

  for (; f; f &= f-1) {
    Ulong m1 = constants::firstBit(f);
    if (order[m1] < order[m])
      m = m1;
  }

  return m;
}


/*
  This function reads in c from b (analogous to readBitMap in bits::SubSet).
*/
void readBitMap(list::List<coxtypes::CoxNbr>& c, const bits::BitMap& b)
{
  c.setSize(b.bitCount());

  bits::BitMap::Iterator i =  b.begin();

  for (Ulong j = 0; j < c.size(); ++j) {
    c[j] = *i;
    ++i;
  }
}

void read_bitmap(CoxNbrList& c, const bits::BitMap& b)
{
  c.reserve(b.bitCount());
  c.clear();

  bits::BitMap::Iterator i =  b.begin();

  for (Ulong j = 0; j < c.size(); ++j)
  {
    c.push_back(*i);
    ++i;
  }
}

};

namespace schubert {

void setBitMap(bits::BitMap& b, const list::List<coxtypes::CoxNbr>& c)

/*
  Reads c into b. It is assumed that b has already been set to the current
  context size.
*/

{
  b.reset();

  for (Ulong j = 0; j < c.size(); j++)
    b.setBit(c[j]);
}

Ulong sum(const Homology& h)

{
  Ulong a = 0;

  for (Ulong j = 0; j < h.size(); ++j) {
    a += h[j];
  }

  return a;
}

};

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
  , d_star(0,2*nStarOps()) // rows are filled on demand
  , d_downset(2*d_rank,bitmap::BitMap(1))
  , d_parity{bitmap::BitMap(1),bitmap::BitMap(1)}
{
  d_parity[0].insert(0);
}


/******** accessors ********************************************************/


/*
  Append to |g| the ShortLex normal form of x. The normal form is
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


/*
  This functions returns the number corresponding to g in the current
  context; returns coxtypes::undef_coxnbr if g is not in the context.
*/
coxtypes::CoxNbr SchubertContext::contextNumber(const coxtypes::CoxWord& g) const
{
  coxtypes::CoxNbr x = 0;

  for (Ulong j = 0; j < g.length(); ++j)
  {
    coxtypes::Generator s = g[j]-1;
    x = rshift(x,s);
    if (not in_context(x))
      return coxtypes::undef_coxnbr;
  }

  return x;
}


// the following function uses a double representation |q|,|elements| for speed
// iteration is over (initial) elements, (old) membership is tested though |q|
void SchubertContext::spread_subset
  (bitmap::BitMap& q, CoxNbrList& elements, coxtypes::Generator s) const
{
  elements.reserve(q.capacity()); // ensure, vitally, that no reallocation will occur
  const auto initial_end = elements.end(); // so that this iterator remains valid.
  for (auto it = elements.begin(); it!=initial_end; ++it) // don't visit new elements
  {
    const coxtypes::CoxNbr candidate = rshift(*it,s);
    if (not q.is_member(candidate)) // actually only excludes old elements
    {
      q.insert(candidate); // may or may not be new
      elements.push_back(candidate);
    }
  }
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
  result.insert(0);
  CoxNbrList elements {0}; // enumeration of |result| for faster iteration

  for (auto s : word(x))
    spread_subset(result,elements,s);

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
    if (not in_context(x1))
      return coxtypes::undef_coxnbr;
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
  (coxtypes::CoxWord& g, coxtypes::CoxNbr x, const bits::Permutation& order) const
{
  g.reset();

  while (x>0) {
    coxtypes::Generator s = firstLDescent(x,order);
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
coxtypes::CoxNbr SchubertContext::extend_context
  (const coxtypes::CoxWord& g)
{
  coxtypes::CoxNbr y = 0;
  bitmap::BitMap q(d_size);
  q.insert(0);
  CoxNbrList elements {0};

  Ulong j = 0;

  CATCH_MEMORY_OVERFLOW = true;

  for (; j < g.length(); ++j) {
    coxtypes::Generator s = g[j]-1;
    auto ys = rshift(y,s);
    if (not in_context(ys))
      break;
    spread_subset(q,elements,s);
    if (ERRNO)
      goto error_handling;
    y = ys;
  }

  for (; j < g.length(); ++j) {
    coxtypes::Generator s = g[j]-1;
    extend_context(q,elements,s); // adapts tables and enlarges |q|, |elements|
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
    if (q.is_member(xs))
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
void SchubertContext::revertSize(Ulong n)
{
  d_length.resize(n);
  d_hasse.resize(n);
  d_descent.resize(n);

  d_shift.shrink(n);
  for (Ulong j = 0; j < 2ul*rank(); ++j)
    d_downset[j].set_capacity(n);

  d_parity[0].set_capacity(n);
  d_parity[1].set_capacity(n);

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
void SchubertContext::increase_size(Ulong n)
{
  Ulong prev_size = size();

  try {
    d_length.reserve(n);
    d_hasse.reserve(n);
    d_descent.reserve(n);
    d_shift.grow(n,coxtypes::undef_coxnbr); // extend with entirely undefined row
    for (Ulong j = 0; j < 2ul*rank(); ++j)
      d_downset[j].set_capacity(n);

    d_parity[0].set_capacity(n);
    d_parity[1].set_capacity(n);

    d_size = n;
  }
  catch(...)
  {
    d_shift.shrink(prev_size);
    for (Ulong j = 0; j < 2ul*rank(); ++j)
      d_downset[j].set_capacity(prev_size);

    d_parity[0].set_capacity(prev_size);
    d_parity[1].set_capacity(prev_size);
    throw;
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

   - extend_context(q,elements,s) : fills in the extension obtained by adding
     the elements xs, x in q; the other methods do part of its work:
   - fill_Hasse(first,s) : fills in the coatom lists of new elements (from
     |first| to the end) added in the extension by |s|;
   - fill_shifts_and_descents(first,s) : fills in the shift tables of elements
     from |first| to end, which were added through a shift by |s|, and also
     complete the setting of their |d_descent| and |d_downset| entries
   - set_shifts_dihedral(x,s) : fills in the shifts in the case where
     x is dihedral;

*****************************************************************************/


/*
  When extending a |SchubertContext| by a shift |s| (which could be at the left
  or the right, althoug in practice it will be on the right), the shift by |s|
  for the new elements has been defined, but no other shifts. The next thing to
  define are their sets of coatoms, which can be computed in terms of shifts by
  |s| only; we do this for all new elements |x|, which start at |first|. Indeed
  the coatoms of the Bruhat interval $[e,x]$ are its elements of length one less
  than that of $x$, which are either |xs|, or (if new) obtained from an element
  of $[e,xs]$ by a shift $s$. In the latter case that element is for length
  reasons necessarily a coatom of |xs| for which |s| is an ascent, and all such
  coatoms give rise to a coatom of $[e,x]$, all distinct.
*/
void SchubertContext::fill_Hasse(Ulong first, coxtypes::Generator s)
{

  for (coxtypes::CoxNbr x = first; x < d_size; ++x)
  {
    coxtypes::CoxNbr xs = shift(x,s);
    assert(xs<x);

    containers::sl_list<coxtypes::CoxNbr> c = {xs}; // start with this singleton

    for (coxtypes::CoxNbr z : d_hasse[xs])
    {
      coxtypes::CoxNbr zs = shift(z,s);
      assert(in_context(zs));
      if (is_ascent(z,s))
	c.push_back(zs);
    }

    assert(d_hasse.size()==x);
    d_hasse.push_back(c.to_vector()); // convert list tightly to vector
  }
}


/*
  Fill in the shift tables of the new elements that were added to the context
  using a shift by |s|, where |first| is the first new element. Although our
  program uses right shifts (because |extend_context| is called with |s| equal
  to a letter of the Coxeter word that we are enlarging the context for), this
  function was written in a way that it works equally well for left shifts (this
  only occasionally make the code a bit larger or more intricate). When we come
  here, the tables have been filled for new elements only as far as coatoms,
  lengths and shifts by |s| are concerned; our task is to do the remainder. We
  use the algorithm deduced form Dyer's theorem referenced (in the form of
  corollary 3.7) in the paper by Fokko mentioned in the introduction.
*/
void SchubertContext::fill_shifts_and_descents
  (coxtypes::CoxNbr first, coxtypes::Generator s)
{
  const coxtypes::Rank two_r = 2*rank();

  coxtypes::CoxNbr x = first;

  // only the initial new element |first| might have length 1; check this case
  if (length(x) == 1)
  { assert(x==shift(0,s)); // we allow it to be on the left or right
    coxtypes::Generator t = (s + d_rank) % two_r; // |s| from opposite side
    d_shift.entry(0,t) = x;
    d_shift.entry(x,t) = 0;
    d_descent[x] |= constants::eq_mask[t];
    d_downset[t].insert(x);
    ++x; // don't treat this case in next loop
  }

  for (; x < d_size; ++x)
  {
    const CoxNbrList& coatoms = d_hasse[x]; // this was already computed

    if (coatoms.size() == 2) // treat the "dihedral case" separately
    {
      const coxtypes::CoxNbr xs = shift(x,s);

      coxtypes::Generator t = // the other dihedral group generator, same side
	s < d_rank ? firstRDescent(xs) : firstLDescent(xs) + d_rank;

      coxtypes::CoxNbr z = // the coatom of |x| that is not |xs|
	 coatoms[0]+coatoms[1]-xs;
      coxtypes::Generator s_op = // on opposite side such that |z=shift(x,s_op)|
	((length(x)%2==0 ? t : s) + rank()) % two_r;

      d_shift.entry(x,s_op) = z;
      d_shift.entry(z,s_op) = x;
      d_descent[x] |= constants::eq_mask[s_op];
      d_downset[s_op].insert(x);

      if (length(x)==d_graph.M(s%rank(),t%rank())) // |x| tops dihedral group?
      { // if so, |z| is also a descent of |x| (on the same side as |s|) by |t|
	d_shift.entry(x,t) = z;
	d_shift.entry(z,t) = x;
	d_descent[x] |= constants::eq_mask[t];
	d_downset[t].insert(x);
	// and |xs| is also opposite side descent of |x|, by complement of |s_op|
	coxtypes::Generator t_op = // same side as |s_op| and different from it
	  (two_r + s + t - s_op)%two_r; // uses of |two_r| to avoid underflow
	d_shift.entry(x ,t_op) = xs;
	d_shift.entry(xs,t_op) = x ;
	d_descent[x] |= constants::eq_mask[t_op];
	d_downset[t_op].insert(x);
      }
    }

    else
      for (coxtypes::Generator t = 0; t < two_r; ++t)
	if (t != s)
	{ // examine shift by t
	  coxtypes::CoxNbr z = coxtypes::undef_coxnbr;
	  for (auto coatom : coatoms)
	    // |coatom<x|, so we inductively know |d_descent[coatom]| is complete
	    if (is_ascent(coatom,t)) // so this condition can be tested
	    { // now |coatom| has an ascent for |t|
	      if (z == coxtypes::undef_coxnbr) // this is the first occurence
		z = coatom; // so record that coatom in |z|
	      else // more than one coatom has an ascent for |t|
		goto next_t; // then this |t| needs no action here, forget |z|
	    } // |for(coatom)|

	  /* if we arrive here, |z| is the unique coatom with ascent for |t|; by
	    the corollary, this implies |shift(x,t)=z|: it ascends to |x| (if
	    |t| is a descent for |y| then the element it descends to is the
	    unique coatom with an ascent for |t|; the corollary says that, after
	    excluding the dihedral case, this does not happen for other |t|).
	  */
	  d_shift.entry(x,t) = z;
	  d_shift.entry(z,t) = x;
	  d_descent[x] |= constants::eq_mask[t];
	  d_downset[t].insert(x);
	next_t:
	  continue; // (an empty statement would do the same)
	} // |for(t)|
  } // |for(x)|

} // |fill_shifts_and_descents|

template<bool left> coxtypes::CoxNbr SchubertContext::falling_star
  (coxtypes::CoxNbr x, GenSet st) const
{
  const GenSet f = left ? ldescent(x) : rdescent(x);
  const GenSet f0 = f  & st;
  const GenSet f1 = f0 ^ st; // its complement in $\{s,t\}$
  if (f0 == 0 or f1 == 0) // must have singleton-intersection with $\{s,t\}$
    return coxtypes::undef_coxnbr;

  const Lflags sided_st = left ? static_cast<Lflags>(st) << rank() : st;
  const coxtypes::CoxNbr x_min = minimize(x,sided_st);
  const coxtypes::Length dl = d_length[x] - d_length[x_min];
  coxtypes::Generator s = constants::firstBit(f0);
  coxtypes::Generator t = constants::firstBit(f1); // the _only_ bits
  const auto m = d_graph.M(s,t);
  if (left)
    s += rank(), t+=rank(); // so shifts below will be left shifts

  if (2*dl<m) // whether star operation will increase length
    return coxtypes::undef_coxnbr;

  unsigned count = 2*dl-m; // number of steps to take

  while (count-->0)
  {
    assert(isDescent(x,s));
    x = shift(x,s);
    assert(x!=coxtypes::undef_coxnbr);
    std::swap(s,t);
  }
  return x;
} // |SchubertContext::falling_star|

/*
  Extend the |d_star| table to have |size()| rows. This function is called by
  the |star| method when the table has no row yet for the requested group
  element, as will certainly be the case the first time that method is called.
  We then fill the table for all currently known elements, so that we avoid
  repetitive resizing of the |d_star| matrix.

  Recall that a star operation is associated to each pair $\{s,t\}$ of
  generators with $2 < m(s,t) < \infty$. The domain of the left star operation
  for this edge is the set of elements whose |ldescent| intersects $\{s,t\}$ in
  in a singleton (in other words, the elements that are neither minimal nor
  maximal in the left coset under the dihedral subgroup for $s$ and $t$.)

  NOTE : a value coxtypes::undef_coxnbr means that either the element is not in
  the domain of the star operation, or that the star operation takes us out of
  context; hence an undefined value may become defined after extension; but this
  will always happen for elements paired up with a new element, so we only have
  to go through the new ones.
*/
void SchubertContext::fill_star_table()
{
  const coxtypes::CoxNbr old_size = d_star.nr_rows();
  d_star.grow(size(),coxtypes::undef_coxnbr);

  const containers::vector<GenSet>& ops = d_graph.finite_edges();
  assert(d_star.nr_cols()==2*ops.size());

  for (coxtypes::CoxNbr x = old_size; x < d_size; ++x)
    for (coxtypes::StarOp j = 0; j < ops.size(); ++j)
    {
      auto y = falling_star<false>(x,ops[j]);
      if (y != coxtypes::undef_coxnbr)
      { assert(y<=x);
	right_star_row(x)[j] = y;
	right_star_row(y)[j] = x;
      }
      y = falling_star<true>(x,ops[j]);
      if (y != coxtypes::undef_coxnbr)
      { assert(y<=x);
	left_star_row(x)[j] = y;
	left_star_row(y)[j] = x;
      }
    }

} // |SchubertContext::fill_star_table|


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
  assumed that |extend_context| labours under the constraint that
  CATCH_MEMORY_OVERFLOW is set; i.e., we don't want to exit brutally if
  we get a memory overflow, losing all previous computations.

  If an overflow error occurs, it is guaranteed that the context stays in
  its original form (except for sizes of varlists.)
*/
void SchubertContext::extend_context
  (bitmap::BitMap& q, CoxNbrList& elements, coxtypes::Generator s)
{
  /* check length overflow */
  coxtypes::CoxNbr y = elements.back(); /* most recent element in q */

  if (d_length[y] == coxtypes::LENGTH_MAX) { /* overflow */
    ERRNO = LENGTH_OVERFLOW;
    return;
  }

  /* determine the size of the extension */
  coxtypes::CoxNbr c = // number of currently undefined shifts by |s|
    std::count_if
      (elements.begin(),elements.end(),
       [s,this](coxtypes::CoxNbr x) { return not in_context(shift(x,s)); }
      );

  /* check for size overflow */
  if (c > coxtypes::COXNBR_MAX - d_size) { // overflow predicited
    ERRNO = COXNBR_OVERFLOW;
    return;
  }

  /* resize context */

  coxtypes::CoxNbr prev_size = d_size;
  increase_size(d_size+c); // enlarge capacity of internal tables; do not fill
  q.set_capacity(d_size+c);

  // fill in lengths, and shifts by |s|
  for (coxtypes::CoxNbr x : elements)
    if (shift(x,s) == coxtypes::undef_coxnbr)
    { // create new entry and fill in only very basic attributes
      const coxtypes::CoxNbr xs = d_length.size();
      d_length.push_back(d_length[x] + 1);
      d_descent.push_back(constants::eq_mask[s]);
      d_shift.entry( x,s) = xs;
      d_shift.entry(xs,s) = x;
      d_parity[d_length[xs]%2].insert(xs);
      d_downset[s].insert(xs);
    }

  // complete information for the new elements
  fill_Hasse(prev_size,s); // now coatoms for closures are defined
  fill_shifts_and_descents(prev_size,s);

  spread_subset(q,elements,s);  // update |q| and |elements| for next round

  if (ERRNO)
    goto revert;
  return;

 revert:
  revertSize(prev_size);
  return;
}

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
  , sp(0) // there will be one state entry to start
  , elements{0}
  , visited(p.size())
{
  visited.insert(0);
  state.reserve(p.max_length()+1);
  for (unsigned i=0; i<=p.max_length(); ++i)
    state.push_back(p.size()); // create all bitmaps to avoid later (de)allocations
  state[0].closure.insert(0);
  state[0].closure_size=1;
  state[0].current = coxtypes::CoxNbr(0);
  state[0].asc = p.rascent(0);

  elements.reserve(p.size()); // ensure, vitally, no reallocation will occur
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

  // find right ascents of |d_current|, possible extensions of the current word
  while(true) // loop exits by |return|, whether past the end or not
  { // |not state.empty()| is loop invariant
    auto& asc = state[sp].asc;
    while (asc!=0) // there are still candidate right ascents for |current()|
    {
      coxtypes::Generator s = constants::first_bit(asc);
      asc &= asc-1; // clear the bit for |s|, whether or not it is productive
      coxtypes::CoxNbr x = p.shift(current(),s);
      if (not p.in_context(x) or visited.is_member(x))
	continue;

      // now |x| is new and obtained by extending |current()| on the right by |s|
      visited.insert(x);
      const auto& prev_closure = closure(); // read this from old |state|
      assert(sp+1<state.size()); // we shoudld not run out of |state| s[ace
      state[++sp].current = x; // push a new node on stack (filled below)
      auto& new_closure = state[sp].closure = prev_closure; // intially copy memory
      state[sp].asc = p.rascent(x);
      const auto end = elements.end(); // must fix current high water mark
      for (auto it = elements.begin(); it<end; ++it) // iterate over old part
      { auto ys = p.rshift(*it,s);
	if (not prev_closure.is_member(ys)) // avoid elements that are not new
	{ // all others will be added as new, and distinct
	  elements.push_back(ys); // no realloc, |it| and |end| still valid
	  new_closure.insert(ys);
	}
      }
      state[sp].closure_size = elements.size();
      return; // we're done incrementing, and still valid
    } // |for(s)|

    // if we get here, there are no extensions of the current word
    if (sp==0)
    { state.clear(); // this signals iteration has terminated
      return; // our iterator is past the end
    }
    --sp; // otherwise, pop stack and drop the no longer present |elements|
    elements.erase(elements.begin()+state[sp].closure_size,elements.end());
  } // |while(true)|
}

};

/****************************************************************************

        Chapter III -- Bruhat order and descent

  This section contains the definitions for the functions dealing with the
  Bruhat order and descent sets:

    - extractInvolutions(p,b) : extracts involutions;
    - shortlex_leq(p,order,x,y) : whether $x\leq y$ in (|order| collated) ShortLex;

 ****************************************************************************/

namespace schubert {


/*
  This function extracts from b the involutions contained in it. First
  we check if L(x) = R(x); if yes, we check if x.x = e.
*/
bool is_involution(const SchubertContext& p,coxtypes::CoxNbr x)
{
  if (p.rdescent(x) != p.ldescent(x))
    return false;

  coxtypes::CoxNbr y = x;
  while (x!=0)
  { // we know |p.rdescent(x) == p.ldescent(y)|, so descent through an |s| in it
    coxtypes::Generator s = p.firstRDescent(x);
    x = p.rshift(x,s);
    y = p.lshift(y,s);
    if (p.rdescent(x) != p.ldescent(y) or p.ldescent(x)!=p.rdescent(y))
      return false;
  }
  return true;
}


void select_maxima_for
  (const SchubertContext& p, Lflags f, bitmap::BitMap& b)
{
  while(f!=0)
  {
    coxtypes::Generator s = constants::first_bit(f); // either on right or left
    b &= p.down_set(s); // retain elements for which |s| is a descent
    f &= f-1;
  }
}



/*
  Test whether $x \leq y$ in the ShortLex order of the normal forms  w.r.t. the given order (as usual, this means that we compare order[] to
  compare generators.) In other words, the result is true iff either length(x)
  < length(y), or lengths are equal and firstterm(x) < firstterm(y), or
  firstterms are equal (to s), and sx <= sy in ShortLex order.
*/
bool shortlex_leq(const SchubertContext& p, const bits::Permutation& order,
		  coxtypes::CoxNbr x, coxtypes::CoxNbr y)
{
  if (x == y) // since this test is fast, get the equality case out of the way
    return true;
  if (p.length(x) != p.length(y))
    return p.length(x) < p.length(y);

  coxtypes::Generator s_x = p.firstLDescent(x,order); // runs over letters of |x|
  coxtypes::Generator s_y = p.firstLDescent(y,order); // runs over letters of |y|

  while (s_x == s_y)
  {
    assert(x>0 and y>0); // since equality and unequal length cases are handled
    s_x = p.firstLDescent(x = p.lshift(x,s_x),order);
    s_y = p.firstLDescent(y = p.lshift(y,s_y),order);
  }

  // now we are at the fisrt difference of letters
  return order[s_x] < order[s_y];
}

};

/****************************************************************************

        Chapter IV -- Input/output

  This section defines the input/output functions declared in schubert.h :

   - print(file,p) : prints the context on a file;
   - printPartition(file,pi,p,I) : prints a partition;

 ****************************************************************************/

namespace schubert {


/*
  This function prints out the contents of the Schubert context.
*/
void print(FILE* file, SchubertContext& p)
{
  fprintf(file,"size : %lu  maxlength : %lu",static_cast<Ulong>(p.size()),
	  static_cast<Ulong>(p.max_length()));
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


/*
  This function prints the partition pi, assumed to hold a partition of the
  context p.
*/

void printPartition
  (FILE* file, const bits::Partition& pi, const SchubertContext& p,
   const interface::Interface& I)
{
  Ulong count = 0;

  for (bits::Partition::iterator pit=pi.begin(); pit; ++pit)
  {
    auto range = *pit;
    fprintf(file,"%lu(%lu):{",count,range.size());
    for (const auto& elt : range) {
      coxtypes::CoxWord g(0);
      p.append(g,elt);
      I.print(file,g);
      if (&elt+1 != &*range.end()) /* there is more to come */
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

   - h = betti(y,p) : puts in h the betti numbers of [e,y];
   - min(c,nfc) : extracts the minimal element from c;
   - extractMaximals(p,c) : extracts from c the list of its maximal elements;
   - first_flagged(f,order) : finds the smallest element in f w.r.t. order;
   - readBitMap(c,b) : reads b into c;
   - setBitMap(b,c) : reads c into b;
   - sum(h) : returns the sum of the terms in h;

 ****************************************************************************/

namespace schubert {


/*
  Return the list of ordinary Betti numbers for |y|, by a simple-minded counting.
  No overflow is possible here.
*/

Homology betti(coxtypes::CoxNbr y, const SchubertContext& p)
{
  auto cl = p.closure(y);

  Homology result(p.length(y)+1,0);

  for (auto it = cl.begin(); it(); ++it)
    ++result[p.length(*it)]; // make length generating polynomial on closere set

  return result;
}

containers::sl_list<Ulong> indices_of_maxima
(const SchubertContext& p, containers::vector<coxtypes::CoxNbr>& c)
{
  bitmap::BitMap seen(p.size());
  containers::sl_list<Ulong> result(p. size());

  for (auto it=c.rbegin(); it!=c.rend(); ++it)
  { coxtypes::CoxNbr x = *it;
    if (not seen.is_member(x))
    { result.push_front(&*it-&c[0]);
      seen |= p.closure(x);
    }
  }

  return result;
}

/*
  Return the set bit position in f for which order is smallest. In practice,
  order is the external numbering of the generators; so this gives the
  internal number of the descent generator with smallest external number.

  NOTE : the return value is BITS(Ulong) if f is empty.
*/
coxtypes::Generator first_flagged(GenSet f, const bits::Permutation& order)
{
  Permutation inv = order.inverse();
  for (coxtypes::Generator s : inv)
    if ((f&constants::eq_mask[s])!=0)
      return s;
  return BITS(Ulong); // better assert(false) or throw an error
}


Ulong sum(const Homology& h)
{
  Ulong s = 0;

  for (auto coef : h)
    s += coef;

  return s;
}

}; // |namespace schubert|

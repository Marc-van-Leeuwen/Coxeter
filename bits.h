/*
  This is bits.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef BITS_H  /* guarantee single inclusion */
#define BITS_H

#include "globals.h"

#include <limits.h>
#include <iterator>
#include <algorithm>
#include <new>
#include <string>
#include <map>

#include "containers.h"
#include "bitmap.h"
#include "io.h"
#include "constants.h"

/******** type declarations *************************************************/

namespace bits {
  class Partition;
  class Permutation;
  class SubSet;
  using SetElt = Ulong;

/******** function declarations *********************************************/

  unsigned bitCount(const Lflags& f);
  bool isRefinement(const Partition& pi1, const Partition& pi2);
  template <class T> void rightRangePermute
    (list::List<T>& r, const Permutation& a);
  template <class T> void sortI(const list::List<T>& r, Permutation& a);
  template <class T> void sortI(const containers::vector<T>& r, Permutation& a);
  template <class T, class C> void sortI(const list::List<T>& r, C& inOrder,
					 Permutation& a);
  template <class T, class F> void sortI_f(const list::List<T>& r, F& f,
					   Permutation& a);
  template <class T> void sortI(const containers::vector<T>& r, Permutation& a);
  template <class T, class C> Permutation inverse_standardization
    (const containers::vector<T>& r, C& less);
  template <class T, class F> void sortI_f
    (const containers::vector<T>& r, F& f, Permutation& a);

/******** type definitions **************************************************/


class Permutation
  : public containers::vector<Ulong>
{
 public:
/* constructors and destructors */
  Permutation();
  Permutation(containers::vector<Ulong> v)
    : containers::vector<Ulong>(std::move(v)) {}
  Permutation(Ulong n); // identity of size |n|
  Permutation(Ulong n, Ulong val) : containers::vector<Ulong>(n,val) {}
/* manipulators */
  Permutation& identity(Ulong n); // assign identity to |*this| (returned)
  Permutation inverse() const;
  Permutation& compose(const Permutation& a);
  Permutation& rightCompose(const Permutation& a);
}; // |class Permutation|


class SubSet
{
 private:
  bitmap::BitMap d_bitmap;
  containers::vector<Ulong> row;
 public:
/* constructors and destructors */
  SubSet() {};
  SubSet(Ulong n) : d_bitmap(n), row() {}
  SubSet(containers::vector<Ulong>&& row,Ulong n)
    : d_bitmap(n), row(std::move(row))
  { for (Ulong entry : row) d_bitmap.insert(entry); }
  SubSet(const SubSet& q):d_bitmap(q.d_bitmap), row(q.row) {}
/* accessors */
  const containers::vector<Ulong>& elements() const { return row; }
  const Ulong operator[] (Ulong j) const { return row[j]; }
  const bitmap::BitMap& bitMap() const { return d_bitmap; }
  Ulong find(const SetElt& x) const
  { return std::lower_bound(row.begin(),row.end(),x)-row.begin(); }
  bool is_member(Ulong n) const { return d_bitmap.is_member(n); }
  Ulong size() const { return row.size(); }
  containers::vector<Ulong>::const_iterator begin() const { return row.begin(); }
  containers::vector<Ulong>::const_iterator end() const {return row.end(); }
/* modifiers */
  Ulong& operator[] (Ulong j) { return row[j]; }
  void add(Ulong n); // add |n| to bitmap and list, unless present
  SubSet& assign(const SubSet& q) { new (this) SubSet(q); return *this; }
  bitmap::BitMap& bitMap() { return d_bitmap; }
  void readBitMap(); // set |row| from |d_bitmap| contents
  void reset();
  void setBitMapSize(Ulong n) { d_bitmap.set_capacity(n); }
  void setListSize(Ulong n) { row.resize(n); }
  void sortList() { return std::sort(row.begin(),row.end()); }
}; // |SubSet|


class Partition
{
private:
  containers::vector<Ulong> classifier;
  Permutation inv_stnd; // inverse standardization
  containers::vector<size_t> cum_class_sizes; // cumulative class sizes
 public:
/* class definitions */
  typedef Ulong result_type;
/* constructors and destructors */
  Partition() : classifier(0), inv_stnd(0), cum_class_sizes() {}
  Partition(containers::vector<Ulong>&& v, Ulong cc); // trust the caller
  // Partition(const Partition& a, const BitMap& b); // partition of subset
  template <class F> Partition(Ulong n, const F& property);
  template <class I, class F> Partition
    (const I& first, const I& last, const F& f);
  template <class I> Partition // from segments of input iteration
    (I first, Ulong count, const containers::vector<Ulong>& class_sizes);
/* accessors */
  // the next method allows treating the Partition as a (classifying) function
  Ulong operator() (Ulong j) const { return classifier[j]; }
  Ulong class_count() const {return cum_class_sizes.size();}
  Ulong size() const { return classifier.size(); }

  Permutation standardization() const; // standardization of |classifier|
  const Permutation& inverse_standardization() const { return inv_stnd; }
  SubSet class_nr(Ulong i) const; // equivalence class with label |n|
  SubSet class_of(Ulong x) const { return class_nr(classifier[x]); }
  Ulong class_size(Ulong i) const
  { return i==0 ? cum_class_sizes[0] : cum_class_sizes[i]-cum_class_sizes[i-1]; }
  containers::vector<Ulong> class_sizes() const;
/* modifiers */
  bool refine(const containers::vector<Partition>& L);
  Partition& normalize() { return *this = Partition(size(),*this); }
  void permute_base(const Permutation& a);
  void permute_range(const Permutation& a);

  using iter = containers::vector<Ulong>::iterator;
  using citer = containers::vector<Ulong>::const_iterator;
  struct range
  { citer a,b;
    citer begin() const { return a; }
    citer end() const { return b; }
    size_t size() const { return b-a; }
    Ulong operator[] (Ulong i) const { return *(a+i); }
  }; // |struct range|

  class iterator
  {
    const Partition& parent;
    Ulong class_nr; // a class number
  public:
/* constructors and destructors */
    iterator(const Partition& pi,Ulong i) : parent(pi), class_nr(i) {}
/* iterator operations */
    operator bool() const { return class_nr<parent.class_count(); }
    bool operator!= (const iterator& e) const { return class_nr!=e.class_nr; }
    bool operator== (const iterator& e) const { return class_nr==e.class_nr; }
    void operator++() { ++class_nr; }
    range operator*() const
    { return range{parent.class_bound(class_nr),parent.class_bound(class_nr+1)}; }
  }; // |iterator|

  iterator begin() const { return iterator(*this,0); }
  iterator end() const { return iterator(*this,class_count()); }

private:
  citer class_bound(Ulong i) const // where |0<=i<=class_count()|
  { return inv_stnd.cbegin()+(i==0 ? 0 : cum_class_sizes[i-1]); }
  iter class_bound(Ulong i) // where |0<=i<=class_count()|
  { return inv_stnd.begin()+(i==0 ? 0 : cum_class_sizes[i-1]); }
  void set_inverse_standardization(Ulong class_count);
}; // |Partition|


/**** Inline implementations **********************************************/

/******** template definitions ***********************************************/


/*
  This constructor defines the partition of [0,r.size()[ defined by f; f
  is supposed to be a function taking arguments of type T. It is also
  assumed that operator<= is defined for the value type of f (so that
  the function insert may be applied.)
*/
template <class F> Partition::Partition
  (Ulong n, const F& property) : classifier(), inv_stnd(), cum_class_sizes()
{
  std::map<decltype(property(0)),Ulong> renumber;

  Ulong class_count = 0;
  for (Ulong i=0; i<n; ++i)
    if (renumber.insert(std::make_pair(property(i),class_count)).second)
      ++class_count; // |property(i)| was inserted (new key); consolidate value

  classifier.reserve(n);
  for (Ulong i=0; i<n; ++i)
    classifier.push_back(renumber[property(i)]);
  set_inverse_standardization(class_count);
}


/*
  A rather general partition constructor. It is assumed that I is an Input
  Iterator (in an informal sense). It is assumed that f is a functor taking
  one argument of type the value type of I, and that operator< is defined
  for the value type of f. A partition is constructed on the range [0,N[,
  where N is the number of iterations from first to last; the class numbers
  are attributed in the order of the values of f on the range.
*/
template <class I, class F> Partition::Partition
  (const I& first, const I& last, const F& f)
    : classifier(), inv_stnd(), cum_class_sizes()
{
  std::map<decltype(f(*first)),Ulong> renumber;

  Ulong size=0;
  for (auto it=first; it!=last; ++it,++size)
    renumber.insert(std::make_pair(f(*it),0)); // form |renumber| as a mere set

  Ulong class_count=0;
  for (auto& p : renumber)
    p.second = class_count++; // fill in class numbering, increasingly

  classifier.reserve(size);
  for (auto it=first; it!=last; ++it)
    classifier.push_back(renumber[f(*it)]);
  set_inverse_standardization(class_count);
}

template <class I> Partition::Partition
  (I first, Ulong count, const containers::vector<Ulong>& class_sizes)
    : classifier(count,-1)
    , inv_stnd(count,-1)
    , cum_class_sizes(class_sizes)
{
  // cumulate left-to-right, giving starting index of each class (for now)
  Ulong sum=0;
  for (Ulong& count : cum_class_sizes)
  { std::swap(count,sum); // current sum replaces count
    sum += count; // which get added to sum
  }

  classifier.reserve(count);
  for (Ulong class_nr = 0; class_nr<class_sizes.size(); ++class_nr)
  {
    auto& index = cum_class_sizes[class_nr]; // index into |inv_stnd.| to use
    for (Ulong i = 0; i<class_sizes[class_nr]; ++i)
    {
      classifier[*first++]=class_nr;
      inv_stnd[index++] = class_nr;
    }
    assert(class_nr+1==class_sizes.size() or index==cum_class_sizes[class_nr+1]);
  }

  assert(classifier.size()==count); // check that |count| is sum of class sizes
}

/*
  Apply the permutation |a| to the range of the list, on the right (this
  amounts to applying the inverse permutation in the usual sense). In
  other words, we have new[j] = old[a[j]].

  We cannot write this directly however, or we would overwrite. So we
  do something similar as with ordinary range permutations.
*/
template <class T> void rightRangePermute
  (list::List<T>& r, const Permutation& a)
{
  bitmap::BitMap b(r.size());

  for (Ulong j = 0; j < a.size(); ++j) {

    if (b.is_member(j))
      continue;

    if (a[j] == j) {
      b.insert(j);
      continue;
    }

    Ulong k = j;
    b.insert(j);

    for (Ulong i = a[j]; i != j; i = a[i]) {
      T buf = r[k];
      r[k] = r[i];
      r[i] = buf;
      k = i;
      b.insert(i);
    }
  }
} // |rightRangePermute|

// a more straightforward version of the above, using move semantics
template <class T>
  void right_permute (containers::vector<T>& r, const Permutation& a)
{
  assert(r.size()==a.size());
  containers::vector<T> result; result.reserve(r.size());
  for (Ulong i=0; i<a.size(); ++i)
  {
    assert(a[i]<r.size());
    result.push_back(std::move(r[a[i]]));
  }
  r = std::move(result);
} // |right_permute|


/*
  General sort function for lists. It is assumed that operator < is defined
  for T.

  Doesn't actually modify r; it only writes down in a the permutation
  s.t. new[j] = old[a[j]].
*/
template <class T> void sortI(const list::List<T>& r, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */

  Ulong h = 1;

  for (; h < r.size()/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < r.size(); ++j) {
      Ulong buf = a[j];
      Ulong i = j;
      for (; (i >= h) && (r[a[i-h]] > r[buf]); i -= h)
	a[i] = a[i-h];
      a[i] = buf;
    }
  }
} // |sortI|

template <class T> void sortI(const containers::vector<T>& r, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */
  Ulong h = 1;
  while (h < r.size()/3)
    h = 3*h+1;

  for (; h > 0; h /= 3)
  { // sort classes modulo |h| of indices in |a|
    for (Ulong j = h; j < a.size(); ++j)
    {
      const Ulong key = a[j]; // an index into |r|
      const T& val = r[key];
      // insertion sort |key| into subarray of |a| at |h|-class of |j| upto |j|
      Ulong i = j; // the index into |a| where a value was last moved from
      for ( ; i>=h; i-=h)
      {
	if (not val < r[a[i-h]])
	  break; // don't use |i-h|, the place is |i|
	a[i] = a[i-h]; // move index |h| places up, leaving |a[i-h]| moved from
      }
      // now |i| is the index into |a| to move |key| to
      a[i] = key;
    }
  }
} // |sortI(r,a)|


/*
  General sort function taking a comparison functor as a parameter.
  It is assumed that inOrder takes two arguments of type T, and returns
  a boolean value, so that the corresponding relation is a total preorder
  relation.

  Doesn't actually modify r; it only writes down in a the permutation
  s.t. new[j] = old[a[j]].
*/
template <class T, class C> void sortI
  (const list::List<T>& r, C& inOrder, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */

  Ulong h = 1;

  for (; h < r.size()/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < r.size(); ++j) {
      Ulong buf = a[j];
      Ulong i = j;
      for (; (i >= h) && !inOrder(r[a[i-h]],r[buf]); i -= h)
	a[i] = a[i-h];
      a[i] = buf;
    }
  }
} // |sortI|

/* The permutation computed by |sortI| is called the inverse standardization in
   combinatorics, because it is the inverse of the "standardization"
   permutation, the one that assigns to each position the rank of the
   corresponding entry in the ordering (so any permutation is its own
   standardization).
*/
template <class T, class C> Permutation inverse_standardization
  (const containers::vector<T>& r, C& less)
{
  Permutation pi(r.size());

  /* set the starting value of h */
  Ulong h = 1;
  while (h < r.size()/3)
    h = 3*h+1;

  for (; h > 0; h /= 3)
  { // sort classes modulo |h| of indices in |a|
    for (Ulong j = h; j < pi.size(); ++j)
    {
      const Ulong key = pi[j]; // an index into |r|
      const T& val = r[key];
      // insertion sort |key| into subarray of |a| at |h|-class of |j| upto |j|
      Ulong i = j; // the index into |a| where a value was last moved from
      for ( ; i>=h; i-=h)
      {
	if (less(r[pi[i-h]],val))
	  break; // don't use |i-h|, the place is |i|
	pi[i] = pi[i-h]; // move index |h| places up, leaving |a[i-h]| moved from
      }
      // now |i| is the index into |a| to move |key| to
      pi[i] = key;
    }
  }
  return pi;
} // |inverse_standardization(r,inOrder)|


/*
  General sort function where the comparison is made using a functor f.
  It is assumed that operator> is defined for the valuetype of f.

  Doesn't actually modify r; it only writes down in a the permutation
  s.t. new[j] = old[a[j]].
*/
template <class T, class F> void sortI_f
  (const list::List<T>& r, F& f, Permutation& a)
{
  a.identity(r.size());

  /* set the starting value of h */

  Ulong h = 1;

  for (; h < r.size()/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < r.size(); ++j) {
      Ulong buf = a[j];
      Ulong i = j;
      for (; (i >= h) && (f(r[a[i-h]]) > f(r[buf])); i -= h)
	a[i] = a[i-h];
      a[i] = buf;
    }
  }

  return;
} // |sortI_f|

}; // |namespace bits|

#endif

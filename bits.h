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
#include "list.h"
#include "bitmap.h"

/******** type declarations *************************************************/

namespace bits {
  class BitMap;
  class Partition;
  class PartitionIterator;
  class Permutation;
  class SubSet;
  typedef unsigned char Flags;
  typedef Ulong SetElt;
  typedef list::List<SetElt> Set;
};

/******** function declarations *********************************************/

#include "io.h"

namespace bits {
  unsigned bitCount(const Lflags& f);
  bool isRefinement(const Partition& pi1, const Partition& pi2);
  void memSet(void *dest, void *source, Ulong size, Ulong count);
  void print(FILE* file, const BitMap& map);
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
};

/******** type definitions **************************************************/

#include "constants.h"

class bits::Permutation:public Set
{
 public:
/* constructors and destructors */
  Permutation();
  Permutation(Ulong n); // identity of size |n|
/* manipulators */
  Permutation& identity(Ulong n); // assign identity to |*this| (returned)
  Permutation inverse() const;
  Permutation& compose(const Permutation& a);
  Permutation& rightCompose(const Permutation& a);
};

class bits::BitMap {
 private:
  list::List<Lflags> d_map;
  Ulong d_size;
 public:
/* constructors and destructors */
  BitMap() {};
  BitMap(Ulong n);
  BitMap(const BitMap& map): d_map(map.d_map), d_size(map.d_size) {}
  BitMap(const bitmap::BitMap& bm); // convert format fast
  ~BitMap(); /* standard destructor */
/* modifiers */
  BitMap& operator=(const BitMap& map) { return assign(map); }
  BitMap& assign(const BitMap& map);
  void clearBit(Ulong n)
    { d_map[n/BITS(Lflags)] &= ~(constants::eq_mask[n%BITS(Lflags)]); }
  void permute(Permutation& q);
  void reset() { d_map.setZero(); }
  void setBit(Ulong n)
    { d_map[n/BITS(Lflags)] |= constants::eq_mask[n%BITS(Lflags)]; }
  void setBit(Ulong n, bool t) { if (t) setBit(n); else clearBit(n); }
  void setSize(Ulong n);
/* operations */
  void operator~ ();
  void operator&= (const BitMap& map);
  void operator|= (const BitMap& map);
  void andnot(const BitMap& map);
/* accessors */
  Ulong bitCount() const;
  Lflags chunk(Ulong m) const { return d_map[m]; }
  Ulong firstBit() const;
  bool isEmpty(Ulong m) const;
  Ulong lastBit() const;
  Lflags lastchunk() const
    { return constants::leq_mask[(size()-1)%BITS(Lflags)]; }
  bool getBit(Ulong n) const
    { return d_map[n/BITS(Lflags)] & constants::eq_mask[n%BITS(Lflags)]; }
  Ulong size() const { return d_size; }
/* iterator */
  class Iterator;
  class ReverseIterator;
  friend class Iterator;
  Iterator begin() const;
  Iterator end() const;
  ReverseIterator rbegin() const { return ReverseIterator(end()); }
  ReverseIterator rend() const   { return ReverseIterator(begin()); }

  class Iterator :
    public std::iterator<std::forward_iterator_tag, Ulong>
  { /* is really a constant iterator */
  private:
    static const Lflags posBits = BITS(Lflags) - 1;  /* BITS(Lflags) should be a
							power of two */
    static const Lflags baseBits = ~posBits;
    const BitMap* d_b;
    const Lflags* d_chunk;
    Ulong d_bitAddress;
  public:
    Iterator();
    Iterator(const BitMap& b);
    ~Iterator() = default;
    Ulong bitPos() const     { return d_bitAddress&posBits; }
    Ulong operator* () const { return d_bitAddress; }
    Iterator& operator++ ();
    Iterator& operator-- ();
    bool operator== (const Iterator& i) const
      { return d_bitAddress == i.d_bitAddress; }
    bool operator!= (const Iterator& i) const
      { return d_bitAddress != i.d_bitAddress; }

  /* friend declaration */
  friend Iterator BitMap::end() const;
  }; // |class Iterator|

  class ReverseIterator {
  private:
    Iterator d_i;
  public:
    ReverseIterator() {};
    explicit ReverseIterator(const Iterator& i):d_i(i) {};
    ~ReverseIterator() {};
    Ulong operator* () const { Iterator tmp(d_i); --tmp; return *tmp; }
    ReverseIterator& operator++ () { --d_i; return *this; }
    ReverseIterator& operator-- () { ++d_i; return *this; }
    bool operator== (const ReverseIterator& i) const { return d_i == i.d_i; }
    bool operator!= (const ReverseIterator& i) const { return d_i != i.d_i; }
  }; // |class ReverseIterator|
}; // |class BitMap|

class bits::SubSet {
 private:
  bitmap::BitMap d_bitmap;
  containers::vector<Ulong> row;
 public:
/* constructors and destructors */
  SubSet() {};
  SubSet(Ulong n) : d_bitmap(n), row() {}
  SubSet(const SubSet& q):d_bitmap(q.d_bitmap), row(q.row) {}
/* accessors */
  const containers::vector<Ulong>& elements() const { return row; }
  const Ulong operator[] (Ulong j) const { return row[j]; }
  const bitmap::BitMap& bitMap() const { return d_bitmap; }
  Ulong find(const SetElt& x) const
  { return std::lower_bound(row.begin(),row.end(),x)-row.begin(); }
  bool is_member(Ulong n) const { return d_bitmap.is_member(n); }
  Ulong size() const { return row.size(); }
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
}; // |bits::SubSet|


class bits::Partition
{
private:
  containers::vector<Ulong> classifier;
  Ulong d_classCount;
 public:
/* class definitions */
  typedef Ulong result_type;
/* constructors and destructors */
  Partition() : Partition(0) {}
  Partition(Ulong n) : classifier(n,0), d_classCount(0) {}
  // Partition(const Partition& a, const BitMap& b); // partition of subset
  template <class F> Partition(Ulong n, const F& property);
  template <class I, class F> Partition
    (const I& first, const I& last, const F& f);
  template <class I> Partition
    (I first, Ulong count, const containers::vector<Ulong>& class_sizes);
/* accessors */
  // the next method allows treating the Partition as a (classifying) function
  Ulong operator() (Ulong j) const { return classifier[j]; }
  Ulong classCount() const {return d_classCount;}
  Ulong size() const { return classifier.size(); }

  Permutation standardization() const; // standardization of |classifier|
  Permutation inverse_standardization() const; // its inverse
  SubSet class_nr(Ulong n) const; // equivalence class with label |n|
  SubSet class_of(Ulong x) const { return class_nr(classifier[x]); }
/* modifiers */
  Ulong& operator[] (Ulong j)     { return classifier[j]; }
  Partition& normalize() { return *this = Partition(size(),*this); }
  void permute(const Permutation& a);
  void permuteRange(const Permutation& a);
  void setClassCount();
  void setClassCount(Ulong count) { d_classCount = count; }
  void setSize(Ulong n) { classifier.resize(n); }
/* input/output */
  void printClassSizes(FILE* file) const;
}; // |bits::Partition|

class bits::PartitionIterator {
  const Partition& d_pi;
  Permutation d_a;
  Set d_class;
  Ulong d_base;
  bool d_valid;
 public:
/* constructors and destructors */
  PartitionIterator(const Partition& pi);
/* iterator operations */
  operator bool() const { return d_valid; }
  void operator++();
  const Set& operator()() const { return d_class; }
}; // |bits::PartitionIterator|

/**** Inline implementations **********************************************/

/******** template definitions ***********************************************/

namespace bits {


/*
  This constructor defines the partition of [0,r.size()[ defined by f; f
  is supposed to be a function taking arguments of type T. It is also
  assumed that operator<= is defined for the value type of f (so that
  the function insert may be applied.)
*/
template <class F> Partition::Partition
  (Ulong n, const F& property) : classifier(), d_classCount(0)
{
  std::map<decltype(property(0)),Ulong> renumber;

  for (Ulong i=0; i<n; ++i)
    if (renumber.insert(std::make_pair(property(i),d_classCount)).second)
      ++d_classCount; // |property(i)| was inserted (new key); consolidate value

  classifier.reserve(n);
  for (Ulong i=0; i<n; ++i)
    classifier.push_back(renumber[property(i)]);
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
  (const I& first, const I& last, const F& f) : classifier(), d_classCount(0)
{
  std::map<decltype(f(*first)),Ulong> renumber;

  Ulong size=0;
  for (auto it=first; it!=last; ++it,++size)
    renumber.insert(std::make_pair(f(*it),0)); // form |renumber| as a mere set

  for (auto& p : renumber)
    p.second = d_classCount++; // fill in values increasingly

  classifier.reserve(size);
  for (auto it=first; it!=last; ++it)
    classifier.push_back(renumber[f(*it)]);
}

template <class I> Partition::Partition
  (I first, Ulong count, const containers::vector<Ulong>& class_sizes)
    : classifier(count,-1)
  , d_classCount(class_sizes.size())
{
  classifier.reserve(count);
  for (Ulong class_nr = 0; class_nr<d_classCount; ++class_nr)
    for (Ulong i = 0; i<class_sizes[class_nr]; ++i)
      classifier[*first++]=class_nr;
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
  BitMap b(r.size());

  for (Ulong j = 0; j < a.size(); ++j) {

    if (b.getBit(j))
      continue;

    if (a[j] == j) {
      b.setBit(j);
      continue;
    }

    Ulong k = j;
    b.setBit(j);

    for (Ulong i = a[j]; i != j; i = a[i]) {
      T buf = r[k];
      r[k] = r[i];
      r[i] = buf;
      k = i;
      b.setBit(i);
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

};

#endif

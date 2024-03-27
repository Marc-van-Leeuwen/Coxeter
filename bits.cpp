/*
  This is bits.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "bits.h"

#include <limits.h>
#include <algorithm>

#include "bitmap.h"
#include "io.h"


/*****************************************************************************

  This module regroups classes and functions for bit-access of structures.
  Subsets of given structures will frequently in this program be represented
  by bitmaps; these may be small bitmaps, fitting inside a single Ulong,
  or larger ones, necessitating a list.

  In this program, a bitmap will always be a list of unsigned chars, ensuring
  portability across all platforms, independently of word size and endianism
  problems. However, to accelerate bitwise operations we always make sure
  that bitmaps are allocated in multiples of sizeof(Ulong) characters,
  so that bitwise operations can be done one Ulong at a time.

 *****************************************************************************/

/*****************************************************************************

        Chapter I -- The Permutation class

  A permutation of a set is just the enumeration of the elements of the set
  (assumed to be running from 0 to N-1) in a certain order; a[j] is
  interpreted as the image of element #j under a.

  We have tried to be consistent in the whole program about this, and the
  way permutations acts on sets and functions; when a acts on the range,
  we somply apply a to the image of a function; when a acts on the domain,
  we compose f with a^{-1} (we are able to do this essentially in-place.)

  The following functions are defined :

    - constructors and destructors :

      - Permutation();
      - Permutation(n);
      - ~Permutation();

    - accessors

      - inverse() : inverse of the permutation;

    - manipulators :

      - compose(a) : increment by a under right composition;
      - leftCompose(a) : same on the left;

 *****************************************************************************/

namespace bits {

  Permutation::Permutation() : containers::vector<Ulong>()
{}

Permutation::Permutation(Ulong n)
  : containers::vector<Ulong>()
{ identity(n);
}



// Set |*this| to the identity permutation of [0,n[.
Permutation& Permutation::identity(Ulong n)
{
  clear();
  reserve(n);

  for (Ulong i = 0; i < n; ++i)
    push_back(i);

  return *this;
}


// Inverse of |*this| : new(old(x)) = x.
Permutation Permutation::inverse() const
{
  Permutation result(size());

  const Permutation& t = *this; // for more readable subscription

  for (SetElt x = 0; x < size(); ++x)
    result[t[x]] = x;

  return result;
}


/*
  Compose |*this| on the right with |a|: |(*this) *= a|
  new(x) = old(a(x)). The same problem occurs as with inverse.
*/

Permutation& Permutation::rightCompose(const Permutation& a)
{
  Permutation c; c.reserve(size());
  Permutation& t = *this;

  for (SetElt x = 0; x < size(); ++x)
    c.push_back(t[a[x]]);

  return t = std::move(c);
}


/*
  Compose |*this| on the left with |a|: |(*this) = a*(*this)|
  new(x) = a(old(x)).
*/
Permutation& Permutation::compose(const Permutation& a)
{
  Permutation& t = *this;
  for (SetElt x = 0; x < size(); ++x)
    t[x] = a[t[x]];

  return *this;
}


};

/****************************************************************************

        Chapter IV -- The Partition class

  A partition of a set is represented by a function of the set to the range
  [0,N[, where N is the number of classes in the partition. Moreover it
  is useful to have direct access to the number of classes.

  The following functions are defined :

   - constructors and destructors :

     - Partition();
     - Partition(n) : constructs a partition with a list of size n;
     - ~Partition() (inlined);

   - accessors :

     - classCount() : returns the number of classes; (inlined)
     - sort(v) : returns a permutation vector of the set, so that classes
       are contiguous and sorted in their original order;
     - writeClass(b,n) : sets b to hold class #n;

   - modifiers :

     - normalize() : normalizes the permutation;
     - normalize(a) : writes the normalizing permutation in a;
     - permute(a) : permutes the partition according to a;
     - setClassCount(n) : sets the number of classes to n;
     - setClassCount() : fill in the number of classes;

   - input/output :

     - printClassSizes(file) : prints the sizes of the classes;

  NOTE : this class really has nothing to do with bits. It should probably
  be moved to another module (maybe sort.h, or sets.h ?)

 ****************************************************************************/

namespace bits {

/******* accessors **********************************************************/

Partition::Partition(containers::vector<Ulong>&& v, Ulong class_count)
: classifier(std::move(v))
, inv_stnd() // filled by below
, cum_class_sizes() // set by |set_inverse_standardization|
{ set_inverse_standardization(class_count); }

/*
  Return the permutation vector to a partition with the same class sizes, but
  for which the classes are contiguous, in increasing order, and each class is
  kept in the enumeration order of the original set. In other words, |new| is
  the weakly increasing stable sorting of the |old| elements, the |a| is such
  that $new[a[j]] = old[j]$ for all $j$ and $a$ is increasing on the positions
  of each equivalence class. This is called the standardization of the
  classifying vector.

  We do this by counting each class, then putting each element directly
  in its right place in a.
*/
Permutation Partition::standardization() const
{
  containers::vector<Ulong> class_counter(class_count(),0);

  for (Ulong class_nr : classifier)
    ++class_counter[class_nr];

  // cumulate left-to-right
  Ulong sum=0;
  for (Ulong& count : class_counter)
  { std::swap(count,sum); // current sum replaces count
    sum += count; // which get added to sum
  }

  Permutation result(size());

  for (Ulong j = 0; j < size(); ++j)
    result[j] = class_counter[classifier[j]]++;

  return result;
}

/*
  Compute the inverse standardization permutation into |inv_stnd|, while setting
  |cum_class_sizes| to the cumulative lass sizes (entry |i| being up to and
  including calls |i|). This permutation being more useful tham its inverse,
  which is the standardization the, |Partition| these dedicated fields to store
  those values; they allow easy traversal of each individial class.

  Here $new[j] = old[a[j]]$ for the same |old| and |new| as in |standardization|.
  The implementation is the same as for |standardization|, except for the body
  of the final loop, where we interchange index to |result| and assigned value.
*/
void Partition::set_inverse_standardization(Ulong class_count)
{
  cum_class_sizes.assign(class_count,0);
  for (Ulong class_nr : classifier)
    ++cum_class_sizes[class_nr];

  // cumulate left-to-right, giving starting index of each class (for now)
  Ulong sum=0;
  for (Ulong& count : cum_class_sizes)
  { std::swap(count,sum); // current sum replaces count
    sum += count; // which get added to sum
  }

  if (inv_stnd.size()!=classifier.size())
    inv_stnd.assign(classifier.size(),-1);
  // traverse |classifier| advancing each class individually
  for (Ulong j = 0; j < size(); ++j)
    inv_stnd[cum_class_sizes[classifier[j]]++] = j;
  // now |cum_class_sizes| entries indicate end index of each class
  assert(class_count==0 or cum_class_sizes.back()==sum);
} // ||

bool Partition::refine(const containers::vector<Partition>& L)
{
  assert(L.size()==class_count());
  containers::vector<Ulong> class_counter(class_count(),0);
  for (Ulong class_nr=0; class_nr<L.size(); ++class_nr)
    class_counter[class_nr] += L[class_nr].class_count()-1; // count new classes

  // now cumulate those values starting from |class_count()|
  Ulong sum=class_count();
  for (auto& entry : class_counter)
  { std::swap(sum,entry); // store value of |sum| before addition
    sum += entry;
  }
  if (sum==class_count())
    return false; // no actual refinement

  containers::vector<Ulong> class_index(class_count(),0); // class-relative pos

  for (Ulong& class_nr : classifier)
  { const auto i = class_index[class_nr]++;
    auto sub_class = L[class_nr].classifier[i];
    if (sub_class>0)
      class_nr = class_counter[class_nr] + sub_class-1; // assign new class
  }

  assert(std::all_of(class_index.begin(),class_index.end(),
		     [&L,&class_index](const Ulong& index) -> bool
		     { return L[&index-&class_index[0]].size()==index; }));
  set_inverse_standardization(sum); // just reconstructing is easiest
  return true;
}


SubSet Partition::class_nr(Ulong i) const
{
  return SubSet(containers::vector<Ulong>(class_bound(i),class_bound(i+1)),size());
}


/******* modifiers **********************************************************/


/*
  Permute set underlying the partition according to |a| (i.e., apply |a| to the
  elements of each part of the partition)
*/
void Partition::permute_base(const Permutation& a)
{
  // modifying |inv_stnd| is straightforward (modulo order withing classes)
  for (auto& entry : inv_stnd)
    entry = a[entry];

  // sort each class, and now adapt |classifier| to permuted positions
  for (Ulong i=0; i<class_count(); ++i)
  { std::sort(class_bound(i),class_bound(i+1));
    for (auto it=class_bound(i); it!=class_bound(i+1); ++i)
      classifier[*it] = i;
  }
}

// Apply the permutation |a| to the values of the classifying function.
void Partition::permute_range(const Permutation& a)
{
  for (SetElt x = 0; x < size(); ++x)
    classifier[x] = a[classifier[x]];

  set_inverse_standardization(class_count());
}


containers::vector<Ulong> Partition::class_sizes() const
{
  containers::vector<Ulong> result; result.reserve(class_count());
  for (Ulong i=0; i<class_count(); ++i)
    result.push_back(class_size(i));
  return result;
}


/****************************************************************************

        Chapter VI -- The SubSet class.

  The SubSet class is designed to hold subsets of enumerated sets. Typically,
  we will want to know which elements are in the subset, and also, to have
  a quick way of deciding whether an arbitrary element is in it.

  We implement this using a list of the elements in the subset, and a bitmap
  describing the subset within the big set. Of course everything could be
  done just from the bitmap, but this would slow down traversal somewhat.
  (Otherwise we would have to write an iterator directly from the bitmap,
  doubtless a worthwile exercise.)

  The following functions are defined :

   - constructors and destructors :

     - SubSet(n) : constructs a subset with a bitmap of size n; (inlined)
     - ~SubSet() : standard destructor (inlined);

   - accessors :

     - isMember(n) : tells if element #n is a member of the subset; (inlined)
     - size() : number of elements in the subset; (inlined)

   - modifiers :

     - add(n) : adds element n to the subset;
     - readBitMap() : reads the contents of the bitmap into the list;
     - reset() : resets the subset to the empty subset;


 *****************************************************************************/


/*
  Add a new element to the subset. It is assumed that n is a legal value
  w.r.t. the bitmap. We do not sort the elements in order; some special
  function should take care of that if required.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.
*/

void SubSet::add(Ulong n)
{
  if (d_bitmap.is_member(n)) /* n is already in there */
    return;

  d_bitmap.insert(n);
  row.push_back(n);

  /* the error OUT_OF_MEMORY may have been set here */
}


/*
  Puts the content of the bitmap in the list in a simple-minded way.
*/
void SubSet::readBitMap()
{
  row.resize(d_bitmap.size());

  Ulong i = 0;
  for (auto e : d_bitmap)
    row[i++] = e;
}

void SubSet::reset()
{
  d_bitmap.reset();
  row.resize(0);
}

};

/*****************************************************************************

        Chapter VII -- Counting bits.

  This section contains functions for counting bits in bitmaps :

  - bitCount(f) : counts the number of set bits in an Lflags;

 *****************************************************************************/


unsigned bits::bitCount(const Lflags& d_f)

/*
  Returns the number of set bits in f.
*/

{
  unsigned count = 0;

  for (Lflags f = d_f; f; f &= f-1)
    count++;	/* see K&R */

  return count;
}


/*****************************************************************************

        Chapter IX -- Input/Output.

  This section contains i/o functions for the classes defined in this
  module :

   - append(l,map) : append the BitMap map to the string l;
   - print(file,map) : prints the map to the file;

 *****************************************************************************/

namespace bits {


};

/*****************************************************************************

        Chapter X -- Utilities.

  This section contains various utility functions :

   - isRefinement(pi1,pi2) : tells whether pi1 is a refinement of pi2;

 *****************************************************************************/

namespace bits {


/*
  Whether pi1 is a refinement of pi2. Both are assumed to be partitions
  of the same range; the condition is that pi2 should be constant on the
  classes of pi1.
*/

bool isRefinement(const Partition& pi1, const Partition& pi2)
{
  for (Partition::iterator it=pi1.begin(); it!=pi1.end(); ++it)
  {
    const auto range = *it;
    Ulong a = pi2(*range.begin());
    for (auto jt=range.begin()+1; jt!=range.end(); ++jt)
      if (pi2(*jt) != a)
	return false;
  }

  return true;
}

};

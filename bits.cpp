/*
  This is bits.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "bits.h"
#include "bitmap.h"
#include "io.h"

#include <limits.h>

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

Permutation::Permutation():List<Ulong>()

{}

Permutation::Permutation(Ulong n)
  : List<Ulong>()
{ identity(n);
}


Permutation& Permutation::identity(Ulong n)

/*
  Sets the permutation to the identity permutation of [0,n[.
*/

{
  setSize(n);

  for (Ulong j = 0; j < size(); ++j) {
    d_ptr[j] = j;
  }

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

Permutation& Permutation::rightCompose(const Permutation& a)

/*
  Increments the current permutation by composition on the right with a :
  new(x) = old(a(x)). The same problem occurs as with inverse.
*/

{
  static Permutation c(0);

  c.setSize(size());
  Permutation& t = *this;

  for (SetElt x = 0; x < size(); ++x)
    c[x] = t[a[x]];

  t.assign(c);

  return t;
}

Permutation& Permutation::compose(const Permutation& a)

/*
  Increments the current permutation by composition on the left with a :
  new(x) = a(old(x)).
*/

{
  for (SetElt x = 0; x < size(); ++x)
    operator[](x) = a[operator[](x)];

  return *this;
}

};

/*****************************************************************************

        Chapter II -- The BitMap class.

  The following functions are defined for BitMap :

  - constructors :

    - BitMap(Ulong n);
    - Bitmap(const BitMap&); // copy constructor
    - Bitmap(const bitmap::BitMap& bm); // conversion
    - ~BitMap() : standard destructor;

  - accessors :

    - bitCount() : returns the number of set bits;
    - firstBit() : returns the bit-address of the first set bit;
    - getBit(n) : checks if bit n is set;
    - isEmpty(m) : checks if [m,->[ is empty;
    - lastBit() : returns the bit-address of the last set bit;

  - modifiers :

    - assign(map) : sets the bitmap equal to map, resizing if necessary;
    - clearBit(n) : clears bit n; (inlined)
    - permute(q) : applies q to the bitmap;
    - setBit(n) : sets bit n; (inlined)
    - setSize(n) : resizes the bitmap to size n;

  - operations :

    - operator~ () : changes into opposite bitmap;
    - operator&= (map) : intersects with map;
    - operator|= (map) : does union with map;
    - andnot(map) : intersects with the negative of map;

 *****************************************************************************/

namespace bits {

/********** constructors and destructors *************************************/


/*
  Constructor for the BitMap class; constructs a bitmap capable of
  holding n bits.
*/
BitMap::BitMap(Ulong n)
  :d_map(n/BITS(Lflags)+(bool)(n%BITS(Lflags))), d_size(n)
{
  d_map.setSize(n/BITS(Lflags)+(bool)(n%BITS(Lflags)));
}

BitMap::BitMap(const bitmap::BitMap& bm)
  : BitMap(bm.capacity())
{
  constexpr auto chunk = BITS(Lflags);
  const size_t limit = d_size/chunk;
  size_t ii=0;
  for (size_t i=0; i<limit; ++i,ii+=chunk)
    d_map[i] = bm.range(ii,chunk);
  if (bm.capacity()>ii)
    d_map[limit] = bm.range(ii,bm.capacity()-ii);
}

/*
  No memory is directly allocated by BitMap.
*/
BitMap::~BitMap()
{}

/********** accessors ********************************************************/

Ulong BitMap::bitCount() const

/*
  Returns the number of set bits in the bitmap.
*/

{
  Ulong count = 0;

  for (Ulong j = 0; j < d_map.size(); ++j)
    count += bits::bitCount(d_map[j]);

  return count;
}

Ulong BitMap::firstBit() const

/*
  Returns the bit position of the first set bit.

  NOTE : this looks totally wrong!
*/

{
  Ulong first = 0;
  Lflags f = (Lflags)1;

  for (Ulong j = 0; j < d_map.size(); ++j) {
    if (d_map[j]) { /* first bit found */
      f = d_map[j];
      return f;
    }
    else
      first += BITS(Lflags);
  }

  return first + constants::firstBit(f);
}

bool BitMap::isEmpty(Ulong m) const

/*
  This function checks whether the intersection of the bitmap with the
  interval [m,size[ is empty
*/

{
  Ulong lsize = d_size/BITS(Lflags)+(bool)(d_size%BITS(Lflags));

  /* look at word containing m */

  Ulong ml = m/BITS(Lflags);
  Ulong mr = m%BITS(Lflags);
  Ulong mc = BITS(Lflags)-1 - mr;
  Lflags f = constants::leq_mask[mc] << mr;

  if (d_map[ml]&f)
    return(false);

  for (Ulong j = ml+1; j < lsize; ++j) {
    if (d_map[j])
      return false;
  }

  return true;
}

Ulong BitMap::lastBit() const

/*
  This function returns the bit-address of the last set bit in b. The
  return value is b.size() if b is empty.
*/

{
  if (d_size == 0)
    return 0;

  Ulong base = (d_size-1)/BITS(Lflags)+1;

  while(base) {
    base--;
    Lflags f = d_map[base];
    if (f)
      return (base*BITS(Lflags)+constants::lastBit(f));
  }

  /* if we reach this point, the bitmap is empty */

  return (d_size);
}

/********** modifiers ********************************************************/

BitMap& BitMap::assign(const BitMap& map)

/*
  Copies the content of map into the current map.
*/

{
  d_map.assign(map.d_map);
  d_size = map.d_size;

  return *this;
}


void BitMap::permute(Permutation& q)

/*
  This function applies the permutation q to the bitmap b. Here b is
  interpreted as a bool-valued table; we want that b[q[x]] hold the
  value previously held by b[x]. This is a range-permutation, as
  explained in kl.cpp.
*/

{
  static BitMap b(0);

  b.setSize(q.size());
  b.reset();

  for (Ulong i = 0; i < d_size; ++i) {
    if (b.getBit(i))
      continue;
    for (Ulong j = q[i]; j != i; j = q[j]) {
      bool t = getBit(j);
      setBit(j,getBit(i));
      setBit(i,t);
      b.setBit(j);
    }
    b.setBit(i);
  }

  return;
}


void BitMap::setSize(Ulong n)

/*
  Resizes the bitmap to hold n bits. If the size grows, it is guaranteed
  that the new bits are set to zero.
*/

{
  d_map.setSize(n/BITS(Lflags) + (bool)(n%BITS(Lflags)));

  if (n > size()) { /* set new bits to zero  */
    Ulong f = size()/BITS(Lflags);  /* word holding first new bit */
    Ulong fb = size()%BITS(Lflags); /* bit address of first new bit in f */
    Lflags old = ((1L << fb) - 1L);     /* flags old bits */
    d_map[f] &= old;
    d_map.setZero(f+1,d_map.size()-f-1);
  }

  d_size = n;
}

/********** operators ********************************************************/

void BitMap::operator~ ()

/*
  Transforms the bitmap into its complement. One has to be careful not to
  exceed the size of the bitmap!
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] = ~d_map[j];

  d_map[d_map.size()-1] &= lastchunk();

  return;
}

void BitMap::operator&= (const BitMap& map)

/*
  Does the bitwise intersection with map. It is assumed that map has at
  least size size().
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] &= map.chunk(j);
  return;
}

void BitMap::operator|= (const BitMap& map)

/*
  Does the bitwise union with map. It is assumed that map has at
  least size size().
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] |= map.chunk(j);
  return;
}

void BitMap::andnot(const BitMap& map)

/*
  Does the bitwise intersection with ~map. It is assumed that map has at
  least size size().
*/

{
  for (Ulong j = 0; j < d_map.size(); ++j)
    d_map[j] &= ~(map.chunk(j));
  return;
}

};

/****************************************************************************

        Chapter III -- The Iterator class.

  This is my first attempt at an STL-style iterator. I'm not trying to
  stick exactly to the stl-notation and requirements, but I hope to be
  in the right spirit. The Iterator class for bitmaps is bidirectional;
  it traverses the set bits for the bitmap.

  The only constructors are for begin() (an iterator pointing to the first
  non-zero element) and end (a non-dereferenceable, past-the-end iterator.)
  It is guaranteed that begin() == end() if the bitmap is empty (i.e.
  bitCount returns zero.)

  The data in the Iterator class have the following meaning. Recall that
  a BitMap is implemented as a list of Lflags. There is a current such
  LFlag, which is pointed by d_chunk; d_bitAddress is the bit address
  of the current set bit. The past-the-end iterator is the one with
  bitAddress equal to the size of the bitmap.

  The following functions are defined :

   - Iterator() construct iterator with null pointer
   - Iterator(const BitMap& m) : constuct |m.begin()|
   - operator* () : returns the position of the current bit;
   - operator++ () : moves to the next set bit, or past-the-end;
   - operator-- () : moves to the previous set bit (valid if iterator is
     not begin());
   - operator== (i) : says if the two iterators point to the same bit-address;
     (inlined);
   - operator!= (i) : the negation of operator==.

 ****************************************************************************/

namespace bits {

BitMap::Iterator::Iterator()
  : d_b(nullptr)
  , d_chunk(nullptr)
  , d_bitAddress(0)
{}

BitMap::Iterator::Iterator(const BitMap& b)
  : d_b(&b)
  , d_chunk(d_b->d_map.ptr())
  , d_bitAddress(0)
{ // do essentially |operator++|, but without skipping the initial set bit
  for ( ; d_bitAddress < d_b->size(); d_bitAddress += BITS(Lflags),++d_chunk)
    if (*d_chunk)
    {
      d_bitAddress += constants::firstBit(*d_chunk);
      break;
    }
  if (d_bitAddress > d_b->size()) /* bitmap was empty */
    d_bitAddress = d_b->size();
}


/*
  Increment operator (prefix notation). Valid if the iterator is
  dereferenceable. Result is dereferenceable or past-the-end.
*/
BitMap::Iterator& BitMap::Iterator::operator++ ()
{
  Lflags f = *d_chunk >> bitPos();
  f >>= 1;

  if (f!=0) // bit found in current chunk
    d_bitAddress += constants::firstBit(f)+1; //
  else { /* go to next chunk */
    d_bitAddress &= baseBits; // back up to notch, |for| starts advancing notch
    ++d_chunk;
    for (d_bitAddress += BITS(Lflags) ; d_bitAddress < d_b->size();
	 d_bitAddress += BITS(Lflags), ++d_chunk)
    {
      if (*d_chunk) {
	d_bitAddress += constants::firstBit(*d_chunk);
	break;
      }
    }
    if (d_bitAddress > d_b->size())
      d_bitAddress = d_b->size();
  }

  return *this;
}

BitMap::Iterator& BitMap::Iterator::operator-- ()

/*
  Decrement operator (prefix notation). Valid if the iterator is past-the-end,
  or is dereferenceable and not equal to begin().
*/

{
  Lflags f = 0;

  if (bitPos()) {
    f = *d_chunk & constants::lt_mask[bitPos()];
  }

  if (f) {
    d_bitAddress &= baseBits;
    d_bitAddress += constants::lastBit(f);
  }
  else { /* go to previous chunk */
    d_bitAddress &= baseBits;
    while (d_bitAddress) {
      d_bitAddress -= BITS(Lflags);
      --d_chunk;
      if (*d_chunk) {
	d_bitAddress += constants::lastBit(*d_chunk);
	break;
      }
    }
  }

  return *this;
}

BitMap::Iterator BitMap::begin() const

/*
  Returns
*/

{
  static Iterator i;
  new(&i) Iterator(*this);
  return i;
}

BitMap::Iterator BitMap::end() const

{
  static Iterator i;

  i.d_b = this;
  i.d_bitAddress = d_size;
  i.d_chunk = d_map.ptr()+d_map.size();
  if (i.bitPos())
    i.d_chunk--;

  return i;
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
  containers::vector<Ulong> class_counter(d_classCount,0);

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
  Return the inverse standardization permutation directly. This is in fact more
  often useful then the standardization, because the inverse standardization
  permutation permits easy traversal of all values in weakly increasing order,
  for instance finding the classes of indices with equal values in a partition.

  Now $new[j] = old[a[j]]$ for the same |old| and |new| as in |standardization|.
  The implementation is the same as for |standardization|, except for the body
  of the final loop, where we interchange index to |result| and assigned value.
*/
Permutation Partition::inverse_standardization() const
{
  containers::vector<Ulong> class_counter(d_classCount,0);

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
    result[class_counter[classifier[j]]++] = j;

  return result;
} // |inverse_standardization|

SubSet Partition::class_nr(Ulong n) const
{
  SubSet result(size());
  for (Ulong j = 0; j < size(); ++j)
    if (classifier[j] == n)
      result.add(j);

  return result;
}


/******* modifiers **********************************************************/


/*
  Permute set underlying the partition according to |a| (i.e., apply |a| to the
  elements of each part of the partition)
*/
void Partition::permute(const Permutation& a)
{
  bitmap::BitMap seen(size());

  for (SetElt x = 0; x < size(); ++x)
  {
    if (not seen.is_member(x))
      for (SetElt y = a[x]; y != x; y = a[y])
      { std::swap(classifier[x],classifier[x]);
	seen.insert(y);
      }
    // we can do without |seen.insert(x)|: we cannot come here again
  }
}


// Apply the permutation |a| to the values of the classifying function.
void Partition::permuteRange(const Permutation& a)
{
  for (SetElt x = 0; x < size(); ++x)
    classifier[x] = a[classifier[x]];

  return;
}

void Partition::setClassCount()
{
  Ulong count = 0;

  for (Ulong j = 0; j < size(); ++j) {
    if (classifier[j] >= count)
      count = classifier[j]+1;
  }

  d_classCount = count;

  return;
}

/******** input/output ******************************************************/


/*
  This function prints out the sizes of the classes in the partition.
*/
void Partition::printClassSizes(FILE* file) const
{
  static list::List<Ulong> count(0);

  count.setSize(d_classCount);
  count.setZero();

  for (Ulong j = 0; j < size(); ++j) {
    count[classifier[j]]++;
  }

  for (Ulong j = 0; j < d_classCount; ++j) {
    fprintf(file,"%lu",count[j]);
    if (j < d_classCount-1)
      fprintf(file,",");
  }

  fprintf(file,"\n");

  return;
}

};

/****************************************************************************

        Chapter V -- The PartitionIterator class.

  This class is intended for the convenient traversal of a partition. At
  each iteration, a new class is provided as a list (maybe this should
  become a SubSet ?). To do this conveniently, a permutation d_a is used
  to sort the partition by contiguous classes. The current class is
  kept in d_class.

 ****************************************************************************/

namespace bits {

PartitionIterator::PartitionIterator(const Partition& pi)
  : d_pi(pi)
  , d_a(pi.size())
  , d_class(0)
  , d_base(0)
  , d_valid(pi.size()!=0)
{
  if (d_valid)
  {
    d_a = pi.inverse_standardization();

    /* load first class */

    Ulong j = 0;

    for (; (j < d_a.size()) && (d_pi(d_a[j]) == d_pi(d_a[d_base])); ++j) {
      d_class.append(d_a[j]);
    }
  }
} // |PartitionIterator::PartitionIterator|


void PartitionIterator::operator++ ()

{
  d_base += d_class.size();

  if (d_base == d_pi.size()) {
    d_valid = false;
    return;
  }

  Ulong j = d_base;
  d_class.setSize(0);

  for (; (j < d_a.size()) && (d_pi(d_a[j]) == d_pi(d_a[d_base])); ++j) {
    d_class.append(d_a[j]);
  }
}

};

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


namespace bits {


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

        Chapter VIII -- Copying memory.

  This section contains functions for copying memory between bitmaps :

  - MemSet(dest,source,size,count) : copies into dest count repetitions
    of the pattern made up by the first size bits in source;

 *****************************************************************************/


void bits::memSet(void *dest, void *source, Ulong size, Ulong count)

/*
  Copies into dest count repetitions of the pattern made up by the first size
  bits in source.
*/

{
  Ulong c;

  if (count == 0)
    return;

  memmove(dest,source,size);
  source = dest;
  dest = (void *)((char *)dest + size);

  for (c = 1; c <= count/2; c *= 2)
    {
      memmove(dest,source,c*size);
      dest = (void *)((char *)dest + c*size);
    }

  memmove(dest,source,(count-c)*size);

  return;
}


/*****************************************************************************

        Chapter IX -- Input/Output.

  This section contains i/o functions for the classes defined in this
  module :

   - append(l,map) : append the BitMap map to the string l;
   - print(file,map) : prints the map to the file;

 *****************************************************************************/

namespace bits {


// append to |l| a string on the alphabet |{0.1}| for the bitmap
std::string& append(std::string& l, const BitMap& map)
{
  for (Ulong j = 0; j < map.size(); ++j)
    l.push_back(map.getBit(j) ? '1' : '0');

  return l;
}

void print(FILE* file, const BitMap& map)
{
  static std::string buf;

  buf.clear();
  append(buf,map);
  io::print(file,buf);

  return;
}

};

/*****************************************************************************

        Chapter X -- Utilities.

  This section contains various utility functions :

   - isRefinement(pi1,pi2) : tells whether pi1 is a refinement of pi2;

 *****************************************************************************/

namespace bits {

bool isRefinement(const Partition& pi1, const Partition& pi2)

/*
  Tells whether pi1 is a refinement of pi2. Both are assumed to be partitions
  of the same range; the condition is that pi2 should be constant on the
  classes of pi1.
*/

{
  for (PartitionIterator i(pi1); i; ++i) {
    const Set& l = i();
    Ulong a = pi2(l[0]);
    for (Ulong j = 1; j < l.size(); ++j)
      if (pi2(l[j]) != a)
	return false;
  }

  return true;
}

};

/*
  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006,2017,2022 Marc van Leeuwen
  part of the Coxeter, and of Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

// Definitions and declarations for the BitMap class.

#ifndef BITMAP_H  /* guard against multiple inclusions */
#define BITMAP_H

#include <vector>
#include <cassert>

#include <limits>
#include <cstdint> // for |uint32_t| and |uint64_t|

#include "constants.h"
#include "containers.h"

namespace bitmap {

  class BitMap;

/******** type definitions **************************************************/

/*
  Container of a large (more than twice the machine word size) set of bits.

  From the point of view of a user of the class, a |BitMap| should be seen as a
  container of |size_t| values, not bits: these unsigned values are the
  addresses of the set (value=1) bits. When the class is used for example to
  flag the noncompact roots in a set of roots, it is most convenient to think of
  it as containing the numbers of the noncompact roots (on a fixed list of all
  roots). The class obeys the semantics of a Forward Container (from the C++
  standard library). Its "size" as a container is the number of |size_t| values
  that it "contains"; that is, the number of bits set to 1 in the bitmap.

  The basic data is in |d_map|, a vector of |chunk_type| integers. Each of
  these integers is a "chunk" (of size |longBits|, presumably the machine word
  length) of the bit map. The capacity (in bits) of the |BitMap| is
  |d_capacity|; the size of the vector |d_map| is |d_capacity/longBits| (plus
  one if there is a remainder in the division).

  We provide bit-addressed access to this map, in other words given a value it
  is fast to test or change its membership status. Also we wish to define an
  iterator class, which traverses the _set_ bits of the bitmap; so that for
  instance, b.begin() would access the first set bit. Dereferencing the iterator
  yields its bit-address, a |size_t| value.
*/
class BitMap
{
  using chunk_type = std::uint64_t;
  static constexpr size_t index_shift       = 6; // between bit and chunk address
  static constexpr size_t chunk_bits        = 1ul<<index_shift; // |0x40 == 64|
  static constexpr size_t position_bits     = chunk_bits-1;     // |0x3F|
  static constexpr size_t chunk_index_bits  = ~position_bits;
  using chunk_iterator = containers::vector<chunk_type>::iterator;
  using chunk_const_iterator = containers::vector<chunk_type>::const_iterator;

  size_t d_capacity;
  containers::vector<chunk_type> d_map;

 public:

// type definitions

  using value_type = size_t; // type of values virtually present in a |BitMap|
  typedef value_type& reference;
  typedef const value_type& const_reference;
  typedef value_type* pointer;
  typedef const value_type* const_pointer;
  typedef ptrdiff_t difference_type;
  using size_type = size_t; // container capacity same type as |value_type|

// iterators
  class const_iterator;
  class iterator;

// constructors and destructors
 BitMap() : d_capacity(0), d_map() {} // create a bitmap without capacity

  /*
    Construct a zero-initialized bitmap with a capacity of n bits.

    Note that the size of the vector |d_map| exceeds |n >> index_shift| by
    one, unless |longBits| exactly divides |n|.
  */
  explicit BitMap(size_t n)
    : d_capacity(n), d_map((n+position_bits)>>index_shift,0)
  {
    assert (d_capacity==n); // that is, conversion to |U| did not lose any bits
  }

  // convert range defined by iterators into a bitmapped_set
  template <typename I> BitMap(size_t n, const I& first, const I& last)
  : d_capacity(n)
  , d_map((d_capacity+position_bits)>>index_shift,0)
  {
    for (I it = first; it!=last; ++it)
      insert(*it);
  }

  // an easier version for a full vector
  template <typename U> // unsigned integral type
  BitMap(size_t n,const containers::vector<U>& v)
    : d_capacity(n), d_map((n+position_bits)>>index_shift)
  {
    for (size_t i=0; i<v.size(); ++i)
      insert(v[i]);
  }

  // Set of offsets into [first,last[ of values (also) found in [fsub,lsub[
  template <typename I, typename J>
    BitMap(const I& first, const I& last, const J& fsub, const J& lsub);


// assignment: don't ask don't tell; everything will be generated implicitly. But
// declaring (defaulted) move versions would implicitly remove all copy versions
// unless those copy versions are also declared (defaulted); just don't bother.

// accessors
  /*
   Number of bits in use in the bitmap. This is the capacity
   of the BitMap as a standard library container, not d_map.size(), which is
   approximately longBits times smaller.
  */
  size_t capacity() const { return d_capacity; }
  size_type size() const; // the number of bits that are set in the bitmap

  // fast comparisons for containers and other standadrd library purposes
  bool operator< (const BitMap& b) const
  { return
      d_capacity==b.d_capacity ? d_map < b.d_map : d_capacity<b.d_capacity ; }
  bool operator== (const BitMap& b) const
  { return d_capacity==b.d_capacity and d_map == b.d_map; }
  bool operator!=(const BitMap& b) const
  { return d_capacity!=b.d_capacity or d_map != b.d_map; }

  bool equivalent (const BitMap& b) const; // inclusion both ways, no capacity

  bool empty() const; // whether |size()==0|
  bool full() const; // whether |size()==capacity()|
  bool none() const { return empty(); } // naming compatibility with |BitSet|
  bool any() const { return not empty(); }

  // Whether bit |n| in the bitmap is set (whether |n| is a member of the set)
  bool is_member(size_t n) const
  {
    if (n >= d_capacity)
      return false;
    return (d_map[n >> index_shift] & constants::eq_mask[n&position_bits])!=0;
  }

  // Value at index |n| if viewed as list of |size_t| values
  size_t n_th(size_t n) const;

  // Number of values |<n| present (set) in the bitmap
  size_t position(size_t n) const;

  // Whether all elements of |b| satisfy |is_member|
  bool contains(const BitMap& b) const;

  // Whether none of the elements of |b| satisfy |is_member|
  bool disjoint(const BitMap& b) const;

  BitMap operator~ () const { return BitMap(*this).take_complement(); }

  BitMap operator& (const BitMap& other) && { return operator&=(other); }

  BitMap operator| (const BitMap& other) && { return operator|=(other); }

  BitMap operator^ (const BitMap& other) && { return operator^=(other); }

  BitMap operator& (BitMap&& other) const & { return other &= *this; }

  BitMap operator| (BitMap&& other) const & { return other &= *this; }

  BitMap operator^ (BitMap&& other) const & { return other ^= *this; }

  BitMap operator& (const BitMap& other) const &
  { BitMap result(*this); result&=other; return result; }

  BitMap operator| (const BitMap& other) const &
  { BitMap result(*this); return result|=other; }

  BitMap operator^ (const BitMap& other) const &
  { BitMap result(*this); return result^=other; }



// manipulators

  // this was called |resize|, but sets |capacity()|, whence the new name
  void set_capacity(size_t n); // any new bits will start out cleared
  void extend_capacity(bool b); // extend capacity by |1|, adding member if |b|

  void fill(); // set all bits
  void reset() { d_map.assign(d_map.size(),0ul); } // clear all bits
  void fill(size_t start,size_t stop); // set consecutive range of bits
  void clear(size_t start,size_t stop); // clear consecutive range of bits
  /*
    Set the bit at position |n| (that is, inserts the value |n| into the set);
    this makes |is_member(n)| hold.
  */
  void insert(size_t n)&
  {
    assert(n<d_capacity);
    d_map[n >> index_shift] |= constants::eq_mask[n & position_bits];
  }

  BitMap&& insert(size_t n)&&
  {
    assert(n<d_capacity);
    d_map[n >> index_shift] |= constants::eq_mask[n & position_bits];
    return std::move(*this);
  }

  // insert a range of values (which need not be listed increasingly)
  template<typename I> void insert(I, I);

  // Clear the bit at position |n|, making |is_member(n)| false.
  void remove(size_t n)&
  {
    assert(n<d_capacity);
    d_map[n >> index_shift] &= ~constants::eq_mask[n & position_bits];
  }

  BitMap&& remove(size_t n)&&
  {
    assert(n<d_capacity);
    d_map[n >> index_shift] &= ~constants::eq_mask[n & position_bits];
    return std::move(*this);
  }

  void flip(size_t n)&
  {
    assert(n<d_capacity);
    d_map[n >> index_shift] ^= constants::eq_mask[n & position_bits];
  }

  BitMap&& flip(size_t n)&&
  {
    assert(n<d_capacity);
    d_map[n >> index_shift] ^= constants::eq_mask[n & position_bits];
    return std::move(*this);
  }

  bool set_to(size_t n,bool b)
  { if (b) insert(n); else remove(n); return b; }

  bool set_mod2(size_t n, unsigned v) { return set_to(n,(v&1)!=0); }



  BitMap& take_complement ();

  // these operators allow argument to have different capacity than |*this|
  // in that case the capacity may change as indicated

  BitMap& operator&= (const BitMap& b); // change capacity to |b|'s if smaller

  BitMap& operator|= (const BitMap& b); // change capacity to |b|'s if greater

  BitMap& operator^= (const BitMap& b); // change capacity to |b|'s is smaller

  BitMap& andnot(const BitMap& b); // remove bits of |b|, keep capacity

  // the following two adjust the capacity by |-delta| respectively |+delta|
  BitMap& operator>>= (size_t delta); // shift bits right (decrease values)
  BitMap& operator<<= (size_t delta); // shift bits left (increase values)

// low level chunk access operations for speed

  // set an interval of bits from those (least significant ones) of source
  void set_range(size_t start, unsigned amount, chunk_type source);

  // get a range of bits as |chunk_type| value; see bitmap.ccp for limitations
  chunk_type range(size_t first, unsigned number) const;

// iteration related methods

  // First value for which |is_member| holds, or |capacity()| if none exist
  size_t first() const;

  // Last value for which |is_member| holds, or |capacity()| if none exist
  size_t last() const;

  // decrement |n| (at least once) until it points to a member; whether success
  bool back_up(size_t& n) const;

  /* Iterator to traverse the \emph{set} bits (members present) of a BitMap.

  Because of the nature of a bitmap, only constant iterators make sense (just
  like for iterators into std::set); one cannot "change the value" of an
  element at a given position because the position is determined by the value
  (values are always increasing; in fact a small value-change could be
  realised by swapping a set bit with neighboring unset bits, but that is of
  course not what a non-constant iterator should allow doing). However, these
  iterators do allow changing the underlying bitmap during traversal, even
  though such changes cannot be performed using only the iterator itself (this
  is unlike iterators over a |bitset::BitSet|, which copy the set of bits into
  their own value and therefore will traverse the bits that were set at the
  time of their construction, usually by the |bitset::BitSet::begin| method).

  The most delicate operation here is the increment, which has to find the
  position of the next set bit, while avoiding falling off the bitmap if there
  is no such. Therefore we included the data for the end of the bitmap in the
  iterator. This also allows the |operator()| internal test for exhaustion.
  */
  class const_iterator
    : public std::iterator<std::forward_iterator_tag, size_t>
  {
    chunk_const_iterator chunk_it; // point to current chunk
    size_t bit_index; // value returned if dereferenced
    size_t at_end; // beyond-the-end bit-address, normally constant

   public:
  // constructors and destructors

    const_iterator
      (const chunk_const_iterator& p, size_t n, size_t c) // from raw data
      : chunk_it(p)
      , bit_index(n)
      , at_end(c)
    {}

  // accessors
    bool operator== (const const_iterator& i) const
      { return bit_index == i.bit_index; } //  |chunk_it| is ignored!
    bool operator!= (const const_iterator& i) const
      { return bit_index != i.bit_index; } //  |chunk_it| is ignored!
    bool operator() () const { return bit_index != at_end; }

    const value_type& operator* () const { return bit_index; }

  // manipulators
    const_iterator& operator++ ();
    const_iterator operator++ (int);

    const_iterator& operator-- (); // provided, but less useful than |++|

  protected:
    chunk_type& target_chunk () const
    { return const_cast<chunk_type&>(*chunk_it); }

  }; // |class const_iterator|

  const_iterator begin() const { return iterator_at(first()); }
  const_iterator end() const
  { return const_iterator(d_map.end(),capacity(),capacity()); }

  const_iterator iterator_at(size_t n) const
  { return const_iterator(d_map.begin()+(n>>index_shift),n,capacity()); }

  class iterator : public const_iterator
  {
  public:
    using const_iterator::const_iterator; // use the same constuctor(s)
  // accessors
    bool operator== (const iterator& i) const
    { return const_iterator::operator==(i); }
    bool operator!= (const iterator& i) const
    { return const_iterator::operator!=(i); }

  // manipulators
    iterator& operator++ () { const_iterator::operator++(); return *this; }
    iterator operator++ (int) { auto tmp=*this; ++(*this); return tmp; }

    iterator& operator-- () { const_iterator::operator--(); return *this; }

    void remove_target() const
    { target_chunk() &=	~constants::eq_mask[*(*this)&position_bits]; }
  }; // |class iterator|

  iterator begin() { return iterator_at(first()); }
  iterator end() { return iterator(d_map.end(),capacity(),capacity()); }

  iterator iterator_at(size_t n)
  { return iterator(d_map.begin()+(n>>index_shift),n,capacity()); }

}; // |class BitMap|


} // |namespace bitmap|

#endif

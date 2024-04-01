/*
  This is bitmap.cpp

  Copyright (C) 2004,2005 Fokko du Cloux
  Copyright (C) 2006,2017,2022 Marc van Leeuwen
  part of Coxeter and of the Atlas of Lie Groups and Representations

  For license information see the LICENSE file
*/

#include "bitmap.h"

#include <algorithm> // for |lower_bound|, |copy|, |copy_backward|
#include <iterator> // for |std::distance|
#include <stdexcept>

#include "constants.h" // for |first_bit|, |last_bit|
#include "bits.h" // for |bitCount|

#include <cassert>

namespace bitmap {

/*****************************************************************************

  The implementation of the BitMap class.

  A BitMap should be seen as a container of |size_t| values. The fact that bits
  are used to signal the presence or absence of numbers is an implementation
  details that this class somewhat hides (though the class name is an obvious
  giveaway). It obeys the semantics of a Forward Container (from the C++
  standard library). The implementation is close to what would be expected the
  special class instance |std::vector<bool>| to be (namely using just one bit
  per stored boolean, which the standard makes special previsions for to allow),
  but the interface it presents is more like the one of |std::set<size_t>|.

  A bitmap is a implemented as a vector of |chunk_size|, each representing a
  "chunk" of bits in the map; we have chosen this unsiged type to be as large as
  possible, which improves speed in searching withing fairly sparse sets, so
  that it is likely to be equal |size_t|, but we maintain the name distinction
  ot not conflate two different roles. We do not wish to provide bit-address
  access to this map as |std::vector<bool>| must in some form. By contrast we do
  define an iterator class, which traverses the _set_ bits of the bitmap; so
  that dereferencing for instance |b.begin()| returns the index of the first set
  bit in the bitmap, and loops can use the |for (auto k : bitmap)| syntax.

******************************************************************************/



/*****************************************************************************

        Chapter I -- The BitMap class

******************************************************************************/


/******** constructors and destructors ***************************************/

/*
  In this constructor template we assume that I and J are iterator types with
  the same value_type. The idea is that [first,last[ is an ordered range, for
  which we can call |std::lower_bound| to do binary search. We construct
  the bitmap which flags the elements from [fsub,lsub[ (not necessarily assumed
  ordered or in range) within that range (renumbered consecutively from 0).
  Therefore |I| should be a RandomAccessIterator type, and binary search
  will be (repeatedly) applied to the range [first,last[; by contrast |J| can be
  any InputIterator. Any values from [fsub,lsub[ not found in [first,last[ will
  not be represented in the constructed |BitMap| at all.
*/
template <typename I, typename J>
  BitMap::BitMap(const I& first, const I& last, const J& fsub, const J& lsub)
    : d_capacity(std::distance(first,last))
    , d_map((d_capacity+position_bits)>>index_shift,0)
{
  for (J j = fsub; j != lsub; ++j)
  { auto it = lower_bound(first,last,*j);
    if (*it != last)
      insert(std::distance(first,*it));
  }
}



/******** accessors **********************************************************/

/*
  Return the number of set bits in the bitmap (this is its size as a
  container of size_t.)

  NOTE: correctness depends on unused bits in the final word being cleared.
*/
size_t BitMap::size() const
{
  size_t count = 0;

  for (chunk_type chunk : d_map)
    count += bits::bitCount(chunk);

  return count;
}

/*
  Whether the bitmap is empty. Thanks to our convention of zeroing
  unused bits, it is enough to check whether all the components of |d_map|
  are zero.
*/
bool BitMap::empty() const
{
  for (auto chunk : d_map)
    if (chunk!=0)
      return false;

  return true;
}

/*
  Tell whether the bitmap is full. This means that all blocks are full,
  except maybe for the last one where we have to look only at the significant
  part.
*/
bool BitMap::full() const
{
  for (const auto& chunk : d_map)
    if (~chunk!=0)
    {
      if (&chunk != &d_map.back() or (capacity()&position_bits)==0)
	return false;
      else return chunk==constants::lt_mask[capacity()&position_bits];
    }

  return true;
}


/*
  Return the index (value as container element) of set bit number |ordinal|. In
  other words, if a |std:set<size_t>| were used as implementation instead, this
  would be |*std::next(begin(),ordinal|; the syntax |b[i]| as alternative for
  our |b.n_th(i)| would have been logical, but we want to avoid the suggestion
  that this is a cheap operation. The method returns |d_capacity| if there is no
  such element (if |ordinal>=size()|) but it coyuld have thrown an error.
  Whenever |0<=i<=size()| one has |position(n_th(i))==i|.
*/
size_t BitMap::n_th(size_t ordinal) const
{
  size_t pos = 0;

  for (chunk_type chunk : d_map)
  {
    size_t chunk_count = bits::bitCount(chunk);
    if (ordinal>=chunk_count)
    {
      ordinal -= chunk_count;
      pos += chunk_bits;
    }
    else // |ordinal<chunk_count|
    { // now the required bit is going to be in |chunk|
      while (ordinal-->0) // clear |ordinal| initial bits in |f|
	chunk &= (chunk-1); // clear lowest set bit
      return pos + constants::first_bit(chunk); // the bit index of our set bit
    }
  }

  // if we reach this point the bit is not found
  return capacity(); // indicate out-of-bounds value
}

/*
  Return the number of set bits in positions |< n|; viewing a bitset |b| as a
  container of |size_t|, this is the number of values |< n| that |b| contains.
  If |is_member(n)| holds, then one has |n==n_th(position(n))|.
*/
size_t BitMap::position(size_t n) const
{
  size_t p = 0;
  size_t b = n >> index_shift;

  for (size_t j = 0; j < b; ++j)
    p += bits::bitCount(d_map[j]);

  if ((n&position_bits)!=0)
    p += bits::bitCount(d_map[b] & constants::lt_mask[n&position_bits]);

  return p;
}

/*
  Whether the current bitmap contains |b|.

  This would amount to |b.andnot(*this).empty()| if |b| were by-value rather
  than by contant reference. Meanwhile, the code below is more efficient too.
*/
bool BitMap::contains(const BitMap& b) const
{
  size_t limit = b.d_map.size();
  if (limit>d_map.size()) // then test for above-our-capacity entries
  {
    size_t n = b.last();
    if (n==b.capacity())
      return true; // since |b.empty()|
    if (n>=capacity())
      return false; // since |not is_member(n)|
    limit = (n >> index_shift)+1;
  }

  for (size_t j = 0; j < limit; ++j)
    if ((b.d_map[j] & ~d_map[j])!=0)
      return false; // since we found a bit set in |b| but not in |*this|

  return true; // since we checked all set bits of |b|
}

// Whether the current bitmap is disjoint from |b|.
bool BitMap::disjoint(const BitMap& b) const
{
  size_t m = std::min(d_map.size(),b.d_map.size());
  for (size_t j = 0; j < m; ++j)
    if ((d_map[j] & b.d_map[j])!=0)
      return false;

  return true;
}

bool BitMap::equivalent (const BitMap& b) const
{
  size_t m = std::min(d_map.size(),b.d_map.size());
  for (size_t j = 0; j < m; ++j)
    if (d_map[j] != b.d_map[j])
      return false;
  return m==d_map.size()
    ? std::all_of(&b.d_map[m],&*b.d_map.end(),[](chunk_type c) { return c==0; })
    : std::all_of(&d_map[m],&*d_map.end(),[](chunk_type c) { return c==0; });
}

/******** manipulators *******************************************************/

/*
  Set the capacity of the bitmap to |n|, shrinking |d_map| if possible

  Does not modify the contents up to the previous size, at least if n is
  larger. The new elements are initialized to zero.
*/
void BitMap::set_capacity(size_t n)
{
  if (n<d_capacity and (n&position_bits)!=0)
    d_map[n>>index_shift] &= constants::lt_mask[n&position_bits]; // 0 partial
  d_map.resize((n+position_bits) >> index_shift,0); // grow zeroing or shrink
  if (d_map.size()<d_map.capacity()) // |d_map| has become too large
    d_map.shrink_to_fit(); // shrink-wrap vector |d_map|
  d_capacity = n;
}

// Add one more place to bitmap (amortised efficiently), set it to value |b|
void BitMap::extend_capacity(bool b)
{
  if ((d_capacity&position_bits) == 0)
    d_map.push_back(0); // allocate space; |std::vector| handles efficiency
  set_to(d_capacity++,b);
}

/*
  Set all the bits in the bitmap.

  As usual we have to be careful to leave the unused bits at the end to zero.
*/
void BitMap::fill()
{
  d_map.assign(d_map.size(),~0ul); // (last word may be overwritten next)

  size_t remainder = d_capacity&position_bits;
  if (remainder!=0) // N.B. also makes |capacity()=0| safe (MvL)
    d_map.back() = constants::lt_mask[remainder];
}

// Set all the bits in positions |i| with |start<=i<stop|.
void BitMap::fill(size_t start, size_t stop)
{
  assert(start<=capacity() and stop<=capacity()); // do not try to grow us
  if (start>=stop) return;
  size_t begin = start >> index_shift;    // chunk index of first bit to set
  size_t end   = (stop-1) >> index_shift; // chunk index of last bit to set
  if (begin==end) // then only one word is affected
    d_map[begin] |= constants::lt_mask[stop-start] << (start&position_bits);
  else
  {
    d_map[begin] |= ~constants::lt_mask[start&position_bits];
    for (size_t i=begin+1; i<end; ++i)
      d_map[i] = ~0ul;
    d_map[end] |= constants::lt_mask[stop&position_bits];
  }
}

// Set all the bits in positions |i| with |start<=i<stop|.
void BitMap::clear(size_t start, size_t stop)
{
  if (start>=stop) return;
  size_t begin = start >> index_shift;
  size_t end   = stop >> index_shift;
  if (begin==end) // then only one word is affected
    d_map[begin] &= ~(constants::lt_mask[stop-start] << (start&position_bits));
  else
  {
    d_map[begin] &= constants::lt_mask[start&position_bits];
    for (size_t i=begin+1; i<end; ++i)
      d_map[i] = 0;
    if ((stop&position_bits)!=0) // protect against out-of-bounds setting of 0 bits
      d_map[end] &= ~constants::lt_mask[stop&position_bits];
  }
}

/*
  Here we assume that I is an iterator whose value_type is convertible to
  |size_t|, and we do the sequence of insertions from the range [first,last[.
*/
template<typename I> void BitMap::insert(I first, I last)
{
  while (first!=last)
  {
    insert(*first);
    ++first;
  }
}


/*
  Replace a bitmap by its complement

  We are careful about the last chunk, leaving the unused bits equal to zero.
*/
BitMap& BitMap::take_complement ()
{
  for (chunk_type& chunk : d_map)
    chunk = ~chunk;

  if ((d_capacity&position_bits)!=0) // also makes |capacity()==0| safe here
    d_map.back() &= constants::lt_mask[d_capacity&position_bits];

  return *this;
}

// Intersect the current bitmap with |b|, return whether result is non-empty.
BitMap& BitMap::operator&= (const BitMap& b)
{
  if (b.capacity()<capacity())
  {
    d_capacity=b.capacity();
    d_map.erase(d_map.begin()+b.d_map.size(),d_map.end());
  }
  for (size_t i = 0; i < d_map.size(); ++i)
    d_map[i] &= b.d_map[i];

  return *this;
}

/*
  Take the current bitmap into its set-difference with |b|, i.e.,
  remove from our bitmap any elements appearing in |b|.
*/
BitMap& BitMap::andnot(const BitMap& b)
{
  size_t limit = std::min(d_map.size(),b.d_map.size());
  for (size_t i = 0; i < limit; ++i)
    d_map[i] &= ~b.d_map[i];

  return *this;
}

// Unite |b| into the current bitmap (which may thereby grow its capacity)
BitMap& BitMap::operator|= (const BitMap& b)
{
  size_t limit = b.d_map.size();
  if (b.capacity()>capacity())
  {
    limit = d_map.size();
    d_capacity=b.capacity();
    d_map.insert(d_map.end(),&b.d_map[limit],&*b.d_map.end());
  }
  for (size_t i = 0; i < limit; ++i)
    d_map[i] |= b.d_map[i];

  return *this;
}

// Perform exclusive or (XOR) of |b| into the current bitmap.
BitMap& BitMap::operator^= (const BitMap& b)
{
  size_t limit = b.d_map.size();
  if (b.capacity()>capacity())
  {
    limit = d_map.size();
    d_capacity=b.capacity();
    d_map.insert(d_map.end(),&b.d_map[limit],&*b.d_map.end());
  }
  for (size_t i = 0; i < limit; ++i)
    d_map[i] ^= b.d_map[i];

  return *this;
}


BitMap& BitMap::operator<<= (size_t delta) // increase position values by |delta|
{
  d_capacity += delta;
  d_map.resize((d_capacity+position_bits)>>index_shift,0);

  size_t delta_rem = delta & position_bits;
  delta >>= index_shift; // must move by |delta| words, and then |delta_rem| bits
  chunk_iterator start=d_map.begin(); // cannot be const_iterator; guess why
  if (delta>0) // useless shifting |0| violates |std::copy_backward| requirements
    start = std::copy_backward(d_map.begin(),d_map.end()-delta,d_map.end());
  if (delta_rem>0 and not d_map.empty())
  {
    chunk_iterator it=d_map.end();
    while (--it!=start) // reverse loop, omit initial
    {
      *it <<= delta_rem; // shift bits up
      *it |= *(it-1) >> (chunk_bits-delta_rem); // bits shifted in
    }
    assert(it==start);
    *it <<= delta_rem; // shift bits up in lowest copied word
  }
  std::fill(d_map.begin(),start,0ul);
  return *this;
}

BitMap& BitMap::operator>>= (size_t delta) // decrease values by |delta|
{
  if (delta>=capacity())
    return *this = BitMap(); // clear |d_capacity| and |d_map|
  assert(not d_map.empty()); // since |capacity()>0| implies this
  const size_t old_rem = d_capacity & position_bits;
  d_capacity -= delta;
  size_t delta_rem = delta & position_bits;
  delta >>= index_shift; // must move by |delta| words, and then |delta_rem| bits
  chunk_const_iterator stop=d_map.end();
  if (delta>0) // useless shifting by |0| violates |std::copy| requirements
    stop = std::copy(d_map.begin()+delta,d_map.end(),d_map.begin());
  if (delta_rem>0)
  {
    --stop; // we need to stop
    chunk_iterator it;
    for (it=d_map.begin(); it!=stop; ++it)
    {
      *it >>= delta_rem; // shift bits down (to the right) inside the word
      *it |= *(it+1) << (chunk_bits-delta_rem); // left part from next
    }
    assert(it==stop);
    if (old_rem>delta_rem)
    {
      *it >>= delta_rem; // shift last bits down internally
      ++stop; // not all infomation was shifted out of |*stop|, keep it
    } // else we've already moved all its bits down, do drop |*tstop|
  }
  d_map.erase(stop,d_map.end());
  return *this;
}



/*
  Return |r| bits from position |n| as a |chunk_type| value

  Precondition: |(n+r-1)%chunk_bits == n%chunk_bits + (r-1)|, in other words
  no multiples $m$ of |chunk_bits| satisfy $n<m<n+r$. This is
  ensured for instance when |r| divides both |chunk_bits| and |n|.

  Thus the bits extracted are found in single element of |d_map|.

  It is required that |n<capacity()|, but not that |n+r<=capacity()|; if the
  latter fails, the return value is padded out with (leading) zero bits.
*/
auto BitMap::range(size_t n, unsigned r) const -> chunk_type
{
  size_t m = n >> index_shift; // index where data is stored

  return (d_map[m] >> (n & position_bits)) & constants::lt_mask[r];

}

/*
  Set |r| bits from position |n| to the first |r| bits in |a|.

  Precondition: |n| and |n+r-1| have same integer quotient by |longBits| (or
  |r==0| in which case nothing happens), so in particular |r<=longBits|.
  This condition is always satisfied if |r| divides |longBits| and |n| is a
  multiple of |r|. It ensures that only a single word in |d_map| is affected.
*/
void BitMap::set_range(size_t n, unsigned r, chunk_type a)
{
  size_t m = n >> index_shift;
  unsigned shift = n & position_bits;

  d_map[m] &= ~(constants::lt_mask[r] << shift);     // set bits to zero
  d_map[m] |= (a & constants::lt_mask[r]) << shift;  // insert the bits from a
}


/*****************************************************************************

        Chapter II -- The BitMap::iterator class

  Because of the nature of a bitmap, only constant iterators make sense (just
  like for iterators into |std::set|); one cannot "change the value" of an
  element at a given position because the position is determined by the value
  (values are always increasing; in fact a small value-change could be realised
  by swapping a set bit with neighboring unset bits, but that is of course not
  what a non-constant iterator should allow doing). On the other hand we may
  modify our bitmap |M| using this iterator |it|, notably to clear the set bit
  |it| points at by |it.remove_target()|, which has the same effect as
  |M.remove(*it)|. Unlike the situation for |std::set|, such a removal does not
  invalidate the iterator |it| itself, so it is not an invariant of the
  |BitMap::iterator| class that it always points at a set bit.

  It is well defined what happens when one inserts new elements into a |BitMap|
  while iterating with an iterator |it|: if the value inserted is greater than
  |*it|, then further incrementing |it| will at some point stop at the inserted
  value. If the value inserted is less than |*it|, then further incrementing of
  |it| will no be affected by the insertion. This behaviour is because no part
  of the bitmap itself is copied into the iterator: for each increment we fetch
  the word containing the bit we were on, but shift away any bits up to and
  including that bit before advancing to the next set bit, if any. It looping
  until |it==end()|, that condition is unaffected by the state of the bits.

  The most delicate operation is |operator++|, which has to move to the the
  next set bit, or stop at the end of the bitmap if there is no such.
  Therefore we included the data for the end of the bitmap in the iterator.

******************************************************************************/

/*
  The address of the first member (set bit) of the bitmap,
  or past-the-end indicator |d_capacity| if there is no such.
*/
size_t BitMap::first() const
{
  for (const chunk_type& chunk : d_map)
    if (chunk!=0)
    {
      auto offset = &chunk - &d_map[0];
      return (offset<<index_shift) + constants::first_bit(chunk);
    }

  return capacity(); // signal that we are empty
} // |first|


// We don't do reverse iterators, however reverse interation is easy

// initial value in reverse iteration, |size_t n=capacity(); return back_up(n),n|
size_t BitMap::last () const
{ // we assume any beyond-capacity bits are cleared
  size_t i = d_map.size();
  while (i-->0)
    if (chunk_type m=d_map[i]) // whether chunk is nonzero
      return (i<<index_shift)+constants::last_bit(m);

  return capacity(); // signal that we are empty
} // |last|

/*
  Decrement |n| until it points to a member of the bitset,
  or if none is found returns |false| (in which case |n| is unchanged)

  One can use |size_t n=b.capacity()| followed by |while(b.back_up(n)) {...}|
  or, when it is is known that |not b.empty()|, (marginally faster):
  |size_t n=b.last()| followed by |do {...} while(b.back_up(n));|
*/
bool BitMap::back_up(size_t& n) const
{
  auto i=n>>index_shift; // the index in |d_map| containing bit |n|
  chunk_type m =
    (n&position_bits)==0 ? 0 // taken aside to avoid dereferencing |d_map[i]|
    : d_map[i]&constants::lt_mask[n&position_bits]; // masks out bit |n| and larger

  while(m==0 and i>0) // search backwards for a nonzero chunk
    m=d_map[--i];

  if (m==0) return false;
  n = (i<<index_shift)+constants::last_bit(m);

  return true;
} // |back_up|


/*
  The incrementation operator; it has to move the bitAddress to the next
  set bit, and move the chunk if necessary.

  This code below assumes that in case of an incomplete last chunk, there are
  no bits set in that chunk beyond the end of the bitmap; if there were,
  |first_bit(f)| below (both instances) could make the iterator advance to such
  a bit when it should have halted at |d_capacity|.
*/

auto BitMap::const_iterator::operator++ () -> const_iterator&
{
  const auto cur_pos = bit_index & position_bits;
  chunk_type f = // current chunk masked to bits after current bit
    *chunk_it & ~constants::leq_mask[cur_pos];

  if (f!=0) { // if there is still some bit set in this chunk, jump to it
    bit_index += constants::first_bit(f)-cur_pos;
    return *this;
  }

  // if not, we're done with this chunk; we'll advance |chunk_it| at least once
  bit_index &= chunk_index_bits; // prepare to advance by whole chunks

  const auto old_chunk = chunk_it;
  const auto limit = chunk_it + // offset to beyond-the-end of |d_map|
    ((at_end+position_bits-bit_index) >> index_shift);

  do
    if (++chunk_it==limit) // advance, test
    { // we've run out of chunks, the value of |chunk_it| no longer matters
      bit_index = at_end; // so set beyond-the-end indication
      return *this;
    }
  while ((f=*chunk_it)==0); // pick up next chunck, repeat if it is entirely 0

  // now the bit to advance to is in the current chunk, and not out of bounds
  bit_index += ((chunk_it-old_chunk) << index_shift) + constants::first_bit(f);
  return *this;
}


/*
  Post-increment operator; it should return the value as it was _before_ the
  incrementation. This operator can mostly by avoided, as |M.remove(*it++)|
  can safely be replaced by |M.remove(*it),++it|
*/
auto BitMap::const_iterator::operator++ (int) -> const_iterator
{
  auto tmp = *this;
  ++(*this);
  return tmp;
}

auto BitMap::const_iterator::operator-- () -> const_iterator&
{
  unsigned rem = bit_index & position_bits;
  bit_index &= chunk_index_bits;
  auto m = rem==0 ? 0 // avoid dereferencing |chunk_it|
    : *chunk_it & constants::lt_mask[rem]; // masks out bit |rem| and larger

  if (m==0) // if not we can spare out some work
  {
    const auto limit = chunk_it - (bit_index >> index_shift);
    while (m==0 and chunk_it>limit)
      m=*--chunk_it;

    if (m==0)
      throw std::runtime_error("BitMap iterator underflow");
    bit_index = (chunk_it - limit) << index_shift;
  }
  bit_index += constants::last_bit(m);
  return *this;
} // |operator--|


// Instantiations


} // |namespace bitmap|

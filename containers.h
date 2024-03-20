/*
  This is containers.h, container type customization

  Copyright (C) 2016,2018,2024 Marc van Leeuwen
  developed as part of Atlas and later of Coxeter

  For license information see the LICENSE file
*/

#ifndef CONTAINERS_H /* guard against multiple inclusions */
#define CONTAINERS_H

#include <cstdint> // for |std::uint32_t|
#include <vector> // for |std::vector|
#include <set>    // for |std::set|
#include <map>    // for |std::map|
#include <algorithm> // for |std::copy|, |std::all|
#include <cassert>

#include "memory.h" // for |containers::allocator|

namespace containers {

template <typename T>
  using vector = std::vector<T,containers::allocator<T> >;

template<typename T, typename Compare = std::less<T> >
  using set = std::set<T,Compare,containers::allocator<T> >;

template<typename Key, typename T, typename Compare = std::less<Key> >
  using map = std::map<Key,T,Compare,containers::allocator<T> >;

template<typename Key, typename T, typename Compare = std::less<Key> >
  using multimap =
    std::multimap<Key,T,Compare,
		  containers::allocator<std::pair<const Key,T> > >;

  template <typename Alloc>
  class allocator_deleter;
template<typename T,typename Alloc = containers::allocator<T> >
  struct sl_node;

template<typename T, typename Alloc = containers::allocator<T> >
  class sl_list_const_iterator;
template<typename T,typename Alloc = containers::allocator<T> >
  class sl_list_iterator;

template<typename T, typename Alloc = containers::allocator<T> >
  class weak_sl_list_const_iterator;
template<typename T, typename Alloc> class
  weak_sl_list_iterator;

template<typename T,typename Alloc = containers::allocator<T> >
  class simple_list;
template<typename T,typename Alloc = containers::allocator<T> >
  class sl_list;

template<typename T,typename Alloc = containers::allocator<T> >
  class mirrored_simple_list; // trivial adapter, to allow use with |std::stack|

template<typename T,typename Alloc = containers::allocator<T> >
  class mirrored_sl_list; // trivial adapter, to allow use with |std::stack|

template<typename T,typename Alloc = containers::allocator<T> > class stack;

template<typename T,typename Alloc = containers::allocator<T> > class queue;

// the following is a 2-dimensional counterpart to |vector|
// it was obtained as a simplification of the Atlas type |Matrix_base|
template<typename C> class matrix
{
  using index_t = std::uint32_t; // limit to $2^{32}$ rows and columns
 // dimensions
  index_t n_rows; // one could make this a |size_t| if need arises
  index_t n_cols;
  vector<C> body;   // vector of elements, concatenated by rows

public:

// constructors
  matrix(): n_rows(0),n_cols(0),body() {}

  matrix(index_t m, index_t n) // leaves all entries unset
    : n_rows(m),n_cols(n),body(static_cast<std::size_t>(m)*n) {}
  matrix(index_t m, index_t n, const C& c)
    : n_rows(m), n_cols(n), body(static_cast<std::size_t>(m)*n,c) {}

#if 0 // wait with implementation until need arises
  matrix (const vector<vector<C> >& b,
	  index_t n_cols); // with explicit #rows in case |b| empty

  template<typename I> // from sequence of columns obtained via iterator
  matrix(I begin, I end, index_t n_cols)
    : matrix(std::distance(begin,end),n_cols)
  {
    auto dest = body.begin();
    for (I p=begin; p!=end; ++p)
    {
      assert(std::distance(p->begin(),p->end())==n_cols);
      dest = std::copy(p->begin(),p->end(),dest);
    }
  }
#endif

  // accessors
  index_t nr_rows() const { return n_rows; }
  index_t nr_cols() const { return n_cols; }

  const C* row(index_t i) const
    { assert(i<n_rows); return &body[static_cast<size_t>(i)*n_cols]; }
  const C& entry(index_t i,index_t j) const { return row(i)[j]; }

  const C* at (size_t i,size_t j) const; // safe version of |entry|

  bool operator== (const matrix<C>&) const;
  bool operator!= (const matrix<C>& m) const { return not(operator==(m)); }

  bool is_zero() const
    { return std::all_of(body.begin(),body.end(),[](C c) { return c==C(0); }); }
  bool is_empty() const { return body.size() == 0; }
// manipulators

  void reserve(index_t n) { body.reserve(static_cast<size_t>(n)*n_cols); }
  void grow(index_t n, C c)
    { assert(n>=n_rows);
      body.insert(body.end(),static_cast<size_t>(n-n_rows)*n_cols,c);
      n_rows=n;
    }
  void shrink(index_t n)
    { assert(n<=n_rows);
      body.resize(static_cast<size_t>(n)*n_cols);
      n_rows=n;
    }

  C* row(index_t i)
    { assert(i<n_rows); return &body[static_cast<size_t>(i)*n_cols]; }
  C& entry(index_t i,index_t j) { return row(i)[j]; }

  C* at (size_t i,size_t j); // safe version of |entry|

  void append_row (const vector<C>& r)
    { assert(r.size()==n_cols); body.insert(body.end(),r.begin(),r.end()); }

  void swap_rows (index_t i0,index_t i1)
  { auto* p=row(i0), * q=row(i1);
    for (index_t j=0; j<n_cols; ++j) std::swap(p[j],q[j]);
  }
}; // |template<typename C> class matrix|



template<typename T, typename Compare = std::less<T> >
  struct bag : public set<T,Compare>
  {
    using Base = set<T,Compare>;
    using Base::Base; // inherit all constructors
    const T* find(const T& x)
    {
      auto p = Base::insert(x);
      return &*p.first; // convert iterator to pointer, ignore whether new
    }
  }; // |template<typename T> struct bag|

} // |namespace containers|

#endif

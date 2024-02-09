/*
  This is containers.h, container type customization

  Copyright (C) 2016,2018 Marc van Leeuwen
  developed as part of coxter

  For license information see the LICENSE file
*/

#ifndef CONTAINERS_H /* guard against multiple inclusions */
#define CONTAINERS_H

#include <vector> // for |std::vector|
#include "memory.h" // for |containers::allocator|

namespace containers {

template <typename T>
  using vector = std::vector<T,containers::allocator<T> >;



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

} // |namespace containers|

#endif

/*
  This is list.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef LIST_H  /* guarantee single inclusion */
#define LIST_H

#include "globals.h"
#include <stdio.h> // for |FILE|
#include <limits.h>
#include <algorithm>

/******** type declarations *************************************************/

namespace list {
  template <class T> class List;
};

/******** constants *********************************************************/

namespace list {
  const Ulong undef_size = ULONG_MAX;
  const Ulong not_found = ULONG_MAX;
}

/******** functions provided by list.h **************************************/

namespace list {
  // find first index of |m| in list, or return |list::not_found|
  template <class T> Ulong find(const List<T>& l, const T& m);

  // insert element |m| into ordered list |l|
  template <class T> Ulong insert(List<T>& l, const T& m);
  template <class T> void print(FILE* file, const List<T>& l);
};

/******** type definitions **************************************************/

#include "memory.h"

namespace list {

template <class T> class List {
 protected:
  T* d_ptr;
  Ulong d_size;
  Ulong d_allocated;

 public:
  typedef T eltType;
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(List<T>));}
  void* operator new(size_t, void* ptr) {return ptr;}
  void operator delete(void* ptr, void* placement) {};
  List() : d_ptr(nullptr), d_size(0), d_allocated(0) {}
  List(const Ulong& n);
  List(const List<T>& r);
  List(const T* p, const Ulong& n);
  template <class I> List(const I& first, const I& last);
  template <class I, class F> List(const I& first, const I& last, F& f);
  ~List();
/* modifiers */
  T& operator[] (Ulong j);                                        /* inlined */
  const List<T>& operator= (const List& r);
  void append(const T& x);
  const List<T>& assign(const List& r);
  void erase(const Ulong& n);
  void reverse();
  T* ptr() {return d_ptr;}
  void setData(const T* source, Ulong first, Ulong r);
  void setData(const T* source, Ulong r);                         /* inlined */
  void setSize(Ulong n); // set |d_size| while ensuring corresponding memory
  void setSizeValue(const Ulong& n);                              /* inlined */
  void setZero(Ulong first, Ulong r);                             /* inlined */
  void setZero(Ulong r);                                          /* inlined */
  void setZero() { setZero(0,d_size); } // zero out all entries, keep size
  void shallowCopy(const List& w);                                /* inlined */
  void shiftPtr(const long& d);                                   /* inlined */
  Ulong& size();                                                  /* inlined */
  void sort();
  template<class C> void sort(C& c);                              /* inlined */
/* accessors */
  const T& operator[] (Ulong j) const;                            /* inlined */
  bool operator== (const List& w) const;
  bool operator!= (const List& w) const;
  bool operator< (const List& w) const;
  const T* ptr() const;                                           /* inlined */
  const Ulong& size() const;                                      /* inlined */
/* iterator */
  typedef T* Iterator;
  typedef const T* ConstIterator;
  Iterator begin();                                               /* inlined */
  Iterator end();                                                 /* inlined */
  ConstIterator begin() const;                                    /* inlined */
  ConstIterator end() const;                                      /* inlined */
};

};

/******** Implementation of inline functions *******************************/

namespace list {

/* class List */

/* modifiers */

template<class T> inline T& List<T>::operator[] (Ulong j)
  {return d_ptr[j];}
template<class T>
inline void List<T>::setData(const T* source, Ulong r)
  {setData(source,0,r);}
template<class T> void List<T>::setSizeValue(const Ulong& n)
  {d_size = n;}
template<class T> inline void List<T>::setZero(Ulong first, Ulong r)
  { std::fill(d_ptr+first,d_ptr+first+r,static_cast<T>(0));}
template<class T> inline void List<T>::setZero(Ulong r) {setZero(0,r);}
template<class T> inline void List<T>::shallowCopy(const List<T>& w)
  { d_ptr = w.d_ptr; d_size=w.d_size; d_allocated=w.d_allocated; }
template<class T> inline void List<T>::shiftPtr(const long& d)
  { d_ptr += d; d_size -= d; d_allocated -= d;}
template<class T> Ulong& List<T>::size() {return d_size;}

/* accessors */

template <class T> inline const T& List<T>::operator[] (Ulong j) const
  {return(d_ptr[j]);}
template<class T> inline bool List<T>::operator!= (const List<T>& w) const
  {return !operator==(w);}
template<class T> const T* List<T>::ptr() const {return d_ptr;}
template <class T> inline const Ulong& List<T>::size() const {return d_size;}

/* iterators */

template <class T> inline typename List<T>::Iterator List<T>::begin()
  {return d_ptr;}
template <class T> inline typename List<T>::Iterator List<T>::end()
  {return d_ptr+d_size;}
template <class T> inline typename List<T>::ConstIterator List<T>::begin()
  const {return d_ptr;}
template <class T> inline typename List<T>::ConstIterator List<T>::end() const
  {return d_ptr+d_size;}

};

#include "list.hpp"

#endif

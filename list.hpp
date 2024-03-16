/*
  This is list.hpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include <new>
#include <algorithm> // for |std::copy|
#include <cassert>
#include <memory> // for |std::uninitialized_copy|

#include "error.h"

/**************************************************************************

  This file provides the definitions of the functions for the List and
  List classes. As far as I can see, List is the basic_string type
  of the standard library; the reason I didn't use the library class is
  mostly that I'm not familiar with it; also that I want to keep the code
  small, but most importantly, I want to keep memory management in my own
  hands. List is a variant of List where the sizes are fixed.

 **************************************************************************/

/**************************************************************************

        Chapter I -- The List class.

  This section provides the following functions :

  - constructors and destructors :

    - List(Ulong);
    - List(List&);
    - List(T*,ulong);
    - List(const I&, const I&);
    - List(const I&, const I&, F&);
    - ~List();

  - accessors :

    - operator== (const List&) : tests for equality;
    - operator< (const List&) : comparison operator;

  - modifiers :

    - append(const T&): appends an element (potentially resizing);
    - assign(const List&): copy constructor;
    - erase(const Ulong&): erases the n-th term in the list;
    - reverse(): reverses the order of the elements;
    - setSize(Ulong): resizes;
    - setData(T*,first,Ulong): sets the data, resizing if necessary;
    - sort() : sorts the list;

 **************************************************************************/

namespace list {


// Allocate |*this| to be able to hold |n| elements, but start out empty.
template <class T> List<T>::List(const Ulong& n)
  : d_ptr(static_cast<T*> (memory::arena().alloc(n*sizeof(T))))
  , d_size(0)
  , d_allocated(memory::arena().allocSize(n,sizeof(T)))
{}


template <class T> List<T>::List(const List<T>& r)

/*
  Copy constructor. Contrary to assignment, constructs fully new objects.
*/

{
  d_ptr = static_cast<T*> (memory::arena().alloc(r.size()*sizeof(T)));
  d_allocated = memory::arena().allocSize(r.size(),sizeof(T));
  for (Ulong j = 0; j < r.size(); ++j) {
    new(d_ptr+j) T(r[j]);
  }
  d_size = r.d_size;
}


template <class T> List<T>::List(const T* p, const Ulong& n)
  : List(n)
{
  std::uninitialized_copy(p,p+n,d_ptr);
}

/*
  A list constructor taking iterators as parameters. It is assumed that
  the value-type of I may be allocated to T.
*/
template <class T> template <class I>
List<T>::List(const I& first, const I& last)
  : List() // first construct as empty
{
  for (I i = first; i != last; ++i) {
    append(*i);
  }
}


/*
  Like the previous one, except that in addition F is a functor taking one
  argument of type I::value_type, and whose value-type may be allocated to T.
*/
template<class T> template<class I, class F>
List<T>::List(const I& first, const I& last, F& f)
  : List() // first construct as empty
{
  for (I i = first; i != last; ++i) {
    append(f(*i));
  }
}

template<class T> List<T>::~List()

/*
  Destructor for the List class : releases the memory.

  NOTE : this has essentially the effect of what delete[] d_ptr would do
  if I could get it to work.
*/

{
  for (Ulong j = 0; j < d_allocated; ++j) {
    d_ptr[j].~T();
  }
  memory::arena().free(d_ptr,d_allocated*sizeof(T));
}

/******** accessors *********************************************************/

template<class T> bool List<T>::operator== (const List<T>& w) const

/*
  Equality operator for lists. Two lists are equal if they have the same
  size, and if the elements in the range [0,size[ are pairwise equal. It
  assumes that operator== is defined for T.
*/

{
  if (d_size != w.size())
    return false;

  for (Ulong j = 0; j < d_size; ++j) {
    if (!(d_ptr[j] == w[j]))
      return false;
  }

  return true;
}

template<class T> bool List<T>::operator< (const List<T>& w) const

/*
  Comparison operator for lists. Comparison is length-first, lexicographical.
  It assumes operator< is defined for T.
*/

{
  if (d_size < w.size())
    return true;
  if (d_size > w.size())
    return false;

  /* if we reach this point, sizes are equal */

  for (Ulong j = 0; j < d_size; ++j) {
    if (d_ptr[j] < w[j])
      return true;
    if (d_ptr[j] > w[j])
      return false;
  }

  /* if we reach this point, lists are equal */

  return false;
}

/******** modifiers *********************************************************/

template <class T>
const List<T>& List<T>::operator= (const List<T>& r)

{
  assign(r);
  return *this;
}

// Append one element to the list, resizing if necessary.
/*
  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.
*/
template <class T> void List<T>::append(const T& x)
{
  // we have to be careful in case |x| points into the structure being
  // resized! calling |setSize| directly could invalidate |x|.

  Ulong c = d_size;

  if (d_allocated < c+1)
  {
    T* new_ptr = static_cast<T*> (memory::arena().alloc((c+1)*sizeof(T)));
    if (new_ptr==nullptr) /* overflow */
      { assert(error::ERRNO!=0); return; }
    std::copy(d_ptr,d_ptr+c,new_ptr); // copy the whole old range of values
    new_ptr[c] = x;
    memory::arena().free(d_ptr,d_allocated*sizeof(T));
    d_ptr = new_ptr;
    d_allocated = memory::arena().allocSize(c+1,sizeof(T));
    d_size = c+1;
    return;
  }

  // if we get here no resizing is necesary

  setSize(c+1);
  d_ptr[c] = x;

  return;
}



/*
  Assign r to the current list, by a one-level copy (in other words,
  treats Lists as BasicStrings).

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_OVERFLOW is set.

  NOTE : the return value in case of error should rather be 0, but this
  would mean that the function should really return a pointer.
*/

template <class T> const List<T>& List<T>::assign(const List<T>& r)
{
  setSize(r.size());
  if (error::ERRNO) /* overflow */
    return *this;
  setData(r.ptr(),r.size());
  return *this;
}

template <class T> void List<T>::erase(const Ulong& n)

/*
  This function erases the n-th term in the list, by shifting up.
*/

{
  memmove(d_ptr+n,d_ptr+n+1,(d_size-n-1)*sizeof(T));
  d_size--;
}

template <class T> void List<T>::reverse()

/*
  Reverses the order of the elements in the list.
*/

{
  for (Ulong j = 0; j < d_size/2; ++j) {
    T a = d_ptr[j];
    d_ptr[j] = d_ptr[d_size-j-1];
    d_ptr[d_size-j-1] = a;
  }

  return;
}

// Check if |*this| will hold |n| nodes of data, and resize it if not.

// Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
template <class T> void List<T>::setSize(Ulong n)
{
  if (d_allocated < n) { /* resize */
    void *p = memory::arena().alloc(n*sizeof(T));
    if (error::ERRNO) /* overflow */
      return;
    auto* pp = static_cast<T*> (p);
    // the next line could use |std::uninitialized_move| from C++17
    std::uninitialized_copy(d_ptr,d_ptr+d_size,pp);
    memory::arena().free(d_ptr,d_allocated*sizeof(T)); // clean up emptied memory
    d_ptr = pp; // henceforth use the new, partially initialized, memory
    d_allocated =
      memory::arena().allocSize(n,sizeof(T)); // use actual space bought
  }

  d_size = n;

  return;
}


/*
  After resizing if necessary, move the first |r| entries of |source| to the
  list, after its |first| elements (and discarding the remaining ones).

  Although |source| could point into our current List, this function should
  not be called in such a manner that it would cause some entries to be
  duplicated, namely with |d_ptr<=source<d_ptr+first|

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set.
*/
template <class T>
void List<T>::setData(const T *source, Ulong first, Ulong r)
{
  // we have to be careful in case |source| points into the structure being
  // resized! calling |setSize| directly would invalidate source.

  if (d_allocated < first+r)
  {
    T* new_ptr = static_cast<T*> (memory::arena().alloc((first+r)*sizeof(T)));
    if (new_ptr==nullptr) /* overflow */
      { assert(error::ERRNO!=0); return; }
    std::copy(d_ptr,d_ptr+first,new_ptr); // copy initial part from |*this|
    std::copy(source,source+r,new_ptr+first); // add final part from |source|
    memory::arena().free(d_ptr,d_allocated*sizeof(T));
    d_ptr = new_ptr;
    d_allocated = memory::arena().allocSize(first+r,sizeof(T));
    d_size = first+r;
    return;
  }

  // if we get here no new memory allocation is necessary

  if (d_size < first+r)
    setSize(first+r);

  std::copy(source,source+r,d_ptr+first); // add final part from |source|

  return;
}

template <class T> void List<T>::sort()

/*
  Sorts the list in the natural order of the elements. It is assumed that
  operator> is defined for T.
*/

{
  /* set the starting value of h */

  Ulong h = 1;

  for (; h < d_size/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < d_size; ++j) {
      T a = d_ptr[j];
      Ulong i = j;
      for (; (i >= h) && (d_ptr[i-h] > a); i -= h)
	d_ptr[i] = d_ptr[i-h];
      d_ptr[i] = a;
    }
  }

  return;
}

template<class T> template<class C> void List<T>::sort(C& c)

/*
  Sorts the list in the order defined by the comparison functor c. It is
  assumed that c takes two arguments of type T, and that c(x,y) is true
  if x <= y (so that the relation x > y is expressed by !c(c,y))
*/

{
  /* set the starting value of h */

  Ulong h = 1;

  for (; h < d_size/3; h = 3*h+1)
    ;

  /* do the sort */

  for (; h > 0; h /= 3) {
    for (Ulong j = h; j < d_size; ++j) {
      T a = d_ptr[j];
      Ulong i = j;
      for (; (i >= h) && !c(d_ptr[i-h],a); i -= h)
	d_ptr[i] = d_ptr[i-h];
      d_ptr[i] = a;
    }
  }

  return;
}

};


/**************************************************************************

        Chapter III -- Insertion and deletion

 A list will often be constructed by first catching the elements in a
 buffer, then copying the buffer onto the final list. For each type
 of coefficient, we will need an insertion function. We provide an
 abstract template.

 **************************************************************************/

namespace list {


/*
  Insert a new element in the (ordered) list, using binary search to find
  the insertion point.

  Forwards the error MEMORY_WARNING if CATCH_MEMORY_ERROR is set. Return
  value is not_found in case of error.

  NOTE :It is assumed that operator<= is defined for the class T.
*/
template <class T> Ulong insert(List<T>& l, const T& d_m)
{
  // this is necessary in the case where |d_m| points into the array being
  // resized. It could be avoided by doing the appendage after the allocation
  // of new memory and before the freeing of the old memory; this entails
  // not being able to use setSize. It hasn't seemed worthwile.
  T m = d_m;

  Ulong j0 = ~0L;
  Ulong j1 = l.size();

  while (j1-j0 > 1)
  {
    Ulong j = (j0+j1)/2;
    if (l[j] == m) /* m was found */
      return j;
    if (l[j] < m)
      j0 = j;
    else
      j1 = j;
  }

  // now |j1 == j0+1| and |j1==l.size() or l[j1]>m|; insertion point is |j1|

  l.setSize(l.size()+1);
  if (error::ERRNO) /* overflow */
    return not_found;
  l.setData(l.ptr()+j1,j1+1,l.size()-j1-1); /* shift tail up by one */
  new(l.ptr()+j1) T(m);

  return j1;
}


/*
  Finds the index of m in the list. If m is not found, returns not_found.
  Uses binary search.
*/

template <class T> Ulong find(const List<T>& l, const T& m)
{
  Ulong j0 = ~0L;

  for (Ulong j1 = l.size(); j1-j0 > 1;) {
    Ulong j = j0 + (j1-j0)/2;
    if (l[j] == m) /* m was found */
      return j;
    if (l[j] < m)
      j0 = j;
    else
      j1 = j;
  }

  return not_found;
}

};

/**************************************************************************

        Chapter IV -- Input/output

  This section defines some i/o functions for lists :

   - print(file,l) : prints l on the file;

 **************************************************************************/

namespace list {

template <class T> void print(FILE* file, const List<T>& l)

/*
  Rudimentary print function for lists. It assumes that print(FILE*,T) is
  defined.
*/

{
  for (Ulong j = 0; j < l.size(); ++j) {
    print(file,l[j]);
    if (j+1 < l.size()) /* more to come */
      fprintf(file,",");
  }
}

};

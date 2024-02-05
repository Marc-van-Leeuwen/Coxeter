/*
  This is dictionary.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef DICTIONARY_H  /* guard against multiple inclusions */
#define DICTIONARY_H

#include <string>
#include "globals.h"
#include "sl_list.h"

/******** type declarations *************************************************/

namespace dictionary {
  template <typename T> class Dictionary;
  template <typename T> class DictCell;
  template <typename T> class dict_const_iterator;
};

#include "memory.h"
#include "io.h"
#include <memory> // for |std::shared_ptr|

/******** function declarations *********************************************/

/* class definitions */

namespace dictionary {

template <typename T> class DictCell {
  friend class Dictionary<T>; // only that class can use our private members

  std::shared_ptr<T> ptr;
  std::unique_ptr<DictCell> left;
  std::unique_ptr<DictCell> right;
  char letter;
public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(DictCell));}
  DictCell(char c, std::shared_ptr<T> v,
	   DictCell *l = nullptr, DictCell *r = nullptr)
    : ptr(std::move(v)), left(l), right(r), letter(c) {};
  ~DictCell();
  void make_complete(); // install completions downwards
  bool has_own_action() const; // whether |ptr| defined, and not as completion
  const T* action() const { return ptr.get(); }
  containers::sl_list<std::string> extensions(std::string name) const;
};


template <typename T> class Dictionary
{
  DictCell<T> root_cell;
public:
/* creators and destructors */
  Dictionary(std::shared_ptr<T> v);
  virtual ~Dictionary() {}
/* modifiers */
  void insert(const std::string& str, std::shared_ptr<T> value);
  void remove(const std::string& str);
  void install_command_completion() { root_cell.make_complete(); };
  T* root_action() { return root_cell.ptr.get(); }
/* accessors */
  const T* find (const std::string& str, bool& absent_action) const;
  const DictCell<T>* findCell (const std::string& str) const;
  const DictCell<T>* root() const { return &root_cell; }

  // embedded class
  class const_iterator
    : public std::iterator<std::forward_iterator_tag, T, unsigned long>
  {
    // data
    DictCell<T>* p;
    containers::stack<DictCell<T>*> stack;

    using self = const_iterator;

  public:
    const_iterator() : p(nullptr), stack() {} // end indicator
    explicit const_iterator(const DictCell<T>* ptr)
      : p(const_cast<DictCell<T>*> (ptr)), stack() {}

    // contents access methods; return |const| ref/ptr for |const_iterator|
    const DictCell<T>& operator*()  const { return *p; }
    const DictCell<T>* operator->() const { return p; }

    // equality testing methods: test addresses of cells pointed to
    bool operator==(const self& x) const { return p == x.p; }
    bool operator!=(const self& x) const { return p != x.p; }

    self operator++();
    self operator++(int) { auto old = *this; operator++; return old; }
  }; // |class const_iterator|

  const_iterator begin() const { return const_iterator(&root_cell); }
  const_iterator end() const   { return const_iterator(); }
}; // |template <typename T> class Dictionary|


}; // namespace dictionary|

#include "dictionary.hpp"

#endif

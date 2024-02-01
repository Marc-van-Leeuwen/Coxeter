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

namespace dictionary {
  using namespace globals;
};

/******** type declarations *************************************************/

namespace dictionary {
  template <typename T> class Dictionary;
  template <typename T> struct DictCell;
};

#include "memory.h"
#include "io.h"
#include <memory> // for |std::shared_ptr|

namespace dictionary {
  using namespace memory;
  using namespace io;
};

/******** function declarations *********************************************/

namespace dictionary {
  template <typename T>
    void printExtensions(FILE* file, DictCell<T>* cell, std::string& name,
			 bool& first, const char* sep = ",");
};

/* class definitions */

namespace dictionary {

template <typename T> struct DictCell {
  std::shared_ptr<T> ptr;
  DictCell *left;
  DictCell *right;
  char letter;
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(DictCell));}
  DictCell(char c, std::shared_ptr<T> v,
	   DictCell *l = nullptr, DictCell *r = nullptr)
    :ptr(v), left(l), right(r), letter(c) {};
  ~DictCell();
  void make_complete(); // install completions dwonwards
  bool has_own_action() const; // whether |ptr| defined, and not as completion
};

template <typename T>
  class dict_const_iterator
  : public std::iterator<std::forward_iterator_tag, T, unsigned long>
{
  // data
  DictCell<T>* p;
  containers::stack<DictCell<T>*> stack;

  using self = dict_const_iterator<T>;

public:
  dict_const_iterator() : p(nullptr), stack() {} // end indicator
  explicit dict_const_iterator(const DictCell<T>* ptr)
    : p(const_cast<DictCell<T>*> (ptr)), stack() {}

  // contents access methods; return |const| ref/ptr for |const_iterator|
  const DictCell<T>& operator*()  const { return *p; }
  const DictCell<T>* operator->() const { return p; }

  // equality testing methods: test addresses of cells pointed to
  bool operator==(const self& x) const { return p == x.p; }
  bool operator!=(const self& x) const { return p != x.p; }

  self operator++()
  { if (p->left!=nullptr)
    {
      stack.push(p);
      p=p->left;
      return *this;
    }
    while (not stack.empty() and p->right==nullptr)
    {
      p=stack.top();
      stack.pop();
    }
    if (p->right==nullptr) // which implies |stack.empty|
    {
      p=nullptr; // indicates we are at the end
      return *this;
    }
    p=p->right;
    return *this;
  }

  self operator++(int)
  { auto old = *this;
    operator++;
    return old;
  }
};

template <typename T> class Dictionary {

 protected:
  DictCell<T>* d_root;
 public:
/* creators and destructors */
  Dictionary();
  virtual ~Dictionary();
/* modifiers */
  void insert(const std::string& str, std::shared_ptr<T> value);
  void remove(const std::string& str);
  void install_command_completion() { d_root->make_complete(); };
/* accessors */
  T* find (const std::string& str, bool& absent_action) const;
  DictCell<T>* findCell (const std::string& str) const;
  DictCell<T>* root() const { return d_root; }

  dict_const_iterator<T> begin() const
  { return dict_const_iterator<T>(d_root); }
  dict_const_iterator<T> end() const
  { return dict_const_iterator<T>(); }
};

};

#include "dictionary.hpp"

#endif

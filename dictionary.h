/*
  This is dictionary.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef DICTIONARY_H  /* guard against multiple inclusions */
#define DICTIONARY_H

#include <string>
#include "globals.h"

namespace dictionary {
  using namespace globals;
};

/******** type declarations *************************************************/

namespace dictionary {
  template <class T> class Dictionary;
  template <class T> struct DictCell;
};

#include "memory.h"
#include "io.h"

namespace dictionary {
  using namespace memory;
  using namespace io;
};

/******** function declarations *********************************************/

namespace dictionary {
  template <class T>
    void printExtensions(FILE* file, DictCell<T>* cell, std::string& name,
			 bool& first, const char* sep = ",");
};

/* class definitions */

namespace dictionary {

template <class T> struct DictCell {
  T *ptr;
  DictCell *left;
  DictCell *right;
  char letter;
  bool fullname;
  bool uniquePrefix;
/* constructors and destructors */
  void* operator new(size_t size) {return arena().alloc(size);}
  void operator delete(void* ptr)
    {return arena().free(ptr,sizeof(DictCell));}
  DictCell() {/* not implemented */};
  DictCell(char c, T* v, bool f, bool u, DictCell *l = 0, DictCell *r = 0)
    :ptr(v), left(l), right(r), letter(c), fullname(f), uniquePrefix(u) {};
  ~DictCell();
/* accessors */
  T* value() const {return ptr;}
};

template <class T> class Dictionary {

 protected:
  DictCell<T>* d_root;
 public:
/* creators and destructors */
  Dictionary();
  virtual ~Dictionary();
/* modifiers */
  void insert(const std::string& str, T* const value);
  void remove(const std::string& str);
/* accessors */
  T* find(const std::string& str) const;
  DictCell<T>* findCell(const std::string& str) const;
  DictCell<T>* root() {return d_root;}
};

};

#include "dictionary.hpp"

#endif

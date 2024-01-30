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
#include <memory> // for |std::shared_ptr|

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
  bool has_own_action() const; // whether |ptr| defined, and not as completion
};

template <class T> class Dictionary {

 protected:
  DictCell<T>* d_root;
 public:
/* creators and destructors */
  Dictionary();
  virtual ~Dictionary();
/* modifiers */
  void insert(const std::string& str, std::shared_ptr<T> value);
  void remove(const std::string& str);
/* accessors */
  T* find(const std::string& str, bool& absent_action) const;
  DictCell<T>* findCell(const std::string& str) const;
  DictCell<T>* root() {return d_root;}
};

};

#include "dictionary.hpp"

#endif

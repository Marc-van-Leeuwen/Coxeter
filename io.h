/*
  This is io.h

  Coxeter version 3.0  Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef IO_H  /* guarantee single inclusion */
#define IO_H

#include <string>
#include <cstring>
#include "globals.h"
#include "list.h"
#include "memory.h"

/******** type definitions **************************************************/

namespace io {

  /* style tags for i/o */

  struct Default{};
  struct GAP {};
  struct LaTeX {};
  struct Pretty {};
  struct Terse {};
  struct TeX {};
};

/******** constants **********************************************************/

namespace io {
  const Ulong LINESIZE = 79;
  const Ulong HALFLINESIZE = 39;
};

/******** function declarations **********************************************/

namespace io {
  int alphabeticDigits(Ulong c, Ulong b);
  std::string& append(std::string& l, const Ulong& n);
  std::string& append(std::string& l, const long& m);
  std::string& append(std::string& l, const int& n);
  std::string& append(std::string& l, const unsigned& n);
  int digits(Ulong c, Ulong b);
  void foldLine(FILE* file, const std::string& str, const Ulong& ls,
		const Ulong& h, const char* hyphens);
  // |getInput| will return |nullptr| in case it hits |EOF| immediately
  const char* getInput(FILE *inputfile, std::string& buf, Ulong len = 0);
  std::string& pad(std::string& l, const Ulong& n);
  void print(FILE* file, const char * str);                      /* inlined */
  void print(FILE* file, const std::string& str);                /* inlined */
  void print(FILE* file, const int *const& v, const Ulong& n);
  void printFile(FILE* file, const char *name);
  void printFile(FILE* file, const char *name, const char *dir_name);
  Ulong skipSpaces(const std::string& l, Ulong p);
  Ulong skipSpaces(const std::string& l, Ulong p);
};

/******** type definitions **************************************************/

/******** Inline definitions ***********************************************/

namespace io {

inline void print(FILE *file, const char* str) {fprintf(file,"%s",str);}
inline void print(FILE *file, const std::string& str)
  { fprintf(file,"%s",str.c_str()); }

};

#endif

/*
  This is coxtypes.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/****************************************************************************

  This file provides the basic type definitions for the Coxeter program.

 ****************************************************************************/

#ifndef COXTYPES_H  /* guard against multiple inclusions */
#define COXTYPES_H

#include "globals.h"
#include "constants.h"
#include "containers.h"

#include "io.h"
#include <limits.h>

/* type declarations and bounds */

namespace coxgroup {
  class CoxGroup; // forward declaration
};

namespace coxtypes {

  typedef unsigned short Rank;
  typedef Ulong CoxSize;         /* should hold at least 32 bits */
  typedef Ulong BettiNbr;        /* should hold at least 32 bits */
  using CoxNbr = unsigned int;   /* should fit into a CoxSize */
  typedef CoxNbr ParSize;          /* this should not be changed */
  typedef unsigned short ParNbr;   /* should fit into a CoxNbr */
  using CoxArr = ParNbr*; // points to an array of knwon size located somewhere
  using CoxLetter = unsigned char; /* for string representations */
  using Generator = CoxLetter;     /* internal representation of generators*/
  typedef unsigned short Length;
  typedef Ulong StarOp;          /* for numbering star operations */

  class CoxWord;

/* constants */

  // const Rank RANK_MAX = 255;          /* to enable string representations */
  // MEDRANK_MAX bits should fit in a Lflags
  const Rank MEDRANK_MAX = BITS(Ulong);
  // 2*SMALLRANK_MAX bits should fit in a Lflags
  const Rank SMALLRANK_MAX = BITS(Ulong)/2;
  const Rank RANK_MAX = SMALLRANK_MAX;              /* temporary restriction */
  const Generator GENERATOR_MAX = RANK_MAX-1;       /* top value is reserved */
  const Generator undef_generator = RANK_MAX;
  const CoxSize COXSIZE_MAX = ULONG_MAX-2;        /* top values are reserved */
  const CoxSize infinite_coxsize = COXSIZE_MAX+1;
  const CoxSize undef_coxsize = COXSIZE_MAX+2;
  const BettiNbr BETTI_MAX = ULONG_MAX-1;           /* top value is reserved */
  const BettiNbr undef_betti = BETTI_MAX+1;
  const CoxNbr COXNBR_MAX = UINT_MAX-1;             /* top value is reserved */
  const CoxNbr undef_coxnbr = COXNBR_MAX+1;
  const ParSize LPARNBR_MAX = COXNBR_MAX-RANK_MAX-1;/* top value is reserved */
  const ParNbr PARNBR_MAX = USHRT_MAX-RANK_MAX-1;   /* top value is reserved */
  const Length LENGTH_MAX = USHRT_MAX-1;            /* top value is reserved */
  const Length undef_length = LENGTH_MAX+1;
  const StarOp STAR_MAX = ULONG_MAX-1;              /* top value is reserved */
  const StarOp undef_starop = STAR_MAX+1;
};

/* functions provided by coxtypes.h */

namespace coxtypes {
  bool operator== (const CoxWord& g, const CoxWord& h);
  bool operator< (const CoxWord& g, const CoxWord& h);
  bool operator> (const CoxWord& g, const CoxWord& h);           /* inlined */
  std::string& append(std::string& str, const CoxNbr& x);
  void print(FILE *outputfile, CoxArr a, Rank l);
};

/******** Implementation ****************************************************/

#include "list.h"

namespace coxtypes {

struct Cox_word : public containers::vector<CoxLetter>
{
  using Base = containers::vector<CoxLetter>;
  using Base::Base;
}; // |class Cox_word|

class CoxWord {
 private:
  list::List<CoxLetter> d_list;
 public:
/* constructors */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(CoxWord));}
  CoxWord():d_list(0) {};  /* empty word */
  CoxWord(const Ulong& n);
  ~CoxWord();
/* accessors */
  const CoxLetter& operator[] (const Length& j) const { return d_list[j]; }
  Length length() const { return d_list.size()-1; }
/* modifiers */
  CoxWord& operator= (const CoxWord& h)
  { d_list.assign(h.d_list); return *this; }
  CoxLetter& operator[] (const Length& j) { return d_list[j]; }
  CoxWord& append(const CoxLetter& a);
  CoxWord& append(const CoxWord& h);
  CoxWord& erase(const Length& j);
  CoxWord& insert(const Length& j, const CoxLetter& a);
  CoxWord& reset();
  void setLength(Length n) { d_list.setSize(n+1); }
  CoxWord& setSubWord(const CoxWord& h, const Length& first,
		      const Length& r);
  Cox_word word () const
  { coxtypes::Cox_word result(d_list.begin(),d_list.end());
    for (auto& letter : result)
      --letter; // undo weird off-by-one encoding
    return result;
  }
};

};

/******** Inline definitions ***********************************************/

namespace coxtypes {

  inline bool operator> (const CoxWord& g, const CoxWord& h) {return h < g;}

  /* class CoxWord */


};

#endif

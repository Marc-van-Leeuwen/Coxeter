/*
  This is interface.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef INTERFACE_H  /* guard against multiple inclusions */
#define INTERFACE_H

#include <vector>
#include <algorithm>
#include "globals.h"

/******** type declarations *************************************************/

namespace interface {
  struct DescentSetInterface;
  struct GroupEltInterface;
  struct ReservedSymbols;
  class Interface;
  class TokenTree;

  struct TokenCell;
  struct ParseInterface;

  typedef unsigned int Token;

  // tags

  struct Hexadecimal {};
  struct HexadecimalFromZero {};
  struct Decimal {};
  struct Alphabetic {};
};

#include "automata.h"
#include "coxtypes.h"
#include "io.h"
#include "list.h"
#include "memory.h"
#include "minroots.h"
#include "transducer.h"


/******** constants **********************************************************/

namespace interface {

  const automata::Letter empty_type = 0;
  const automata::Letter generator_type = 1;
  const automata::Letter prefix_type = 2;
  const automata::Letter postfix_type = 3;
  const automata::Letter separator_type = 4;
  const automata::Letter modifier_type = 5;
  const automata::Letter grouping_type = 6;
  const automata::Letter number_type = 7;

};

/******** function declarations **********************************************/

namespace interface {
  const std::string* alphabeticSymbols(Ulong n);
  std::string& append(std::string& str, const coxtypes::CoxWord& g, const GroupEltInterface& GI);
  std::string& append(std::string& buf, const GenSet& f, const Interface& I);
  std::string& appendSymbol(std::string& str, const coxtypes::Generator& s, const Interface& I);
  std::string& appendTwosided(std::string& buf, const Lflags& f, const Interface& I);
  const std::string* checkLeadingWhite(const GroupEltInterface& GI);
  bool checkRepeated(const GroupEltInterface& GI);
  const std::string* checkReserved(const GroupEltInterface& GI, const Interface& I);
  const std::string* decimalSymbols(Ulong n);
  Ulong descentWidth(const Lflags& f, const Interface& I);
  const std::string* hexSymbols(Ulong n);
  const std::string* hexSymbolsFromZero(Ulong n);
  const bits::Permutation& identityOrder(Ulong n);
  bool isBeginGroup(const Token& tok);
  bool isContextNbr(const Token& tok);
  bool isDenseArray(const Token& tok);
  bool isEndGroup(const Token& tok);
  bool isInverse(const Token& tok);
  bool isLongest(const Token& tok);
  bool isModifier(const Token& tok);
  bool isPower(const Token& tok);
  void print(FILE *file, const coxtypes::CoxWord& g, const GroupEltInterface& I);
  void print(FILE *file,
	     const coxtypes::Cox_word& g, const GroupEltInterface& I);
  void print(FILE *file, const GenSet& f, const Interface& I);
  void print(FILE *file,
	     const GenSet& f, const DescentSetInterface& DI,
	     const GroupEltInterface& GI);
  void printSymbol(FILE *file, const coxtypes::Generator& s, const Interface& I);
  void printTwosided(FILE *file, const Lflags& f,
		     const DescentSetInterface& DI,
		     const GroupEltInterface& GI, const coxtypes::Rank& l);
  void printTwosided(FILE *file, const Lflags& f, const Interface& I);
                                                                 /* inlined */
  coxtypes::CoxNbr readCoxNbr(interface::ParseInterface& P, Ulong size);
  automata::Letter tokenType(const Token& tok);
  const std::string* twohexSymbols(Ulong n);
};

/******** type definitions ***************************************************/

namespace interface {

struct ParseInterface {
  std::string str;
  Ulong nestlevel;
  list::List<coxtypes::CoxWord> a;
  coxtypes::CoxWord c;
  automata::State x;
  Ulong offset;
/* constructors and destructors */
  ParseInterface();
  ~ParseInterface();
  void reset();
  coxtypes::Cox_word first_word() const { return a[0].word(); }
};

struct TokenCell {
  Token val;
  char letter;
  TokenCell *left;
  TokenCell *right;
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TokenCell));}
  TokenCell() {};
  ~TokenCell();
};

class TokenTree {
 private:
  TokenCell *d_root;
 public:
/* constructors and destructors */
  TokenTree();
  TokenTree(TokenCell *cell):d_root(cell) {};
  ~TokenTree();
/* manipulators */
  void insert(const std::string& str, const Token& val);
/* accessors */
  Ulong find(std::string& str, Token& val) const;
  Ulong find(const std::string& str, const Ulong& n, Token& val) const;
  TokenCell *root() {return d_root;}
};

struct DescentSetInterface {
  std::string prefix;
  std::string postfix;
  std::string separator;
  std::string twosidedPrefix;
  std::string twosidedPostfix;
  std::string twosidedSeparator;
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void* operator new(size_t size, void* ptr) {return ptr;}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(DescentSetInterface));}
  void operator delete(void* p1, void* p2) {};
  DescentSetInterface();
  DescentSetInterface(io::GAP);
  ~DescentSetInterface();
  void setPostfix(const std::string& str);
  void setPrefix(const std::string& str);
  void setSeparator(const std::string& str);
  void setTwosidedPrefix(const std::string& str);
  void setTwosidedPostfix(const std::string& str);
  void setTwosidedSeparator(const std::string& str);
};

struct GroupEltInterface {
  std::vector<std::string> symbol;
  std::string prefix;
  std::string postfix;
  std::string separator;
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GroupEltInterface));}
  GroupEltInterface();
  GroupEltInterface(const coxtypes::Rank& l);
  GroupEltInterface(const coxtypes::Rank& l, Alphabetic);
  GroupEltInterface(const coxtypes::Rank& l, Decimal);
  GroupEltInterface(const coxtypes::Rank& l, io::GAP);
  GroupEltInterface(const coxtypes::Rank& l, Hexadecimal);
  GroupEltInterface(const coxtypes::Rank& l, HexadecimalFromZero);
  ~GroupEltInterface();
  void setPostfix(const std::string& a);
  void setPrefix(const std::string& a);
  void setSeparator(const std::string& a);
  void setSymbol(const coxtypes::Generator& s, const std::string& a);
  void print(FILE* file) const;
};

struct ReservedSymbols {
  std::string beginGroup;
  std::string endGroup;
  std::string longest;
  std::string inverse;
  std::string power;
  std::string contextnbr;
  std::string densearray;
  ReservedSymbols();
  ReservedSymbols(io::Default);
  ~ReservedSymbols();
};

class Interface {
 protected:
  bits::Permutation d_order;
  TokenTree d_symbolTree;
  automata::Automaton const* d_tokenAut;
  GroupEltInterface* d_in;
  GroupEltInterface* d_out;
  DescentSetInterface* d_descent;
  std::string d_beginGroup;
  std::string d_endGroup;
  std::string d_longest;
  std::string d_inverse;
  std::string d_power;
  std::string d_contextNbr;
  std::string d_denseArray;
  std::string d_parseEscape;
  std::vector<std::string> d_reserved;
  coxtypes::Rank d_rank;
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(Interface));}
  Interface(const type::Type& x, const coxtypes::Rank& l);
  virtual ~Interface();
/* manipulators */
  void readSymbols();
  void setAutomaton();
  void setDescent(io::Default);
  void setDescent(io::GAP);
  virtual void setIn(const GroupEltInterface& i);                /* inlined */
  void setInPostfix(const std::string& a);                            /* inlined */
  void setInPrefix(const std::string& a);                             /* inlined */
  void setInSeparator(const std::string& a);                          /* inlined */
  void setInSymbol(const coxtypes::Generator& s, const std::string& a);         /* inlined */
  void setOrder(const bits::Permutation& gen_order);
  virtual void setOut(const GroupEltInterface& i);               /* inlined */
  void setOutPostfix(const std::string& a);                           /* inlined */
  void setOutPrefix(const std::string& a);                            /* inlined */
  void setOutSeparator(const std::string& a);                         /* inlined */
  void setOutSymbol(const coxtypes::Generator& s, const std::string& a);        /* inlined */
/* accessors */
  const DescentSetInterface& descentInterface() const;           /* inlined */
  Ulong getToken(interface::ParseInterface& P, Token& tok) const;         /* inlined */
  Ulong in(const Ulong& j) const;                            /* inlined */
  const GroupEltInterface& inInterface() const;                  /* inlined */
  const std::string& inPostfix() const;                               /* inlined */
  const std::string& inPrefix() const;                                /* inlined */
  const std::string& inSeparator() const;                             /* inlined */
  const std::string& inSymbol(const coxtypes::Generator& s) const;              /* inlined */
  bool isReserved(const std::string& str) const;                      /* inlined */
  const bits::Permutation& order() const;                              /* inlined */
  Ulong out(const Ulong& j) const;                           /* inlined */
  const GroupEltInterface& outInterface() const;                 /* inlined */
  const std::string& outPostfix() const;                              /* inlined */
  const std::string& outPrefix() const;                               /* inlined */
  const std::string& outSeparator() const;                            /* inlined */
  const std::string& outSymbol(const coxtypes::Generator& s) const;             /* inlined */
  bool parseCoxWord(interface::ParseInterface& P, const minroots::MinTable& T) const;
  coxtypes::Rank rank() const;                                             /* inlined */
  bool readCoxElt(interface::ParseInterface& P) const;
  const TokenTree& symbolTree() const;
// i/o
  virtual std::string& append(std::string& str, const coxtypes::CoxWord& g) const;
  virtual void print(FILE* file, const coxtypes::CoxWord& g) const
  { interface::print(file,g,*d_out); }
  void print(FILE* file, const coxtypes::Cox_word& g) const
  { interface::print(file,g,*d_out); }
}; // |class Interface|

}; // |namespace interface|

/******** inline implementations *******************************************/

namespace interface {

inline std::string& append(std::string& str, const coxtypes::CoxWord& g, const Interface& I)
  {return append(str,g,I.outInterface());}
inline std::string& appendSymbol(std::string& str, const coxtypes::Generator& s,
			    const Interface& I)
  { return str.append(I.outSymbol(s)); }
inline void print(FILE *file, const GenSet& f, const Interface& I)
  {return print(file,f,I.descentInterface(),I.outInterface());}
inline void printSymbol(FILE *file, const coxtypes::Generator& s, const Interface& I)
  {io::print(file,I.outSymbol(s));}
inline void printTwosided(FILE *file, const Lflags& f, const Interface& I)
  {return printTwosided(file,f,I.descentInterface(),I.outInterface(),
			I.rank());}

inline Ulong TokenTree::find(std::string& str, Token& val) const
  {return find(str.c_str(),str.length(),val);}

inline const DescentSetInterface& Interface::descentInterface() const
  {return *d_descent;}
inline Ulong Interface::getToken(interface::ParseInterface& P, Token& tok) const
  {return d_symbolTree.find(P.str,P.offset,tok);}
inline const GroupEltInterface& Interface::inInterface() const
  {return *d_in;}
inline const std::string& Interface::inPostfix() const {return d_in->postfix;}
inline const std::string& Interface::inPrefix() const {return d_in->prefix;}
inline const std::string& Interface::inSeparator() const {return d_in->separator;}
inline const std::string& Interface::inSymbol(const coxtypes::Generator& s) const
  {return d_in->symbol[s];}
inline bool Interface::isReserved(const std::string& str) const
 {return std::find(d_reserved.begin(),d_reserved.end(),str)!=d_reserved.end();}
inline const GroupEltInterface& Interface::outInterface() const
  {return *d_out;}
inline const bits::Permutation& Interface::order() const {return d_order;}
inline const std::string& Interface::outPostfix() const {return d_out->postfix;}
inline const std::string& Interface::outPrefix() const {return d_out->prefix;}
inline const std::string& Interface::outSeparator() const {return d_out->separator;}
inline const std::string& Interface::outSymbol(const coxtypes::Generator& s) const
  {return d_out->symbol[s];}
inline coxtypes::Rank Interface::rank() const {return d_rank;}
inline const TokenTree& Interface::symbolTree() const {return d_symbolTree;}

inline void Interface::setInPostfix(const std::string& a) {d_in->setPostfix(a);}
inline void Interface::setInPrefix(const std::string& a) {d_in->setPrefix(a);}
inline void Interface::setInSeparator(const std::string& a) {d_in->setSeparator(a);}
inline void Interface::setInSymbol(const coxtypes::Generator& s, const std::string& a)
  {return d_in->setSymbol(s,a);}
inline void Interface::setOutPostfix(const std::string& a) {d_out->setPostfix(a);}
inline void Interface::setOutPrefix(const std::string& a) {d_out->setPrefix(a);}
inline void Interface::setOutSeparator(const std::string& a)
  {d_out->setSeparator(a);}
inline void Interface::setOutSymbol(const coxtypes::Generator& s, const std::string& a)
  {return d_out->setSymbol(s,a);}

inline std::string& Interface::append(std::string& str, const coxtypes::CoxWord& g) const
  {return interface::append(str,g,*d_out);}

};

#endif

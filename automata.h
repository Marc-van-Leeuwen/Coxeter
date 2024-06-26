/*
  This is automata.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef AUTOMATA_H  /* guard against multiple inclusions */
#define AUTOMATA_H

#include "globals.h"

/******** type declarations **************************************************/

namespace automata {
  class Automaton;
  class ExplicitAutomaton;
  typedef unsigned Letter;
  typedef unsigned State;
};

/******** type definitions ***************************************************/

#include "bits.h"

namespace automata {

class Automaton {
 public:
/* accessors */
  virtual State act(State x, Letter a) const = 0;
  virtual State initialState() const = 0;
  virtual bool isAccept(State x) const = 0;
  virtual bool isFailure(State x) const = 0;
  virtual Ulong rank() const = 0;
  virtual Ulong size() const = 0;
};

class ExplicitAutomaton:public Automaton {
 private:
  State **d_table;
  bitmap::BitMap d_accept;
  State d_failure;
  State d_initial;
  Ulong d_rank;
  Ulong d_size;
 public:
/* constructors and destructors */
  ExplicitAutomaton(Ulong n, Ulong m);
  virtual ~ExplicitAutomaton();
/* manipulators */
  void setAccept(State x);                                        /* inlined */
  void setFailure(State x);                                       /* inlined */
  void setInitial(State x);                                       /* inlined */
  void setTable(State x, Letter a, State xa);                     /* inlined */
/* accessors */
  State act(State x, Letter a) const;                             /* inlined */
  State initialState() const;                                     /* inlined */
  bool isAccept(State x) const;                                   /* inlined */
  bool isFailure(State x) const;                                  /* inlined */
  Ulong rank() const;                                           /* inlined */
  Ulong size() const;                                           /* inlined */
};

};

/******** inline implementations ******************************************/

namespace automata {

  inline void ExplicitAutomaton::setAccept(State x) {d_accept.insert(x);}
  inline void ExplicitAutomaton::setFailure(State x) {d_failure = x;}
  inline void ExplicitAutomaton::setInitial(State x) {d_initial = x;}
  inline void ExplicitAutomaton::setTable(State x, Letter a, State xa)
    {d_table[x][a] = xa;}

  inline State ExplicitAutomaton::act(State x, Letter a) const
    {return d_table[x][a];}
  inline State ExplicitAutomaton::initialState() const {return d_initial;}
  inline bool ExplicitAutomaton::isAccept(State x) const
    {return d_accept.is_member(x);}
  inline bool ExplicitAutomaton::isFailure(State x) const
    {return x == d_failure;}
  inline Ulong ExplicitAutomaton::rank() const {return d_rank;}
  inline Ulong ExplicitAutomaton::size() const {return d_size;}
};

#endif

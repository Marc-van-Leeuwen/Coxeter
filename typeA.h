/*
  This is typeA.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/****************************************************************************

  This module declares some special features for groups of type A. Currently
  this is just at the level of i/o, being able to input and output elements
  as permutations.

 ****************************************************************************/

#ifndef TYPEA_H  /* guard against multiple inclusions */
#define TYPEA_H

#include "globals.h"
#include "fcoxgroup.h"

/******** type declarations *************************************************/

namespace typeA {
  class TypeAInterface;
  class TypeACoxGroup;
  class TypeABigRankCoxGroup;
  class GeneralTypeABRCoxGroup;
  class TypeAMedRankCoxGroup;
  class GeneralTypeAMRCoxGroup;
  class TypeASmallRankCoxGroup;
  class GeneralTypeASRCoxGroup;
  class TypeASmallCoxGroup;
  class GeneralTypeASCoxGroup;
};

/******** function declarations *********************************************/

namespace typeA {
  void coxWordToPermutation(coxtypes::CoxWord& a, const coxtypes::CoxWord& g);
  void permutationToCoxWord(coxtypes::CoxWord& g, const coxtypes::CoxWord& a);
};

/******** type definitions **************************************************/

namespace typeA {

class TypeACoxGroup:public fcoxgroup::FiniteCoxGroup {
  TypeAInterface* d_typeAInterface;
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TypeACoxGroup));}

  TypeACoxGroup(const coxtypes::Rank& l);
  virtual ~TypeACoxGroup();
// accessors
  void coxWordToPermutation(coxtypes::CoxWord& a, const coxtypes::CoxWord& g) const;
  bool hasPermutationInput() const;                              /* inlined */
  bool hasPermutationOutput() const;                             /* inlined */
  void permutationToCoxWord(coxtypes::CoxWord& g, const coxtypes::CoxWord& a) const;
  const TypeAInterface& typeAInterface() const;                  /* inlined */
// manipulators
  void setPermutationInput(bool b);                              /* inlined */
  void setPermutationOutput(bool b);                             /* inlined */
  TypeAInterface& typeAInterface();                              /* inlined */
// i/o
  virtual bool parseGroupElement(interface::ParseInterface& P) const;
};

class TypeABigRankCoxGroup:public TypeACoxGroup {
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TypeABigRankCoxGroup));}

  TypeABigRankCoxGroup(const coxtypes::Rank& l):TypeACoxGroup(l) {};
  virtual ~TypeABigRankCoxGroup() {};
};

class GeneralTypeABRCoxGroup:public TypeABigRankCoxGroup { // leaf class
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralTypeABRCoxGroup));}

  GeneralTypeABRCoxGroup(const coxtypes::Rank& l):TypeABigRankCoxGroup(l) {};
  ~GeneralTypeABRCoxGroup() {};
};

class TypeAMedRankCoxGroup:public TypeACoxGroup {
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TypeAMedRankCoxGroup));}

  TypeAMedRankCoxGroup(const coxtypes::Rank& l);
  virtual ~TypeAMedRankCoxGroup();
};

class GeneralTypeAMRCoxGroup:public TypeAMedRankCoxGroup { // leaf class
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralTypeAMRCoxGroup));}

  GeneralTypeAMRCoxGroup(const coxtypes::Rank& l):TypeAMedRankCoxGroup(l) {};
  ~GeneralTypeAMRCoxGroup() {};
};

class TypeASmallRankCoxGroup:public TypeAMedRankCoxGroup {
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TypeASmallRankCoxGroup));}

  TypeASmallRankCoxGroup(const coxtypes::Rank& l):TypeAMedRankCoxGroup(l) {};
  virtual ~TypeASmallRankCoxGroup() {};
};

class GeneralTypeASRCoxGroup:public TypeASmallRankCoxGroup { // leaf class
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralTypeASRCoxGroup));}

  GeneralTypeASRCoxGroup(const coxtypes::Rank& l):TypeASmallRankCoxGroup(l) {};
  ~GeneralTypeASRCoxGroup() {};
};

class TypeASmallCoxGroup: public TypeASmallRankCoxGroup {
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TypeASmallCoxGroup));}

  TypeASmallCoxGroup(const coxtypes::Rank& l):TypeASmallRankCoxGroup(l) {};
  virtual ~TypeASmallCoxGroup() {};
// accessors
  int prodD(coxtypes::CoxWord& g, const fcoxgroup::DenseArray& d_x) const;
// i/o
  bool parseDenseArray(interface::ParseInterface& P) const;
  virtual bool parseGroupElement(interface::ParseInterface& P) const;
};

class GeneralTypeASCoxGroup:public TypeASmallCoxGroup { // leaf class
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralTypeASCoxGroup));}

  GeneralTypeASCoxGroup(const coxtypes::Rank& l):TypeASmallCoxGroup(l) {};
  ~GeneralTypeASCoxGroup() {};
};

class TypeAInterface:public interface::Interface {
  interface::Interface* d_pInterface;
  bool d_hasPermutationInput;
  bool d_hasPermutationOutput;
 public:
// constructors and destructors
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(TypeAInterface));}

  TypeAInterface(const coxtypes::Rank& l);
  virtual ~TypeAInterface();
// accessors
  bool hasPermutationInput() const;                              /* inlined */
  bool hasPermutationOutput() const;                             /* inlined */
  bool parsePermutation(interface::ParseInterface& P) const;
// manipulators
  virtual void setIn(const interface::GroupEltInterface& i);
  virtual void setOut(const interface::GroupEltInterface& i);
  void setPermutationInput(bool b);                              /* inlined */
  void setPermutationOutput(bool b);                             /* inlined */
// i/o
  virtual std::string& append(std::string& str, const coxtypes::CoxWord& g) const;
  virtual void print(FILE* file, const coxtypes::CoxWord& g) const;
};

};

/******** inline definitions ************************************************/

namespace typeA {

inline bool TypeAInterface::hasPermutationInput() const
  {return d_hasPermutationInput;}
inline bool TypeAInterface::hasPermutationOutput() const
  {return d_hasPermutationOutput;}
inline void TypeAInterface::setPermutationInput(bool b)
  {d_hasPermutationInput = b;}
inline void TypeAInterface::setPermutationOutput(bool b)
  {d_hasPermutationOutput = b;}

inline bool TypeACoxGroup::hasPermutationInput() const
  {return d_typeAInterface->hasPermutationInput();}
inline bool TypeACoxGroup::hasPermutationOutput() const
  {return d_typeAInterface->hasPermutationOutput();}
inline void TypeACoxGroup::setPermutationInput(bool b)
  {typeAInterface().setPermutationInput(b);}
inline void TypeACoxGroup::setPermutationOutput(bool b)
  {typeAInterface().setPermutationOutput(b);}
inline const TypeAInterface& TypeACoxGroup::typeAInterface() const
  {return *d_typeAInterface;}
inline TypeAInterface& TypeACoxGroup::typeAInterface()
  {return *d_typeAInterface;}

};

#endif

/*
  This is general.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef GENERAL_H  /* guard against multiple inclusions */
#define GENERAL_H

#include "globals.h"
#include "coxgroup.h"

/******** type declarations **************************************************/

namespace general {
  class GeneralCoxGroup;
  class BigRankCoxGroup;
  class GeneralBRCoxGroup;
  class MedRankCoxGroup;
  class GeneralMRCoxGroup;
  class SmallRankCoxGroup;
  class GeneralSRCoxGroup;
};

/********* type definitions **************************************************/

namespace general {

class GeneralCoxGroup : public coxgroup::CoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralCoxGroup));}

  GeneralCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~GeneralCoxGroup();
/* accessors */
  virtual coxtypes::CoxSize order() const;                                 /* inlined */
};

class BigRankCoxGroup : public GeneralCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(BigRankCoxGroup));}
  BigRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~BigRankCoxGroup();
};

 class GeneralBRCoxGroup:public BigRankCoxGroup { /* leaf class */
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralBRCoxGroup));}
  GeneralBRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralBRCoxGroup();
};

class MedRankCoxGroup : public GeneralCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(MedRankCoxGroup));}
  MedRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~MedRankCoxGroup();
};

 class GeneralMRCoxGroup:public MedRankCoxGroup { /* leaf class */
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralMRCoxGroup));}
  GeneralMRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralMRCoxGroup();
};

class SmallRankCoxGroup : public MedRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(SmallRankCoxGroup));}
  SmallRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~SmallRankCoxGroup();
};

 class GeneralSRCoxGroup:public SmallRankCoxGroup { /* leaf class */
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralSRCoxGroup));}
  GeneralSRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralSRCoxGroup();
};

};

/******** inline definitions *************************************************/

namespace general {

inline coxtypes::CoxSize GeneralCoxGroup::order() const
  {return coxtypes::undef_coxsize;}

};

#endif

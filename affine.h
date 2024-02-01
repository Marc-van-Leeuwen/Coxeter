/*
  This is affine.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef AFFINE_H  /* guarantee single inclusion */
#define AFFINE_H

#include "globals.h"
#include "coxgroup.h"

/******** type declarations *************************************************/

namespace affine {
  class AffineCoxGroup;
  class AffineBigRankCoxGroup;
  class GeneralABRCoxGroup;
  class AffineMedRankCoxGroup;
  class GeneralAMRCoxGroup;
  class AffineSmallRankCoxGroup;
  class GeneralASRCoxGroup;
};

/******** type definitions **************************************************/

namespace affine {

class AffineCoxGroup : public coxgroup::CoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(AffineCoxGroup));}

  AffineCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~AffineCoxGroup();
/* accessors */
  coxtypes::CoxSize order() const;                                         /* inlined */
};

class AffineBigRankCoxGroup : public AffineCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(AffineBigRankCoxGroup));}
  AffineBigRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~AffineBigRankCoxGroup();
};

class GeneralABRCoxGroup:public AffineBigRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralABRCoxGroup));}
  GeneralABRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralABRCoxGroup();
};

class AffineMedRankCoxGroup : public AffineCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(AffineMedRankCoxGroup));}
  AffineMedRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~AffineMedRankCoxGroup();
};

class GeneralAMRCoxGroup:public AffineMedRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralAMRCoxGroup));}
  GeneralAMRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralAMRCoxGroup();
};

class AffineSmallRankCoxGroup : public AffineMedRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(AffineSmallRankCoxGroup));}
  AffineSmallRankCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  virtual ~AffineSmallRankCoxGroup();
};

class GeneralASRCoxGroup:public AffineSmallRankCoxGroup {
 public:
/* constructors and destructors */
  void* operator new(size_t size) {return memory::arena().alloc(size);}
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(GeneralASRCoxGroup));}
  GeneralASRCoxGroup(const type::Type& x, const coxtypes::Rank& l);
  ~GeneralASRCoxGroup();
};

};

/******** Inline implementations ******************************************/

namespace affine {

inline coxtypes::CoxSize AffineCoxGroup::order() const {return coxtypes::infinite_coxsize;}

};

#endif

/*
  This is hecke.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef HECKE_H  /* guard against multiple inclusions */
#define HECKE_H

#include "globals.h"
#include "list.h"

/******** type declarations *************************************************/

namespace hecke {
  template<class P> class HeckeMonomial;
  template<class P> struct NFCompare;
  template<class P> class HeckeIterator;
  template<class P> class ToCoxNbr;
};

/******** function declarations *********************************************/

#include "interface.h"
#include "schubert.h"

namespace hecke {
  template<class P>
  void append(std::string& str, const HeckeMonomial<P>& m,
	      const schubert::SchubertContext& p,
	      const interface::Interface& I);
  template<class P>
  void prettyPrint(FILE* file, const list::List<HeckeMonomial<P> >& h,
		   const bits::Permutation& a, const schubert::SchubertContext& p,
		   const interface::Interface& I, const coxtypes::Length& l,
		   const Ulong &ls = io::LINESIZE);
  template<class P>
  void printBasis(FILE* f, const list::List<HeckeMonomial<P> >& h,
		  const interface::Interface& I);
  template<class P>
  void singularStratification(list::List<HeckeMonomial<P> >& hs,
			      const list::List<HeckeMonomial<P> >& h,
			      const schubert::SchubertContext& p);
};

/******** type definitions **************************************************/

namespace hecke {

template<class P> class HeckeMonomial {
 private:
  coxtypes::CoxNbr d_x;
  const P* d_pol;
 public:
  typedef P PolType;
/* constructors and destructors */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(HeckeMonomial));}
  HeckeMonomial() {};
  HeckeMonomial(const coxtypes::CoxNbr& x, const P* pol);
  ~HeckeMonomial();
/* accessors */
  bool operator> (const HeckeMonomial& m);
  const P& pol() const;                                          /* inlined */
  coxtypes::CoxNbr x() const;                                              /* inlined */
/* manipulators */
  void setData(const coxtypes::CoxNbr& x, const P* pol);
};

template<class P> class ToCoxNbr {
 private:
  const list::List<HeckeMonomial<P> >* d_h;
 public:
  ToCoxNbr(const list::List<HeckeMonomial<P> >* h):d_h(h) {};
  ~ToCoxNbr() {};
  coxtypes::CoxNbr operator() (const Ulong& j) {return (*d_h)[j].x();}
};

template<class P> struct NFCompare {
  const schubert::SchubertContext& p;
  const bits::Permutation& order;
  NFCompare
    (const schubert::SchubertContext& q, const bits::Permutation& generator_ordering)
    :p(q),order(generator_ordering) {};
  ~NFCompare() {};
  bool operator()(const HeckeMonomial<P>& a, const HeckeMonomial<P>& b) const
    {return shortLexOrder(p,a.x(),b.x(),order);}
};

};

/******** inline definitions ************************************************/

namespace hecke {

template<class P>
inline bool HeckeMonomial<P>::operator> (const HeckeMonomial<P>& m)
  {return d_x > m.d_x;}
template<class P> inline const P& HeckeMonomial<P>::pol() const
  {return *d_pol;}
template<class P> inline coxtypes::CoxNbr HeckeMonomial<P>::x() const
  {return d_x;}
template<class P>
inline void HeckeMonomial<P>::setData(const coxtypes::CoxNbr& x, const P* pol)
  {d_x = x; d_pol = pol;}

};

#include "hecke.hpp"

#endif

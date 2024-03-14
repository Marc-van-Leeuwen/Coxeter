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
		   const bits::Permutation& a,
		   const schubert::SchubertContext& p,
		   const interface::Interface& I, const coxtypes::Length& l,
		   const Ulong &ls = io::LINESIZE);
  template<class P>
  void printBasis(FILE* f, const list::List<HeckeMonomial<P> >& h,
		  const interface::Interface& I);
  template<class P>
  containers::vector<HeckeMonomial<P> > singular_stratification
    (const schubert::SchubertContext& p,
     const containers::vector<HeckeMonomial<P> >& h);
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
  bool operator> (const HeckeMonomial& m) const { return d_x > m.d_x; }
  bool operator< (const HeckeMonomial& m) const { return d_x < m.d_x; }
  const P& pol() const       { return *d_pol; }
  coxtypes::CoxNbr x() const { return d_x; }
/* manipulator */
  void setData(const coxtypes::CoxNbr& x, const P* pol) { d_x = x; d_pol = pol; }
};

template<class P> class ToCoxNbr {
 private:
  const containers::vector<HeckeMonomial<P> >* d_h;
 public:
  ToCoxNbr(const containers::vector<HeckeMonomial<P> >* h):d_h(h) {};
  ~ToCoxNbr() {};
  coxtypes::CoxNbr operator() (const Ulong& j) {return (*d_h)[j].x();}
};

template<class P> struct NFCompare {
  const schubert::SchubertContext& p;
  const bits::Permutation& order;
  NFCompare
    (const schubert::SchubertContext& q,
     const bits::Permutation& generator_ordering)
    :p(q),order(generator_ordering) {};
  ~NFCompare() {};
  bool operator()(const HeckeMonomial<P>& a, const HeckeMonomial<P>& b) const
    {return shortLexOrder(p,a.x(),b.x(),order);}
}; // |template<class P> struct NFCompare|

}; // |namespace hecke|

#include "hecke.hpp"

#endif

/*
  This is hecke.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#include "polynomials.h"

namespace hecke {

/*****************************************************************************

  This module provides some definitions and functions (mostly output-oriented)
  for Hecke algebra elements. We have refrained in this program from actually
  implementing the Hecke algebra structure, which could easily lead to
  immense computations. Our main purpose is to have convenient containers
  for communicating data.

  (For this reason also our coefficients are the ordinary k-l polynomials
  in q, instead of Laurent polynomials in q^{1/2}, as they should be.)

  Things are written as templates because various types of coefficients could
  (and do, in this program) occur.

******************************************************************************/

namespace {

  template<class P> struct PPtrF {
    typedef const P* valueType;
    valueType operator() (const HeckeMonomial<P>& m) const {return &m.pol();}
  };

  template<class P>
  void appendStar(std::string& str, const HeckeMonomial<P>& m,
		  const schubert::SchubertContext& p, const coxtypes::Length& l);
  template<class P>
  Ulong maxLength(const list::List<HeckeMonomial<P> >& h, const schubert::SchubertContext& p,
		    const interface::Interface& I, const coxtypes::Length& l);
  template<class P>
  void oneColumnPrint(FILE* file, const list::List<HeckeMonomial<P> >& h,
		      const bits::Permutation& a, const schubert::SchubertContext& p,
		      const interface::Interface& I, const coxtypes::Length& l, const Ulong& ls);
  template<class P>
  void twoColumnPrint(FILE* file, const list::List<HeckeMonomial<P> >& h,
		      const bits::Permutation&a, const schubert::SchubertContext& p,
		      const interface::Interface& I, const coxtypes::Length& l, const Ulong& ls);
}; // |namespace|

/*****************************************************************************

        Chapter II -- The HeckeMonomial class.

  HeckeMonomial's are simply building blocks for HeckeElt's. They contain
  a context number and a polynomial reference.

  The following functions are defined :

   - constructors and destructors :

     - HeckeMonomial(x,pol);
     - ~HeckeMonomial();

   - accessors :

   - manipulators :

     - order(c) : orders according to c;

 *****************************************************************************/

template <class P>
HeckeMonomial<P>::HeckeMonomial(const coxtypes::CoxNbr& x, const P* pol)
  :d_x(x), d_pol(pol)
{}

template<class P> HeckeMonomial<P>::~HeckeMonomial()
{}


/*****************************************************************************

        Chapter III -- Utilities.

  This section defines some utility functions declared in hecke.h :

   - append(str,m,p,I) : appends m to str using I;
   - appendStar(str,m,p,l) : appends a star to str if there is a
     mu-coefficient;
   - maxLength(h,p,I,l) : computes the maximal length of an output line;
   - oneColumnPrint(h,p,I,l,ls) : does one-column output;
   - prettyPrint(f,h,I,l) : pretty-prints h on f, using I;
   - singular_stratification(hs,h,p) : puts in hs the singular stratification
     of h;
   - singularLocus(hs,h,p) : puts in hs the rational singular locus of h;
   - twoColumnPrint(h,p,I,l,ls) : does two-column output;

 *****************************************************************************/


template<class P>
void append(std::string& str, const HeckeMonomial<P>& m, const schubert::SchubertContext& p,
	    const interface::Interface& I)

/*
  Outputs m to str.
*/

{
  p.append(str,m.x(),I);
  str.append(" : ");
  polynomials::append(str,m.pol(),"q");

  return;
}

namespace {

template<class P>
void appendStar(std::string& str, const HeckeMonomial<P>& m,
		  const schubert::SchubertContext& p, const coxtypes::Length& l)

{
  coxtypes::Length lx = p.length(m.x());

  if (static_cast<long>(2*m.pol().deg()) == static_cast<long>(l-lx-1))
    str.append(" *");

  return;
}

template<class P>
Ulong maxLength(const list::List<HeckeMonomial<P> >& h, const schubert::SchubertContext& p,
		  const interface::Interface& I, const coxtypes::Length& l)

/*
  Returns the length of the longest line that would be printed out by
  oneColumnPrint(file,h,I,l). This is a preliminary to prettyprinting.
*/

{
  static std::string buf;

  Ulong maxl = 0;

  for (Ulong j = 0; j < h.size(); ++j) {
    buf.clear();
    const HeckeMonomial<P>& m = h[j];
    hecke::append(buf,m,p,I);
    appendStar(buf,m,p,l);
    if (maxl < buf.length())
      maxl = buf.length();
  }

  return maxl;
}

template<class P>
void oneColumnPrint(FILE* file,
		    const list::List<HeckeMonomial<P> >& h,
		    const bits::Permutation& a,
		    const schubert::SchubertContext& p,
		    const interface::Interface& I,
		    const coxtypes::Length& l,
		    const Ulong& ls)

/*
  This function prints out the row in one-column format, trying to fold long
  lines decently. The width of the column is given by ls.
*/

{
  static std::string buf;

  for (Ulong j = 0; j < h.size(); ++j) {
    buf.clear();
    hecke::append(buf,h[a[j]],p,I);
    appendStar(buf,h[a[j]],p,l);
    io::foldLine(file,buf,ls,4,"+");
    fprintf(file,"\n");
  }
}

};

template<class P>
void prettyPrint(FILE* file,
		 const list::List<HeckeMonomial<P> >& h,
		 const bits::Permutation& a,
		 const schubert::SchubertContext& p,
		 const interface::Interface& I,
		 const coxtypes::Length& l,
		 const Ulong& ls)

/*
  This function does the prettyprinting of h to the file. The formatting
  of the output is optimized for screen viewing. This means that if two
  entries fit on a line, we will do two-column output. Otherwise, we do
  one-column output, and moreover we try to fold long lines decently.

  The parameter l is needed to determine the non-zero mu-coefficients.
*/

{
  static std::string buf;

  Ulong maxl = maxLength(h,p,I,l);
  Ulong hl = (ls-1)/2;

  if (maxl > hl)
    return oneColumnPrint(file,h,a,p,I,l,ls);
  else
    return twoColumnPrint(file,h,a,p,I,l,ls);

  return;
}


/*
  This function extracts the "rational singular stratification". By this we
  mean that we sort by Kazhdan-Lusztig polynomials, and then consider maximal
  elements (for the Bruhat ordering) in each class.

  Geometrically, when the Bruhat ordering comes from the stratification
  of a Schubert variety cl(X_y), and the row is the extremal row for y,
  this means that we are looking at a version of "equisingularity" (the
  Kazhdan-Lusztig polynomial P_{x,y} being a measure of the failure of
  smoothness along the subvariety X_x), and taking maximal elements amounts
  to taking components of the equisingular locus. Note that as P_{x,y} is
  constant on each orbit in [e,y] under the descent set of y, these
  components always correspond to extremal elements.

  This is the set of data that has been popularized by Goresky in the
  files on his website. It is printed out in printRow.

  Note that from Irving (Ann. ENS ...) it is known that P_{z,y} <= P_{x,y}
  coefficientwise when x <= z, so P_{x,y} is a decreasing function of x,
  in the case of finite Weyl groups; presumably this is also known for
  general crystallographic Coxeter groups (= Weyl groups of Kac-Moody
  algebras).

  It is assumed that row is sorted in ShortLex order. The row is also
  returned sorted in ShortLex order.
*/
template<class P>
containers::vector<HeckeMonomial<P> > singular_stratification
  (const schubert::SchubertContext& p,
   const containers::vector<HeckeMonomial<P> >& h)
{
  // sort row by kl-polynomial
  PPtrF<P> f;
  bits::Partition pi(h.begin(),h.end(),f);

  // find maximal elements in each class
  containers::vector<HeckeMonomial<P> > result; // final size cannot be predicted
  for (bits::PartitionIterator pit(pi); pit; ++pit)
  {
    const bits::Set& pi_class = pit();
    Ulong m = pi_class[0];
    if (h[m].pol().deg() == 0) // polynomial is one
      continue;
    containers::vector<coxtypes::CoxNbr> c;
    c.reserve(pi_class.size());
    for (auto elt : pi_class)
      c.push_back(h[elt].x());
    for (auto i : indices_of_maxima(p,c))
      result.push_back(h[pi_class[i]]);
  }
  return result;
}

namespace {

template<class P>
void twoColumnPrint(FILE* file,
		    const list::List<HeckeMonomial<P> >& h,
		    const bits::Permutation& a,
		    const schubert::SchubertContext& p,
		    const interface::Interface& I,
		    const coxtypes::Length& l,
		    const Ulong& ls)

/*
  This function prints out the row in two-column format, on lines of length
  ls. It is assumed that it has been checked (using maxLength for instance)
  that the maximum size of an output line in print(file,kl,row) is at most
  (ls-1)/2 (so that there is room for at least one unit of whitespace
  in-between columns.)
*/

{
  static std::string buf;

  Ulong hl = (ls-1)/2; /* width of output column */
  Ulong fl = h.size()/2; /* number of full lines */
  Ulong i = 0;

  for (Ulong j = 0; j < fl; ++j) { /* print out a full line */
    buf.clear();
    hecke::append(buf,h[a[i]],p,I);
    appendStar(buf,h[a[i]],p,l);
    io::pad(buf,ls-hl);
    i++;
    hecke::append(buf,h[a[i]],p,I);
    appendStar(buf,h[a[i]],p,l);
    i++;
    io::print(file,buf);
    fprintf(file,"\n");
  }

  if (h.size()%2) { /* print out a half line */
    buf.clear();
    hecke::append(buf,h[a[i]],p,I);
    appendStar(buf,h[a[i]],p,l);
    io::print(file,buf);
    fprintf(file,"\n");
  }

  return;
}

}; // |namespace|

}; // |namspace hecke|

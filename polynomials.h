/*
  This is polynomials.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef POLYNOMIALS_H  /* include guard */
#define POLYNOMIALS_H

#include "globals.h"
#include <limits>

#include "containers.h"
#include "io.h"

namespace polynomials {

/******** type declarations **************************************************/

  typedef Ulong Degree;
  typedef long SDegree;
  template <class T> class Polynomial;
  template <class T> class LaurentPolynomial;

/******** constants **********************************************************/

  static constexpr Degree undef_degree = ~0;
  static constexpr Degree DEGREE_MAX = std::numeric_limits<Ulong>::max()-1;

  static constexpr SDegree SDEGREE_MAX = std::numeric_limits<long>::max();
  static constexpr SDegree SDEGREE_MIN = std::numeric_limits<long>::min()+1;

/******** function definitions ***********************************************/

  template <class T>
  std::string& append(std::string& str, const Polynomial<T> &p, const char *x);
  template <class T>
  std::string& append(std::string& str, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x);
  template <class T>
  std::string& append(std::string& str, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,io::GAP);
  template <class T>
  std::string& append(std::string& str, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,io::Terse);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const char *x);
  template <class T>
  void print(FILE* file, const LaurentPolynomial<T>& p, const char *x);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,io::GAP);
  template <class T>
  void print(FILE* file, const Polynomial<T>& p, const Degree& d,
	     const long& m, const char *x,io::Terse);
  template <class T>
  SDegree sumDegree(const LaurentPolynomial<T>& p,
		    const LaurentPolynomial<T>& q);
  template <class T>
  SDegree sumValuation(const LaurentPolynomial<T>& p,
		       const LaurentPolynomial<T>& q);

/******** type definitions ***************************************************/


template <class T> class Polynomial {

  containers::vector<T> v;
 public:
  typedef struct {} const_tag;
/* constructors and destructors */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(Polynomial<T>));}
  Polynomial<T>():v() {} // polynomial with empty coefficient vector (is zero)
  Polynomial<T>(Degree d) : v(d+1,0) {} // prepare for degree |d| but zero coefs
  Polynomial<T>(const Polynomial<T>& q) = default;
  Polynomial<T>(T* const& ptr, const Degree& d) : v(ptr,ptr+d+1) {} // fill up
  Polynomial<T>(const T& c, const_tag): v{c} {} // degree 0 polynomial, coef |c|

/* manipulators */
  T& operator[] (const Ulong& j) { return v[j]; }
  void reduceDeg() { while (not v.empty() and v.back()==T(0)) v.pop_back(); }
  void setDeg(const Degree& d)  { v.resize(d+1,T(0)); }
  void setDegValue(const Degree& d) { v.resize(d+1,T(0)); } // no difference
  void setVect(const T *source, const Ulong& n) { v.assign(source,source+n); }
  void setZero() { std::fill(v.begin(),v.end(),0); } // empty coefficient vector
  void setZero(const Ulong& r) { std::fill(&v[0],&v[r],0); } // zero out
  void setZero(const Ulong& first, const Ulong& r)
    { std::fill(&v[first],&v[first+r],0); }

/* accessors */
  const T& operator[] (const Ulong& j) const { return v[j]; }
  Ulong deg() const { return v.size()-1; }
  bool isZero() const { return v.empty(); } // not true after |setZero|!

/* operators and operations */
  Polynomial<T>& operator= (const Polynomial<T>& q) = default;
  Polynomial<T>& operator+= (const Polynomial<T>& q);
  Polynomial<T>& operator-= (const Polynomial<T>& q);
  Polynomial<T>& operator*= (const T& a);
  Polynomial<T>& operator*= (const Polynomial<T>& q);
  Polynomial<T>& operator/= (const Polynomial<T>& q);

  bool operator== (const Polynomial<T>& q) const { return v==q.v; }
  bool operator!= (const Polynomial<T>& q) const { return v!=q.v; }
  bool operator<  (const Polynomial<T>& q) const
  { return v.size()==q.v.size() ? v < q.v : v.size()<q.v.size(); }
  bool operator>  (const Polynomial<T>& q) const { return q < *this; }
  bool operator<= (const Polynomial<T>& q) const { return not (q < *this); }
  bool operator>= (const Polynomial<T>& q) const { return not (*this < q); }
}; // template |class Polynomial<T>|

template <class T> class LaurentPolynomial {
  // interpreting |coef| as polynomial in $X$, represents $coef*X^{d_valuation}$
  // convention: zero polynomial has no coefficients and (ignored) valuation 0
  containers::vector<T> coef;
  SDegree d_valuation; /* degree of first non-zero coefficient */
 public:
/* constructors and destructors */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(LaurentPolynomial<T>));}
  LaurentPolynomial<T>() : coef(), d_valuation(0) {}
  LaurentPolynomial<T>(const SDegree& degree, const SDegree& offset = 0)
    : coef(degree-offset+1,T(0)),d_valuation(offset)
  {}
  LaurentPolynomial<T>(containers::vector<T>&& c, const SDegree& offset = 0)
    : coef(std::move(c)),d_valuation(offset) {}

/* accessors */
  SDegree deg() const { return coef.size()-1+d_valuation; }
  bool isZero() const { return coef.empty(); }
  SDegree val() const { return d_valuation; }

  const T& operator[] (const SDegree& j) const // subscript: shifted by |val()|
  { return coef[j-d_valuation]; }

  bool operator== (const LaurentPolynomial& q) const
  { return isZero ? q.isZero() : val()==q.val() and coef==q.coef; }
  bool operator!= (const LaurentPolynomial& q) const{ return not operator==(q); }
  bool operator< (const LaurentPolynomial& q) const;
  bool operator>  (const LaurentPolynomial<T>& q) const { return q < *this; }
  bool operator<= (const LaurentPolynomial<T>& q) const
  { return not (q < *this); }
  bool operator>= (const LaurentPolynomial<T>& q) const
  { return not (*this < q); }

/* manipulators */
  T& operator[] (const SDegree& j) { return coef[j-d_valuation]; }

}; // template |class LaurentPolynomial<T>|


}; // |namespace polynomials|

#include "polynomials.hpp"

#endif

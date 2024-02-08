/*
  This is polynomials.h

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

#ifndef POLYNOMIALS_H  /* include guard */
#define POLYNOMIALS_H

#include "globals.h"
#include <limits>
#include <vector>

#include "io.h"
#include "vector.h"

namespace polynomials {

/******** type declarations **************************************************/

  typedef Ulong Degree;
  typedef long SDegree;
  template <class T> class Polynomial;
  template <class T> class LaurentPolynomial;

/******** constants **********************************************************/

  static constexpr Degree undef_degree = ~0;
  static constexpr Degree DEGREE_MAX = std::numeric_limits<Ulong>::max()-1;

  static constexpr SDegree undef_valuation = std::numeric_limits<long>::min();
  static constexpr SDegree SDEGREE_MAX = std::numeric_limits<long>::max();
  static constexpr SDegree SDEGREE_MIN = std::numeric_limits<long>::min()+1;

/******** function definitions ***********************************************/

  template <class T>
  bool operator== (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator!= (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator<= (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator>= (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator< (const Polynomial<T>& p, const Polynomial<T>& q);
  template <class T>
  bool operator> (const Polynomial<T>& p, const Polynomial<T>& q);
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
 protected:
  vector::Vector<T> v;
 public:
  typedef struct {} const_tag;
/* constructors and destructors */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(Polynomial<T>));}
  Polynomial<T>(){};
  Polynomial<T>(Degree d):v(d+1) {};
  Polynomial<T>(const Polynomial<T>& q):v(q.v) {};
  Polynomial<T>(T* const& ptr, const Degree& d):v(ptr,d+1) {};
  Polynomial<T>(const T& c, const_tag):v(1) {v[0] = c; setDegValue(0);}
  ~Polynomial<T>();
/* manipulators */
  T& operator[] (const Ulong& j);                                /* inlined */
  void reduceDeg();                                              /* inlined */
  void setDeg(const Degree& d);                                  /* inlined */
  void setDegValue(const Degree& d);                             /* inlined */
  void setVect(const T *source, const Ulong& n);                 /* inlined */
  void setZero();                                                /* inlined */
  void setZero(const Ulong& r);                                  /* inlined */
  void setZero(const Ulong& first, const Ulong& r);              /* inlined */
  vector::Vector<T>& vect();                                             /* inlined */
/* accessors */
  const T& operator[] (const Ulong& j) const;                    /* inlined */
  Ulong deg() const;                                             /* inlined */
  bool isZero() const;                                           /* inlined */
  const vector::Vector<T>& vect() const;                                 /* inlined */
/* operators and operations */
  Polynomial<T>& operator= (const Polynomial<T>& q);             /* inlined */
  Polynomial<T>& operator+= (const Polynomial<T>& q);
  Polynomial<T>& operator-= (const Polynomial<T>& q);
  Polynomial<T>& operator*= (const T& a);
  Polynomial<T>& operator*= (const Polynomial<T>& q);
  Polynomial<T>& operator/= (const Polynomial<T>& q);
}; // template |class Polynomial<T>|

template <class T> class LaurentPolynomial {
  // interpreting |coef| as polynomial in $X$, represents $coef*X^{d_valuation}$
  std::vector<T> coef;
  SDegree d_valuation; /* degree of first non-zero coefficient */
 public:
/* constructors and destructors */
  void operator delete(void* ptr)
    {return memory::arena().free(ptr,sizeof(LaurentPolynomial<T>));}
  LaurentPolynomial<T>() : coef(), d_valuation(undef_valuation) {};

  LaurentPolynomial<T>(const SDegree& degree, const SDegree& offset = 0)
    : coef(degree-offset+1,T(0)),d_valuation(offset)
  {}
  LaurentPolynomial<T>(std::vector<T>&& c, const SDegree& offset = 0)
    : coef(std::move(c)),d_valuation(offset) {}
  ~LaurentPolynomial<T>() {}
/* accessors */
  SDegree deg() const { return coef.size()-1+d_valuation; }
  bool isZero() const { return coef.empty(); }
  SDegree val() const { return d_valuation; }

  const T& operator[] (const SDegree& j) const { return coef[j-d_valuation]; }

  bool operator== (const LaurentPolynomial& p) const;
  bool operator!= (const LaurentPolynomial& p) const{ return not operator==(p); }
  bool operator<= (const LaurentPolynomial& p) const;
  bool operator>= (const LaurentPolynomial& p) const;
  bool operator< (const LaurentPolynomial& p) const { return not operator>=(p); }
  bool operator> (const LaurentPolynomial& p) const { return not operator<=(p); }

/* manipulators */
  T& operator[] (const SDegree& j) { return coef[j-d_valuation]; }

  void setBounds(const SDegree& deg, const SDegree& val)
  {  coef.resize(deg-val+1,T(0)); d_valuation = val; }
  void setDeg(const SDegree& n) { coef.resize(n-d_valuation+1,T(0)); }
  void setDegValue(const SDegree& n) { coef.resize(n-d_valuation+1); }
  void setVal(const SDegree& n) // prepare for having valuation |n|, old degree
  { d_valuation=n; coef.resize(coef.size()-n,T(0)); } // ignores old |coef|s
  void setValValue(const SDegree& n) { d_valuation = n;}
  void setZero() { coef.resize(0); }
                                                 /* inlined */
}; // template |class LaurentPolynomial<T>|


/******** inline definitions **************************************************/


template <class T>
inline bool operator!= (const Polynomial<T>& p, const Polynomial<T>& q)
  {return !(p == q);}
template <class T>
inline bool operator< (const Polynomial<T>& p, const Polynomial<T>& q)
  {return !(p >= q);}
template <class T>
inline bool operator> (const Polynomial<T>& p, const Polynomial<T>& q)
  {return !(p <= q);}


template<class T> inline T& Polynomial<T>::operator[] (const Ulong& j)
  {return v[j];}
template<class T> inline void Polynomial<T>::reduceDeg() {v.reduceDim();}
template<class T> inline void Polynomial<T>::setDeg(const Degree& d)
  {v.setDim(d+1);}
template<class T> inline void Polynomial<T>::setDegValue(const Degree& d)
  {v.setDimValue(d+1);}
template<class T> inline void Polynomial<T>::setVect(const T *source,
						     const Ulong& n)
  {v.setVect(source,n);}
template<class T> inline void Polynomial<T>::setZero()
  {v.dim() = 0;}
template<class T> inline void Polynomial<T>::setZero(const Ulong& r)
  {v.setZero(r);}
template<class T> inline void Polynomial<T>::setZero(const Ulong& first,
						     const Ulong& r)
  {v.setZero(first,r);}
template<class T> inline vector::Vector<T>& Polynomial<T>::vect() {return v;}

template<class T> inline const T& Polynomial<T>::operator[] (const Ulong& j)
  const {return v[j];}
template<class T>
inline Polynomial<T>& Polynomial<T>::operator= (const Polynomial<T>& q)
  {v = q.v; return *this;}

template<class T> inline Ulong Polynomial<T>::deg() const {return v.dim()-1;}
template<class T> inline bool Polynomial<T>::isZero() const
  {return deg() == undef_degree;}
template<class T> inline const vector::Vector<T>& Polynomial<T>::vect() const
  {return v;}

}; // |namespace polynomials|

#include "polynomials.hpp"

#endif

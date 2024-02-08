/*
  This is polynomials.cpp

  Coxeter version 3.0 Copyright (C) 2002 Fokko du Cloux
  See file LICENSE for full copyright notice
*/

/****************************************************************************

  This module defines the polynomials used in this program.

  A polynomial is in fact just a vector, suitably reinterpreted. The
  degree of the polynomial is the (reduced) dimension of the vector
  minus one; we always make sure that the degree is exactly right,
  i.e., it is undef_degree = -1 for the zero polynomial, and otherwise we
  have v[deg] != 0.

 ****************************************************************************/

/****************************************************************************

        Chapter I --- The Polynomial class

 ****************************************************************************/

namespace polynomials {

template <class T> Polynomial<T>::~Polynomial()

/*
  Destruction is automatic.
*/

{}

/******** operators *********************************************************/

template <class T>
Polynomial<T>& Polynomial<T>::operator+= (const Polynomial<T>& q)

/*
  Adds q to the current polynomial.
*/

{
  /* check for zero */

  if (q.isZero())  /* do nothing */
    return *this;

  if (isZero()) { /* set p to q */
    *this = q;
    return *this;
  }

  v += q.v;
  v.reduceDim();  /* set the degree to the correct value */

  return *this;
}


template <class T>
Polynomial<T>& Polynomial<T>::operator-= (const Polynomial<T>& q)

/*
  Adds q to the current polynomial.
*/

{
  Degree j;

  /* check for zero */

  if (q.isZero())  /* do nothing */
    return *this;

  if (isZero()) { /* set p to -q */
    *this = q;
    for (j = 0; j <= deg(); j++)
      v[j] = -v[j];
    return *this;
  }

  v -= q.v;
  v.reduceDim();  /* set the degree to the correct value */

  return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator*= (const Polynomial<T>& q)

/*
  Multiplies the current polynomial by q.

  NOTE : we have used tmp as workspace for simplicity, although this can
  be eliminated altogether.
*/

{
  static vector::Vector<T> tmp;

  if ( isZero() || q.isZero() )  /* result is 0 */
    {
      setDeg(undef_degree);
      return *this;
    }

  tmp.setDim(deg()+q.deg()+1);

  Degree i, j;

  for (i = 0; i <= deg(); i++)
    for (j = 0; j <= q.deg(); j++)
      tmp[i+j] += v[i]*q[j];

  setDeg(deg()+q.deg()+1);
  v = tmp;

  return *this;
}


template <class T> Polynomial<T>& Polynomial<T>::operator*= (const T& a)

/*
  Scalar multiplication by a.
*/

{
  if ( isZero() || (a == 0) /* should be isZero<T>(a) */)  /* result is 0 */
    {
      setDeg(undef_degree);
      return *this;
    }

  v *= a;

  return *this;
}

template <class T>
Polynomial<T>& Polynomial<T>::operator/= (const Polynomial<T>& q)

/*
  Euclidian division by q; we assume that the leading coefficient of
  q is one. It turns out that everything can be done within p's
  memory space (even when p = q!).
*/

{
  if (deg() < q.deg()) { /* result is 0 */
    setDeg(undef_degree);
    return *this;
  }

  Degree i, j;

  for (j = deg()-q.deg()+1; j;)
    {
      j--;
      for (i = 0; i < q.deg(); i++)
        v[i+j] -= v[q.deg()+j]*q[i];
    }

  v.setVect(v.ptr()+q.deg(),deg()-q.deg()+1);
  setDeg(deg()-q.deg());

  return *this;
}

};

/*****************************************************************************

        Chapter II -- Comparison operators.

  This section defines comparison operators for polynomials, assuming
  that the corresponding comparison operators for type T have beeb defined.

  We define the following :

   - operator== (p,q);
   - operator!= (p,q); (inlined)
   - operator< (p,q);
   - operator> (p,q);

******************************************************************************/

namespace polynomials {

template <class T>
bool operator== (const Polynomial<T>& p, const Polynomial<T>& q)

/*
  Equality operator for polynomials. Two polynomials are equal, if both
  are zero, or if their degrees are equal and the corresponding coefficients
  are equal. It is assumed that operator != is defined for T.
*/

{
  if (p.isZero())
    return q.isZero();

  if (p.deg() != q.deg())
    return false;

  for (Ulong j = 0; j <= p.deg(); ++j) {
    if (p[j] != q[j])
      return false;
  }

  return true;
}

template <class T>
bool operator<= (const Polynomial<T>& p, const Polynomial<T>& q)

/*
  We have p <= q if either p = 0, or deg(p) < deg(q), or degrees are equal and
  coefficients are in order, starting from the top. It is assumed that
  operator< is defined for T
*/

{
  if (p.deg() < q.deg())
    return true;
  if (p.deg() > q.deg())
    return false;

  /* now degrees are equal */

  for (Ulong j = p.deg()+1; j;) { /* loop is empty if p = 0 */
    --j;
    if (p[j] > q[j])
      return false;
    if (p[j] < q[j])
      return true;
  }

  /* now polynomials are equal */

  return true;
}

template <class T>
bool operator>= (const Polynomial<T>& p, const Polynomial<T>& q)

/*
  We have p >= q if either deg(p) > deg(q), or degrees are equal and coefficients
  are in order, starting from the top. It is assumed that operator< is defined
  for T
*/

{
  if (p.deg() > q.deg())
    return true;
  if (p.deg() < q.deg())
    return false;

  /* now degrees are equal */

  for (Ulong j = p.deg()+1; j;) { /* loop is empty if p = 0 */
    --j;
    if (p[j] < q[j])
      return false;
    if (p[j] > q[j])
      return true;
  }

  /* now polynomials are equal */

  return true;
}

};

/**************************************************************************

        Chapter III --- Input/output

  This section defines i/o functions for polynomials :

   - append(str,p,x) : appends the representation of p to the string l (also
     defined for Laurent polynomials);
   - append(str,p,d,m,x) : same, substituting X^d and shifting degrees by m;
   - append(str,p,d,m,x,GAP) : same as the previous one, in GAP style;
   - append(str,p,d,m,x,Terse) : same as the previous one, in Terse style;
   - print(file,p,x) : appends the representation of p to the file (also
     defined for Laurent polynomials);
   - print(str,p,d,m,x) : same, substituting X^d and shifting degrees by m;
   - print(str,p,d,m,x,GAP) : same as the previous one, in GAP style;
   - print(str,p,d,m,x,Terse) : same as the previous one, in Terse style;

 **************************************************************************/

namespace polynomials {

template <class T>
std::string& append(std::string& str, const Polynomial<T> &p, const char *x)

/*
  Appends the string representation of p to l, using the string x to
  represent the indeterminate.
*/

{
  if (p.isZero()) {
    str.append("0");
    return str;
  }

  int firstcoeff = 1;
  Degree j = p.deg()+1;

  while (j) {
    j--;
    if (p[j] == 0)
      continue;
    if (firstcoeff)
      firstcoeff = 0;
    else
      if (p[j] > 0)
	str.append("+");
    switch (j)
      {
      case 0:
	io::append(str,p[j]);
	break;
      default:
	if ((p[j] != 1) && (p[j] != (T)(-1)))
	  io::append(str,p[j]);
	else if (p[j] == (T)(-1))
	  str.append("-");
	break;
      };
    switch (j)
      {
      case 0:
	break;
      case 1:
	str.append(x);
	break;
      default:
	str.append(x);
	str.append("^");
	io::append(str,j);
	break;
      };
  }

  return str;
}

template <class T>
std::string& append(std::string& str, const LaurentPolynomial<T> &p, const char *x)

/*
  Appends the string representation of p to l, using the string x to
  represent the indeterminate.
*/

{
  if (p.isZero()) {
    str.append("0");
    return str;
  }

  int firstcoeff = 1;

  for (long j = p.val(); j <= p.deg(); ++j) {
    if (p[j] == 0)
      continue;
    if (firstcoeff)
      firstcoeff = 0;
    else
      if (p[j] > 0)
	str.append("+");
    switch (j)
      {
      case 0:
	io::append(str,p[j]);
	break;
      default:
	if ((p[j] != 1) && (p[j] != (T)(-1)))
	  io::append(str,p[j]);
	else if (p[j] == (T)(-1))
	  str.append("-");
	break;
      };
    switch (j)
      {
      case 0:
	break;
      case 1:
	str.append(x);
	break;
      default:
	str.append(x);
	str.append("^");
	io::append(str,j);
	break;
      };
   }

  return str;
}

template <class T>
std::string& append(std::string& str, const Polynomial<T> &p, const Degree& d,
	       const long& m, const char *x)

/*
  Appends the string representation of p to str, using the string x to
  represent the indeterminate. In this version, X^d is first substituted
  in the polynomial, and afterwards the whole thing is shifted by m.
*/

{
  if (p.deg() == undef_degree) {
    str.append("0");
    return str;
  }

  int firstcoeff = 1;
  Degree j = p.deg()+1;

  while (j) {
    j--;
    if (p[j] == 0)
      continue;
    if (firstcoeff)
      firstcoeff = 0;
    else
      if (p[j] > 0)
	str.append("+");
    long a = j*d + m;
    switch (a) {
    case 0:
      append(str,p[j]);
      break;
    default:
      if ((p[j] != (T)1) && (p[j] != (T)(-1)))
	append(str,p[j]);
      else if (p[j] == (T)(-1))
	str.append("-");
      break;
    };
    switch (a) {
    case 0:
      break;
    case 1:
      str.append(x);
      break;
    default:
      str.append(x);
      str.append("^");
      io::append(str,a);
      break;
    };
  }

  return str;
}

template <class T>
std::string& append(std::string& str, const Polynomial<T> &p, const Degree& d,
	       const long& m, const char *x, io::GAP)

/*
  Appends the GAP representation of p to str, using the string x to
  represent the indeterminate. In this version, X^d is first substituted
  in the polynomial, and afterwards the whole thing is shifted by m.

  The only difference with the ordinary print is that a * is required between
  then coefficient and the indeterminate.
*/

{
  if (p.deg() == undef_degree) {
    str.append("0");
    return str;
  }

  int firstcoeff = 1;
  Degree j = p.deg()+1;

  while (j) {
    j--;
    if (p[j] == 0)
      continue;
    if (firstcoeff)
      firstcoeff = 0;
    else
      if (p[j] > 0)
	str.append("+");
    long a = j*d + m;
    switch (a) {
    case 0:
      append(str,p[j]);
      break;
    default:
      if ((p[j] != (T)1) && (p[j] != (T)(-1))) {
	append(str,p[j]);
	str.append("*");
      }
      else if (p[j] == (T)(-1))
	str.append("-");
      break;
    };
    switch (a) {
    case 0:
      break;
    case 1:
      str.append(x);
      break;
    default:
      str.append(x);
      str.append("^");
      io::append(str,a);
      break;
    };
  }

  return str;
}

template <class T>
std::string& append(std::string& str, const Polynomial<T> &p, const Degree& d,
	       const long& m, const char *x, io::Terse)

/*
  Appends the Terse representation of p to str, using the string x to
  represent the indeterminate. In this version, X^d is first substituted
  in the polynomial, and afterwards the whole thing is shifted by m.

  In terse style, polynomials are represented as comma-separated lists
  of coefficients, enclosed by parentheses. Whenever either d is different
  from one, or m is different from zero, the polynomial is preceded by
  a pair (d,m), also enclosed in parentheses.
*/

{
  if (p.deg() == undef_degree) {
    str.append("()");
    return str;
  }

  if ((d != 1) || (m != 0)) {
    str.append("(");
    io::append(str,d);
    str.append(",");
    io::append(str,m);
    str.append(")");
  }

  str.append("(");

  for (Ulong j = 0; j <= p.deg(); ++j) {
    str.append(p[j]);
    if ((j+1) <= p.deg()) /* there is more to come */
      str.append(",");
  }

  str.append(")");

  return str;
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const char* x)

{
  static std::string buf;

  buf.clear();
  append(buf,p,x);
  io::print(file,buf);

  return;
}

template <class T> void print(FILE* file, const LaurentPolynomial<T>& p,
			      const char* x)

{
  static std::string buf;

  buf.clear();
  append(buf,p,x);
  io::print(file,buf);

  return;
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const Degree& d, const long& m, const char* x)

/*
  Prints the polynomial with x^d substituted for x, and shifted by m (so that
  we may actually be printing a Laurent polynomial.)
*/

{
  static std::string buf;

  buf.clear();
  append(buf,p,d,m,x);
  io::print(file,buf);

  return;
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const Degree& d, const long& m, const char* x,
			      io::GAP)

/*
  Prints the polynomial with x^d substituted for x, and shifted by m (so that
  we may actually be printing a Laurent polynomial.)
*/

{
  static std::string buf;

  buf.clear();
  append(buf,p,d,m,x,io::GAP());
  io::print(file,buf);

  return;
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const Degree& d, const long& m, const char* x,
			      io::Terse)

/*
  Prints the polynomial in terse style.
*/

{
  static std::string buf;

  buf.clear();
  append(buf,p,d,m,x,io::Terse());
  io::print(file,buf);

  return;
}

/*****************************************************************************

        Chapter IV -- The LaurentPolynomial class

  The LaurentPolynomial class is just an adaptor for a polynomial; all
  operations are done at the level of an underlying polynomial, and there
  is just a shift to be taken into account. The degree (valuation) of a Laurent
  polynomial is the position of its largest (smallest) non-zero coefficient.

  The following functions are defined :

   - constructors and destructors :

     - LaurentPolynomial(const SDegree&, const SDegree&) : constructs a
       Laurent polynomial capable of holding the given degree and valuation;
     - ~LaurentPolynomial();

   - accessors :

     - operator== (const LaurentPolynomial<T>&);
     - operator<= (const LaurentPolynomial<T>&);
     - operator>= (const LaurentPolynomial<T>&);

   - manipulators :

     - setBounds(const SDegree&, const SDegree&) : sets both the degree and
       the valuation, enlarging if necessary;
     - setDeg(const SDegree&) : sets the degree, enlarging if necessary;
     - setVal(const SDegree&) : sets the valuation, enlarging if necessary;

  The following comparison operators are defined :

   - operator== : equality of valuations and polynomials;
   - operator!= : the negation of ==; (inlined)
   - operator<= : inequality of valuations and of polynomials; (inlined)
   - operator>= : inequality of valuations and polynomials; (inlined)
   - operator< : negation of <=; (inlined)
   - operator> : negation of >=; (inlined)

 *****************************************************************************/


/******** accessors *********************************************************/

template<class T>
bool LaurentPolynomial<T>::operator== (const LaurentPolynomial<T>& p) const
{
  if (isZero())
    return p.isZero();
  if (p.isZero())
    return false;

  if (d_valuation != p.d_valuation)
    return false;

  return coef == p.coef;
}


/*
  Comparison operator for Laurent polynomials. Zero is larger than any
  polynomial; otherwise comparison is valuation-first.
*/
template<class T>
bool LaurentPolynomial<T>::operator<= (const LaurentPolynomial<T>& p) const
{
  if (p.isZero())
    return true;
  if (isZero())
    return false;
  if (d_valuation < p.d_valuation)
    return true;

  return  coef <= p.coef;
}

template<class T>
bool LaurentPolynomial<T>::operator>= (const LaurentPolynomial<T>& p) const

/*
  Comparison operator for Laurent polynomials. Zero is larger than any
  polynomial; otherwise comparison is valuation-first.
*/

{
  if (isZero())
    return true;
  if (p.isZero())
    return false;
  if (d_valuation > p.d_valuation)
    return true;

  return  coef >= p.coef;
}

/******** manipulators ******************************************************/


};

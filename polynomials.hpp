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


/******** operators *********************************************************/

template <class T>
Polynomial<T>& Polynomial<T>::operator+= (const Polynomial<T>& q)
{
  if (q.isZero())  /* do nothing */
    return *this;

  if (isZero())
    return *this = q;

  if (deg()>q.deg())
    for (unsigned i=0; i<q.v.size(); ++i)
      v[i]+=q.v[i];
  else if (deg()<q.deg())
  {
    for (unsigned i=0; i<v.size(); ++i)
      v[i]+=q.v[i];
    v.insert(v.end(),q.v.begin()+v.size(),q.v.end()); // copy high coefficients
  }
  else
  {
    for (unsigned i=0; i<v.size(); ++i)
      v[i]+=q.v[i];
    snap_degree();  /* set the degree to the correct value */
  }
  return *this;
}


template <class T>
Polynomial<T>& Polynomial<T>::operator-= (const Polynomial<T>& q)
{
  if (q.isZero())  /* do nothing */
    return *this;

  if (isZero())
  {
    v = q.v;
    for (auto& x : v)
      x = -x;
    return *this;
  }

  if (deg()>q.deg())
    for (unsigned i=0; i<q.v.size(); ++i)
      v[i]-=q.v[i];
  else if (deg()<q.deg())
  {
    v.resize(q.v.size(),0);
    for (unsigned i=0; i<v.size(); ++i)
      v[i]-=q.v[i];
  }
  else
  {
    for (unsigned i=0; i<v.size(); ++i)
      v[i]-=q.v[i];
    snap_degree();  /* set the degree to the correct value */
  }
  return *this;
}

template <class T> Polynomial<T>& Polynomial<T>::operator*= (const T& a)
{
  if (a == T(0))
  {
    v.clear();
    return *this;
  }

  for (auto& x : v)
    x *= a;

  return *this;
}


/*
  Multiply our polynomial by |q|.

  NOTE : Fokko noted that (supposedly after resizing) all can be done in-place
  in the vector |v|, which is therefore what we implemented below.
*/
template <class T>
Polynomial<T>& Polynomial<T>::operator*= (const Polynomial<T>& q)
{
  if (isZero())
    return *this;
  if (q.isZero())
    v.clear();
  else
  {
    v.resize(v.size()+q.deg(),0);
    for (Degree i=v.size(); i-->0; )
    {
      auto c = v[i]; // save coefficient that will be overwritten
      Degree j=i;
      auto it=q.v.begin();
      for (v[j++] = c * *it++; it!=q.v.end(); ++j,++it)
	v[j] += c * *it;
    }
  }
  return *this;
}




/*
  Euclidian division by |q|; we assume that the leading coefficient of
  |q| is one. It turns out that everything can be done within p's
  memory space (even when p = q!).
*/
template <class T>
Polynomial<T>& Polynomial<T>::operator/= (const Polynomial<T>& q)
{
  assert(q.v.back()==T(1)); // we only implement division by a monic polynomial
  if (deg() < q.deg()) { // quotient is 0, remainder |*this| is thrown away
    v.clear();
    return *this;
  }

  const Degree d=q.deg();
  for (Degree j = v.size(); j-->d;)
  {
    auto it = &v[j];
    auto c = *it; // this leading coefficient of remainder now enters quotient
    for (Degree i = q; i-->0; ) // use all non-leading coefs of |q|
      *--it -= c * q[i]; // subtract product from tommost coefs of remainder
  }

  v.erase(&v[0],&v[d]); // drop size |d| remainder, and shift quotient in place

  return *this;
}


/*****************************************************************************

        Chapter II -- Comparison operators. Everything is inlined

******************************************************************************/

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

  bool firstcoeff = true;

  for (Degree j = p.deg()+1; j-->0; )
  {
    if (p[j] == 0)
      continue; // skip zero terms
    if (firstcoeff) // then skip printing "+"
      firstcoeff = false;
    else if (p[j] > 0)
      str.append("+");

    if (j==0)
      io::append(str,p[j]); // print the constant term, even if $\pm1$
    else if (p[j] != T(1) and p[j] != T(-1)) // deg>0, non-unit coefficient
      io::append(str,p[j]); // print it with possible minus sign
    else if (p[j] == T(-1)) // for coef $-1$ just print '-' before |x|
      str.append("-");

    switch (j) // print monomial in |x|
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
    } // |switch|
  }

  return str;
}


/*
  Appends the string representation of p to l, using the string x to
  represent the indeterminate.
*/
template <class T>
  std::string& append
    (std::string& str, const LaurentPolynomial<T> &p, const char *x)
{
  if (p.isZero()) {
    str.append("0");
    return str;
  }

  bool firstcoeff = true;

  for (long j = p.val(); j <= p.deg(); ++j) // print Laurent polynomials rising
  {
    if (p[j] == 0)
      continue;
    if (firstcoeff)
      firstcoeff = false;
    else if (p[j] > 0)
      str.append("+");

    if (j==0)
      io::append(str,p[j]); // print the constant term, even if $\pm1$
    else if (p[j] != T(1) and p[j] != T(-1)) // deg>0, non-unit coefficient
      io::append(str,p[j]); // print it with possible minus sign
    else if (p[j] == T(-1)) // for coef $-1$ just print '-' before |x|
      str.append("-");

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


/*
  Appends the string representation of p to str, using the string x to
  represent the indeterminate. In this version, the substitution $X:=X^d$
  is first made, followed by a multiplication by $X^m$.
*/
template <class T>
std::string& append
  (std::string& str, const Polynomial<T> &p,
   const Degree& d, const long& m, const char *x)
{
  if (p.isZero()) {
    str.append("0");
    return str;
  }

  bool firstcoeff = true;

  for (Degree j = p.deg()+1; j-->0; )
  {
    if (p[j] == 0)
      continue;
    if (firstcoeff)
      firstcoeff = false;
    else if (p[j] > 0)
      str.append("+");

    auto a = j*d + m; // effective degree
    if (a==0)
      io::append(str,p[j]); // print the constant term, even if $\pm1$
    else if (p[j] != T(1) and p[j] != T(-1)) // deg>0, non-unit coefficient
      io::append(str,p[j]); // print it with possible minus sign
    else if (p[j] == T(-1)) // for coef $-1$ just print '-' before |x|
      str.append("-");

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


/*
  Append the GAP representation of p to str, using the string x to
  represent the indeterminate. In this version, X^d is first substituted
  in the polynomial, and afterwards the whole thing is shifted by m.

  The only difference with the ordinary print is that a * is required between
  then coefficient and the indeterminate.
*/
template <class T>
std::string& append(std::string& str, const Polynomial<T> &p, const Degree& d,
	       const long& m, const char *x, io::GAP)
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


/*
  Appends the Terse representation of p to str, using the string x to
  represent the indeterminate. In this version, X^d is first substituted
  in the polynomial, and afterwards the whole thing is shifted by m.

  In terse style, polynomials are represented as comma-separated lists
  of coefficients, enclosed by parentheses. Whenever either d is different
  from one, or m is different from zero, the polynomial is preceded by
  a pair (d,m), also enclosed in parentheses.
*/
template <class T>
std::string& append(std::string& str, const Polynomial<T> &p, const Degree& d,
	       const long& m, const char *x, io::Terse)
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
  std::string buf;
  append(buf,p,x);
  io::print(file,buf);
}

template <class T> void print(FILE* file, const LaurentPolynomial<T>& p,
			      const char* x)
{
  std::string buf;
  append(buf,p,x);
  io::print(file,buf);
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const Degree& d, const long& m, const char* x)
{
  std::string buf;
  append(buf,p,d,m,x);
  io::print(file,buf);
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const Degree& d, const long& m, const char* x,
			      io::GAP)
{
  std::string buf;
  append(buf,p,d,m,x,io::GAP());
  io::print(file,buf);
}

template <class T> void print(FILE* file, const Polynomial<T>& p,
			      const Degree& d, const long& m, const char* x,
			      io::Terse)
{
  std::string buf;
  append(buf,p,d,m,x,io::Terse());
  io::print(file,buf);
}

/*****************************************************************************

        Chapter IV -- The LaurentPolynomial class

  The LaurentPolynomial class is just an adaptor for a polynomial; all
  operations are done at the level of an underlying polynomial, and there
  is just a shift to be taken into account. The degree (valuation) of a Laurent
  polynomial is the position of its largest (smallest) non-zero coefficient.


*/

// comparison is there just for the sake of ordering containers

template<class T>
bool LaurentPolynomial<T>::operator< (const LaurentPolynomial<T>& p) const
{ if (isZero())
    return false; // zero is maximum
  else if (p.isZero())
    return true;
  else if (val()!=p.val())
    return p.val()<val(); // compare opposite to valuation
  return coef<p.coef;
}

};
